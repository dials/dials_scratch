from __future__ import division
from __future__ import print_function
from libtbx.phil import parse
from math import sqrt, exp, pi, acos
from scitbx import matrix
from scitbx.linalg import eigensystem

def generate_start(values, offset):
  assert len(values) == len(offset)
  start = [values]
  for j, o in enumerate(offset):
    next = values.deep_copy()
    next[j] += o
    start.append(next)
  return start

class simple_simplex(object):
  def __init__(self, values, offset, evaluator, max_iter):
    from scitbx import simplex
    self.n = len(values)
    self.x = values
    self.starting_simplex = generate_start(values, offset)
    self.fcount = 0

    optimizer = simplex.simplex_opt(dimension=self.n,
                                    matrix=self.starting_simplex,
                                    evaluator=evaluator,
                                    tolerance=1e-7,
                                    max_iter=max_iter)

    self.x = optimizer.get_solution()
    return

  def get_solution(self):
    return self.x

  def target(self, vector):
    print("SCORE")
    score = scorify(vector)
    return score


def compute_mean_plane(mu, sigma, s0_length):
  z = matrix.col((0, 0, 1))
  axis = z.cross(mu.normalize())
  angle = acos(z.dot(mu.normalize()))
  R = axis.axis_and_angle_as_r3_rotation_matrix(angle=angle, deg=False)
  sigma_1 = R.transpose() * sigma * R
  mu_1 = R.transpose() * mu
  sigma_11 = matrix.sqr((
    sigma_1[0], sigma_1[1],
    sigma_1[3], sigma_1[4]))
  sigma_12 = matrix.col((sigma_1[2], sigma_1[5]))
  sigma_21 = matrix.col((sigma_1[3], sigma_1[7])).transpose()
  sigma_22 = sigma_1[8]
  mu1 = matrix.col((mu_1[0], mu_1[1]))
  mu2 = mu_1[2]
  mu_new_1 = mu1 + sigma_12 * (1/sigma_22) * (s0_length - mu2)
  v = matrix.col((
    mu_new_1[0],
    mu_new_1[1],
    s0_length)).normalize() * s0_length
  x_new = R * v
  return x_new

def generate_observations(experiments, reflections, sigma):

  A = matrix.sqr(experiments[0].crystal.get_A())
  s0 = matrix.col(experiments[0].beam.get_s0())

  s1_obs = flex.vec3_double()
  s2_obs = flex.vec3_double()
  for i in range(len(reflections)):

    h = matrix.col(reflections[i]['miller_index'])

    r = A * h
    s2 = s0 + r

    s1 = compute_mean_plane(s2, sigma, s0.length())

    s1_obs.append(s1)
    s2_obs.append(s2)

  return s1_obs, s2_obs


def generate_observations2(experiments, reflections, sigma):
  from dials_scratch.jmp.stills.potato import PotatoOnEwaldSphere
  from numpy.random import poisson

  A = matrix.sqr(experiments[0].crystal.get_A())
  s0 = matrix.col(experiments[0].beam.get_s0())


  s1_obs = flex.vec3_double()
  s2_obs = flex.vec3_double()
  I_obs = flex.double()
  for i in range(len(reflections)):

    h = matrix.col(reflections[i]['miller_index'])

    r = A * h
    s2 = s0 + r

    model = PotatoOnEwaldSphere(1/s0.length(), s2, sigma)
    s1 = model.centre_of_mass_on_ewald_sphere()
    scale = model.scale_factor()

    I_obs.append(poisson(scale * 1000))
    s1_obs.append(s1)
    s2_obs.append(s2)

  return s1_obs, s2_obs, I_obs


class ProfileRefiner(object):

  def __init__(self,
               experiment,
               reflections):
    from dials.algorithms.refinement.parameterisation.crystal_parameters \
      import CrystalUnitCellParameterisation
    from dials.algorithms.refinement.parameterisation.crystal_parameters \
      import CrystalOrientationParameterisation
    from dials_scratch.jmp.stills import Model
    from dials.array_family import flex
    from scitbx import simplex
    from math import sqrt

    # Store the input
    self.experiment = experiment
    self.reflections = reflections

    self.crystal = self.experiment.crystal

    # Get the current values and generate some offsets
    values = flex.double((
      0.001,
      0, 0.001,
      0, 0, 0.001))
    offset = flex.double(
      [0.001  for v in values])

    # The optimization history
    self.history = []

    # Get the initial cell and initial score
    M = matrix.sqr((
      values[0], 0, 0,
      values[1], values[2], 0,
      values[3], values[4], values[5]))
    initial_sigma = M*M.transpose()
    initial_score = self.target(values)

    self.sigma = initial_sigma

    # Perform the optimization
    optimizer = simple_simplex(values, offset, self, 2000)
    result = optimizer.get_solution()
    M = matrix.sqr((
      result[0], 0, 0,
      result[1], result[2], 0,
      result[3], result[4], result[5]))
    self.sigma = M*M.transpose()
    print('Initial sigma:', initial_sigma)
    print('Final sigma:  ', self.sigma)

    # Compute the eigen decomposition of the covariance matrix
    eigen_decomposition = eigensystem.real_symmetric(self.sigma.as_flex_double_matrix())
    Q = matrix.sqr(eigen_decomposition.vectors())
    L = matrix.diag(eigen_decomposition.values())
    print(Q)
    print(L)

    # Compute RMSD
    # xcal, ycal, _ = self.model.observed().parts()
    # xobs, yobs, _ = self.model.predicted().parts()
    # rmsd_x = sqrt(flex.sum((xcal-xobs)**2) / len(xcal))
    # rmsd_y = sqrt(flex.sum((ycal-yobs)**2) / len(ycal))
    # print 'RMSD X, Y (px): %f, %f' % (rmsd_x, rmsd_y)

  def target(self, vector):
    '''
    The target function

    '''
    from dials_scratch.jmp.stills.potato import PotatoOnEwaldSphere
    from dials.array_family import flex
    from math import sqrt, log

    M = matrix.sqr((
      vector[0], 0, 0,
      vector[1], vector[2], 0,
      vector[3], vector[4], vector[5]))
    sigma = M*M.transpose()

    # Generate observed positions
    s1_cal, s2_cal, I_cal = generate_observations2([self.experiment], self.reflections, sigma)
    s0 = matrix.col((self.experiment.beam.get_s0()))
    I_obs = self.reflections['intensity.sum.value']
    lnL = 0
    for s2, I in zip(s2_cal, I_obs):
      model = PotatoOnEwaldSphere(1/s0.length(), s2, sigma)
      lnL += model.log_likelihood()# * I

    # Compute the rmsd between observed and calculated
    score = -lnL

    # Append to the history
    self.history.append((sigma, score))

    # Print some info
    print('Sigma: %.7f %.7f %.7f %.7f %.7f %.7f %.7f %.7f %.7f; Score: %f' % (
      tuple(sigma) + tuple((score,))))
    return score



phil_scope = parse('')


if __name__ == '__main__':

  from dials.array_family import flex
  from dials.util.options import OptionParser
  from dials.util.options import flatten_experiments
  import libtbx.load_env
  import sys

  # The script usage
  usage = "usage: %s [options] [param.phil] "\
          "{datablock.json | image1.file [image2.file ...]}" \
          % libtbx.env.dispatcher_name

  # Initialise the base class
  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    read_experiments=True)

  # Parse the command line
  params, options = parser.parse_args(show_diff_phil=False)

  # Ensure we have a data block
  experiments = flatten_experiments(params.input.experiments)
  experiments[0].scan.set_oscillation((0, 1), deg=True)

  # The predicted reflections
  reflections = flex.reflection_table.from_predictions_multi(experiments, padding=1)

  # Select only those within 1 deg
  x, y, z = reflections['xyzcal.px'].parts()
  selection = flex.abs(z) < 1
  reflections = reflections.select(selection)


  sigma = matrix.sqr((0.0001, 0, 0,
                      0, 0.0002, 0,
                      0, 0, 0.0003))

  # Generate observed positions
  s1_obs, s2_obs, I_obs = generate_observations2(experiments, reflections, sigma)

  selection = I_obs > 0
  s1_obs = s1_obs.select(selection)
  s2_obs = s2_obs.select(selection)
  I_obs = I_obs.select(selection)
  reflections = reflections.select(selection)

  print(len(selection), len(reflections))

  from matplotlib import pylab
  pylab.hist(I_obs, bins=50)
  pylab.show()


  # from matplotlib import pylab
  # s0 = matrix.col(experiments[0].beam.get_s0())
  # D = [matrix.col(s2).length()/s0.length() for s2 in s2_obs]
  # pylab.hist(D, bins=50)
  # pylab.show()

  angles = []
  for s1, s2 in zip(s1_obs, s2_obs):
    a = matrix.col(s1).angle(matrix.col(s2), deg=True)
    angles.append(a)
  print("Mean angle between s1 and s2 %f degrees " % (sum(angles) / len(angles)))

  # Do the ray intersection
  reflections['intensity.sum.value'] = I_obs
  reflections['s1_obs'] = s1_obs
  reflections['s1'] = s2_obs
  reflections['xyzobs.px'] = flex.vec2_double([experiments[0].detector[0].get_ray_intersection_px(s1) for s1 in s1_obs])


  # Do the refinement
  refiner = ProfileRefiner(experiments[0], reflections)
