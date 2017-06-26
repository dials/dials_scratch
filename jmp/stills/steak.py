from __future__ import division


def golden_section_search(f, a, b, tol=1e-3):
  from math import sqrt
  gr = (sqrt(5) + 1) / 2
  c = b - (b - a) / gr
  d = a + (b - a) / gr
  while abs(c - d) > tol:
    print a, b
    if f(c) < f(d):
      b = d
    else:
      a = c

    c = b - (b - a) / gr
    d = a + (b - a) / gr

  return (b + a) / 2

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
                                    tolerance=1e-10,
                                    max_iter=max_iter)

    self.x = optimizer.get_solution()
    return

  def get_solution(self):
    return self.x

  def target(self, vector):
    print "SCORE"
    score = scorify(vector)
    return score


class CrystalRefiner(object):

  def __init__(self,
               experiments,
               reflections,
               mosaicity):

    self.experiments = experiments
    self.reflections = reflections
    self.mosaicity = mosaicity

  def refine(self):
    from dials.algorithms.refinement.parameterisation.crystal_parameters \
      import CrystalUnitCellParameterisation, \
      CrystalOrientationParameterisation
    from dials.array_family import flex
    from scitbx import simplex

    crystal = self.experiments[0].crystal

    self.cucp = CrystalUnitCellParameterisation(crystal)
    self.cop = CrystalOrientationParameterisation(crystal)

    # 0-point and deltas
    values = flex.double(self.cucp.get_param_vals() +
                         self.cop.get_param_vals())
    offset = flex.double([0.01 * v for v in self.cucp.get_param_vals()] +
                         [0.1, 0.1, 0.1])

    initial = crystal.get_unit_cell()
    self.cells = []
    initial_score = self.target(values)
    doohicky = simple_simplex(values, offset, self, 2000)
    best = doohicky.get_solution()
    print 'Initial cell:', initial
    print 'Final cell:  ', crystal.get_unit_cell()
    print 'Score change', initial_score, self.target(best)
    self.best = best

  def target(self, vector):
    from dials_scratch.jmp.stills import Model
    from dials.array_family import flex
    cell_parms = self.cucp.get_param_vals()
    orientation_parms = self.cop.get_param_vals()
    assert len(vector) == len(cell_parms) + len(orientation_parms)
    tst_cell = vector[:len(cell_parms)]
    tst_orientation = vector[len(cell_parms):len(cell_parms) +
                             len(orientation_parms)]

    self.cucp.set_param_vals(tst_cell)
    self.cop.set_param_vals(tst_orientation)

    from scitbx import matrix
    experiment = self.experiments[0]

    print 'Cell: %.3f %.3f %.3f %.3f %.3f %.3f' % \
      tuple(experiment.crystal.get_unit_cell().parameters())
    print 'Phi(1,2,3): %.3f %.3f %.3f' % tuple(tst_orientation)



    data = experiment.imageset.get_raw_data(0)[0]
    mask = experiment.imageset.get_mask(0)[0]

    model = Model(
      experiment.beam,
      experiment.detector,
      experiment.crystal,
      self.mosaicity,
      5,
      data,
      mask,
      self.reflections)


    score = 0

    for i in range(len(model)):

      B = model.background(i)
      I = model.intensity(i)
      V = model.variance(i)
      S = model.scale(i)
      T = model.success(i)
      XYZOBS = model.observed(i)
      XYZCAL = model.predicted(i)

      if T == True:
        score += (matrix.col(XYZOBS) - matrix.col(XYZCAL)).length_sq()

    return score


class ProfileRefiner(object):

  def __init__(self,
               experiments,
               bandwidth=0.035,
               mosaicity=0.01):

    self.experiments = experiments
    self.bandwidth = bandwidth
    self.mosaicity = mosaicity

  def find_spots(self, min_spot_size=2, max_spot_size=100):
    from dials.algorithms.spot_finding.threshold import XDSThresholdStrategy
    from dials.model.data import PixelList
    from dials.model.data import PixelListLabeller
    from dials.array_family import flex
    from scitbx import matrix

    image = self.experiments[0].imageset.get_raw_data(0)[0]
    mask = self.experiments[0].imageset.get_mask(0)[0]

    threshold_image = XDSThresholdStrategy()

    threshold_mask = threshold_image(image, mask)
    plist = PixelList(0, image, threshold_mask)

    pixel_labeller = PixelListLabeller()
    pixel_labeller.add(plist)

    creator = flex.PixelListShoeboxCreator(
      pixel_labeller, 0, 0, True, min_spot_size, max_spot_size, False)
    shoeboxes = creator.result()

    # turns out we need to manually filter the list to get a sensible answer
    size = creator.spot_size()
    big = size > max_spot_size
    small = size < min_spot_size
    bad = big | small
    shoeboxes = shoeboxes.select(~bad)


    centroid = shoeboxes.centroid_valid()
    intensity = shoeboxes.summed_intensity()
    observed = flex.observation(shoeboxes.panels(), centroid, intensity)

    self.reflections = flex.reflection_table(observed, shoeboxes)

    miller_index = flex.miller_index()
    UB = matrix.sqr(self.experiments[0].crystal.get_A())
    UBi = UB.inverse()

    self.qobs = []

    s0 = matrix.col(self.experiments[0].beam.get_s0())
    selection = flex.size_t()
    for index, refl in enumerate(self.reflections):
      x, y, z = refl['xyzobs.px.value']
      p = matrix.col(self.experiments[0].detector[0].get_pixel_lab_coord((x, y)))
      q = p.normalize() * s0.length() - s0
      self.qobs.append(q)
      hkl = UBi * q
      ihkl = [int(round(h)) for h in hkl]

      if (ihkl[0]-hkl[0])**2 + (ihkl[0]-hkl[0])**2 + (ihkl[0]-hkl[0])**2 < 0.3:
        selection.append(index)

      miller_index.append(ihkl)


    self.reflections['miller_index'] = miller_index

    print len(selection), len(self.reflections)

    self.reflections = self.reflections.select(selection)

    # from matplotlib import pylab
    # X, Y, _ = self.reflections['xyzobs.px.value'].parts()
    # pylab.scatter(X,Y)
    # pylab.show()

  def refine(self):

    self.find_spots()


    #self.sigma_m = golden_section_search(self.target, 0.01, 0.2, tol=1e-5)
    self.sigma_m = 0.0674634452691

    print self.sigma_m

    refiner = CrystalRefiner(self.experiments, self.reflections, self.sigma_m)
    experiments = refiner.refine()

    # import numpy as np
    # X = []
    # Y = []
    # for m in np.linspace(0.030, 0.2, 50):
    #   s = self.target(m)
    #   print m, s
    #   X.append(m)
    #   Y.append(s)

    # from matplotlib import pylab
    # pylab.plot(X, Y)
    # pylab.show()

  def target(self, mosaicity):
    from dials_scratch.jmp.stills import Model
    from dials.array_family import flex

    experiment = self.experiments[0]

    data = experiment.imageset.get_raw_data(0)[0]
    mask = experiment.imageset.get_mask(0)[0]

    model = Model(
      experiment.beam,
      experiment.detector,
      experiment.crystal,
      mosaicity,
      5,
      data,
      mask,
      self.reflections)

    #print flex.max(model.image_data())

    # from matplotlib import pylab
    # pylab.imshow(model.image_data().as_numpy_array(), interpolation='none')
    # pylab.show()
    # pylab.imshow(model.image_mask().as_numpy_array(), interpolation='none')
    # pylab.show()
    # pylab.imshow(model.image_pred().as_numpy_array(), interpolation='none')
    # pylab.show()

    if False:
      score = self.least_squares_score(model)
    else:
      score = self.maximum_likelihood_score(model)

    return score

  def least_squares_score(self, model):
    from math import sqrt

    score = 0

    for i in range(len(model)):

      B = model.background(i)
      I = model.intensity(i)
      V = model.variance(i)
      S = model.scale(i)
      T = model.success(i)

      if T == True and I > 0 and S > 0:# and S > 1e-10 and V > 0 and I / sqrt(V) > 2:

        P = I / S
        c = model.data(i)
        m = model.mask(i)
        p = model.pred(i)

        score += sum((ci - (B + P*pi))**2 for ci, pi, mi in zip(c, p, m) if m & 5)

    return score

  def maximum_likelihood_score(self, model):
    from math import log

    score = 0

    for i in range(len(model)):

      B = model.background(i)
      I = model.intensity(i)
      V = model.variance(i)
      S = model.scale(i)
      T = model.success(i)

      if T == True and I > 0 and S > 0:# and S > 0 and V > 0:

        P = I / S
        c = model.data(i)
        m = model.mask(i)
        p = model.pred(i)

        score += sum((ci*log(B+P*pi) - B - P*pi)  for ci, pi, mi in zip(c, p, m) if m & 5)

    return -score


if __name__ == '__main__':

  from dxtbx.model.experiment_list import ExperimentListFactory
  import sys

  experiments_filename = sys.argv[1]

  # Read the experiments
  experiments = ExperimentListFactory.from_json_file(experiments_filename)

  # Setup the profile refiner
  refiner = ProfileRefiner(
    experiments,
    bandwidth=0.035,
    mosaicity=0.1)

  # Do the refinement
  refiner.refine()
