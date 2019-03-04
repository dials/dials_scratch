from __future__ import division
from __future__ import print_function
from libtbx.phil import parse

def vonmises_fisher(mu, kappa):
  '''
  Sample from the von-mises fisher distribution using the algorithm of Wood
  1994.

  '''
  from numpy.random import beta, uniform
  from math import sqrt, log, acos, pi
  from scitbx import matrix
  p = len(mu)
  assert p == 3

  b = (-2*kappa + sqrt(4*kappa**2 + (p-1)**2)) / (p-1)
  x0 = (1-b)/(1+b)
  c = kappa * x0 + (p-1)*log(1-x0**2)
  while True:
    z = beta((p-1)/2.0,(p-1)/2.0)
    u = uniform(0,1)
    w = (1 - (1+b)*z) / (1 - (1-b)*z)
    if kappa*w + (p-1)*log(1+x0*w) -c >= log(u):
      break

  # Construct a uniform unit vector on a p-1 sphere
  v = matrix.col((
    uniform(-1, 1),
    uniform(-1, 1))).normalize()

  y = sqrt(1-w**2)

  # Construct the vmF distributed vector about (0, 0, 1)
  x = matrix.col((y*v[0], y*v[1], w))

  # Rotate the vmF distributed vector around mu
  a = matrix.col((0, 0, 1))
  b = mu.normalize()

  axis = a.cross(b)
  angle = acos(a.dot(b))

  if angle < 1e-7:
    return x
  elif abs(angle - pi) < 1e-7:
    return -x

  return x.rotate_around_origin(axis=axis, angle=angle) * mu.length()



def simulate_reflection(experiment,
                        reflection,
                        wavelength_spread = 0.0,
                        rlp_covariance    = None,
                        angular_spread    = 0.0,
                        beam_divergence   = 0.0):
  '''
  Simulate a count from a single reflection

  '''
  from scitbx import matrix
  from dials.algorithms.spot_prediction import ScanStaticRayPredictor
  from numpy.random import normal, multivariate_normal

  # Get the models from the experiment
  beam = experiment.beam
  detector = experiment.detector
  goniometer = experiment.goniometer
  scan = experiment.scan
  crystal = experiment.crystal

  # Get some stuff from the models
  A = matrix.sqr(crystal.get_A())
  s0 = matrix.col(beam.get_s0()).normalize()
  wavelength0 = beam.get_wavelength()
  m2 = matrix.col(goniometer.get_rotation_axis())
  fixed_rotation = matrix.sqr(goniometer.get_fixed_rotation())
  setting_rotation = goniometer.get_setting_rotation()
  dphi = scan.get_oscillation_range(deg=False)

  # Generate random wavelengths from a normal distribution
  if wavelength_spread > 0:
    wavelength = normal(wavelength0, wavelength_spread)
  else:
    wavelength = wavelength0

  # Create the incident beam vector
  s0 = s0 / wavelength

  # draw beam directions from a von mises fisher distriubtion
  if beam_divergence > 0:
    s0 = vonmises_fisher(s0, 1.0 / beam_divergence)

  # Get the miller index and entering flag
  h = reflection['miller_index']
  entering = reflection['entering']

  # Get the reciprocal lattice point
  r0 = fixed_rotation * (A * matrix.col(h))

  # Draw reciprocal lattice vectors from a multi variate normal distribution
  if rlp_covariance is not None:
    r0 = matrix.col(multivariate_normal(r0, matrix.sqr(rlp_covariance).as_list_of_lists()))

  # Draw reciprocal lattice vectors from a von mises distribution (spherical cap)
  if angular_spread > 0:
    r0 = vonmises_fisher(r0, 1.0 / angular_spread)

  # Create the ray predictor
  predictor = ScanStaticRayPredictor(
    s0,
    m2,
    fixed_rotation,
    setting_rotation,
    dphi)

  # Predict the rays for the reciprocal lattice point
  rays = predictor(r0)

  # Select the ray that matches the input reflection
  s1 = None
  phi = None
  for ray in rays:
    if ray.entering == entering:
      s1 = ray.s1
      phi = ray.angle

  assert s1 != None and phi != None

  # Get the ray intersection and image number
  x, y = detector[0].get_ray_intersection_px(s1)
  z = scan.get_array_index_from_angle(phi, deg=False)
  return x, y, z


def simulate(experiment,
             reflections,
             intensity=1000,
             wavelength_spread=0.0,
             rlp_covariance=None,
             angular_spread=0.0,
             beam_divergence=0.0):
  '''
  Simulate diffraction image from a number of reflections

  '''

  # Get the image size and number of images
  xsize, ysize = experiment.detector[0].get_image_size()
  zsize = experiment.scan.get_num_images()

  # Allocate the image
  image = flex.double(flex.grid(zsize, ysize, xsize))

  # Loop through the reflections
  for i in range(len(reflections)):

    print("Reflection %d / %d" % (i, len(reflections)))

    # Skip the 0,0,0 reflection
    h = reflections[i]['miller_index']
    if h == (0,0,0):
      continue

    # Simulate a number of counts from each reflection
    for j in range(intensity):

      # Try to simulate a count
      try:
        x, y, z = simulate_reflection(
          experiment,
          reflections[i],
          wavelength_spread = wavelength_spread,
          rlp_covariance    = rlp_covariance,
          angular_spread    = angular_spread,
          beam_divergence   = beam_divergence)
      except AssertionError:
        continue

      # Add count to the image
      if x >= 0 and y >= 0 and z >= 0 and x < xsize and y < ysize and z < zsize:
        image[int(z), int(y), int(x)] += 1

  # Dump the image
  import cPickle as pickle
  pickle.dump(image, open("image.pickle", "w"))


def display():
  import cPickle as pickle
  from dials.array_family import flex
  image = pickle.load(open("image.pickle"))

  from matplotlib import pylab, cm
  for k in range(image.all()[0]):
    data = image[k:k+1,:,:]
    print(flex.max(data))
    pylab.imshow(data.as_numpy_array()[0,:,:], cmap=cm.Greys_r)
    pylab.show()


phil_scope = parse('''

  wavelength_spread = 0.0
    .type = float(value_min=0.0)
    .help = "Wavelength sigma"

  beam_divergence = 0.0
    .type = float(value_min=0.0)
    .help = "Beam divergence"

  angular_spread = 0.0
    .type = float(value_min=0.0)
    .help = "Angular spread"

  rlp_covariance = None
    .type = floats(9)
    .help = "Reciprocal lattice point covariance matrix"

  intensity = 1000
    .type = int(value_min=0)
    .help = "Intensity of each spot"

''')


if __name__ == '__main__':

  from dials.array_family import flex
  from dials.util.options import OptionParser
  from dials.util.options import flatten_experiments, flatten_reflections
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
    read_experiments=True,
    read_reflections=True)

  # Parse the command line
  params, options = parser.parse_args(show_diff_phil=False)

  # Ensure we have a data block
  experiments = flatten_experiments(params.input.experiments)
  reflections = flatten_reflections(params.input.reflections)

  # Predict reflections if none are given
  if reflections is None or len(reflections) == 0:
    print("Predicting reflections")
    reflections = flex.reflection_table.from_predictions(experiments[0])
  else:
    reflections = reflections[0]

  # Siulate the reflections
  simulate(
    experiments[0],
    reflections,
    intensity         = params.intensity,
    wavelength_spread = params.wavelength_spread,
    rlp_covariance    = params.rlp_covariance,
    angular_spread    = params.angular_spread,
    beam_divergence   = params.beam_divergence)

  #display()
