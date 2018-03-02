from __future__ import division
from libtbx.phil import parse
from math import sqrt, exp, pi

def normal_pdf(x, mu, sigma):
  '''
  Multi variate normal

  '''
  N = len(mu)
  sigma_inv = sigma.inverse()
  A = (1.0 / (sqrt((2*pi)**N *sigma.determinant())))
  B = ((x-mu).transpose()*sigma_inv*(x-mu))[0]
  return A * exp(-0.5 * B)



def simulate_pixel(i, j, pixel_to_miller_index, simulator):
  '''
  Simulate the counts in a single pixel

  '''
  from scitbx import matrix
  from numpy.random import normal, multivariate_normal
  from math import floor

  # The coordinates in the pixel
  ic = i+0.5
  jc = j+0.5

  # Index the pixel
  h_frac = pixel_to_miller_index.h(0, ic, jc)
  h = (
    int(floor(h_frac[0]+0.5)),
    int(floor(h_frac[1]+0.5)),
    int(floor(h_frac[2]+0.5)))

  return simulator.integrate_pixel(h, 0, i, j)


def simulate(experiment,
             intensity=1000,
             wavelength_spread=0.0,
             rlp_mosaicity=[0,0,0,0,0,0],
             angular_spread=[0,0,0]):
  '''
  Simulate diffraction image from a number of reflections

  '''
  from dials.algorithms.spot_prediction import PixelToMillerIndex
  from dials_scratch.jmp.stills import StillsSimulator
  from numpy.random import poisson
  from scitbx import matrix

  # Get the image size and number of images
  xsize, ysize = experiment.detector[0].get_image_size()

  # Allocate the image
  image = flex.double(flex.grid(ysize, xsize))
  model = flex.double(flex.grid(ysize, xsize))

  L = matrix.sqr((
    0, 0, 0,
    0, 0, 0,
    0, 0, wavelength_spread))
  Sigma_L = L * L.transpose()

  # Construct the matrix for the angular spread
  W = matrix.sqr((
    angular_spread[0], 0, 0,
    angular_spread[1], angular_spread[2], 0,
    0, 0, 0))
  Sigma_W = W * W.transpose()

  M = matrix.sqr((
    rlp_mosaicity[0], 0, 0,
    rlp_mosaicity[1], rlp_mosaicity[2], 0,
    rlp_mosaicity[3], rlp_mosaicity[4], rlp_mosaicity[5]))
  Sigma_M = M * M.transpose()

  simulator = StillsSimulator(
    experiment.beam,
    experiment.detector,
    experiment.crystal,
    Sigma_L,
    Sigma_W,
    Sigma_M,
    10)

  pixel_to_miller_index = PixelToMillerIndex(
    experiment.beam,
    experiment.detector,
    experiment.crystal)

  for j in range(ysize):
    print j
    for i in range(xsize):
      model[j,i] = simulate_pixel(i, j, pixel_to_miller_index, simulator)
      image[j,i] = poisson(model[j,i]*intensity)

  # Dump the image
  import cPickle as pickle
  pickle.dump(model, open("model.pickle", "w"))
  pickle.dump(image, open("image.pickle", "w"))


def display():
  import cPickle as pickle
  from dials.array_family import flex
  image = pickle.load(open("model.pickle"))

  from matplotlib import pylab, cm
  data = image
  print flex.max(data)
  pylab.imshow(data.as_numpy_array(), cmap=cm.Greys_r, vmax=1)
  pylab.show()

  # image = pickle.load(open("image.pickle"))

  # data = image
  # print flex.max(data)
  # pylab.imshow(data.as_numpy_array(), cmap=cm.Greys_r)
  # pylab.show()


phil_scope = parse('''

  wavelength_spread = 0.0
    .type = float(value_min=0.0)
    .help = "Wavelength sigma"

  angular_spread = 0,0,0
    .type = floats(size=3)
    .help = "Angular spread lower triangular matrix"

  rlp_mosaicity = 0,0,0,0,0,0
    .type = floats(6)
    .help = "Reciprocal lattice point covariance lower triangular matrix"

  intensity = 10
    .type = int(value_min=0)
    .help = "Intensity of each spot"

''')


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

  # Siulate the reflections
  simulate(
    experiments[0],
    intensity         = params.intensity,
    wavelength_spread = params.wavelength_spread,
    rlp_mosaicity    = params.rlp_mosaicity,
    angular_spread    = params.angular_spread)

  display()
