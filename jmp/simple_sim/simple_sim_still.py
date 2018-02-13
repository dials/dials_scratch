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



def simulate_pixel(i, j,
                   experiment,
                   wavelength_spread = 0.0,
                   rlp_mosaicity    = [0,0,0,0,0,0],
                   angular_spread    = [0,0,0]):
  '''
  Simulate the counts in a single pixel

  '''
  from scitbx import matrix
  from dials.algorithms.spot_prediction import PixelToMillerIndex
  from numpy.random import normal, multivariate_normal
  from math import floor

  # Get some stuff from the experiment
  detector = experiment.detector
  panel = detector[0]
  beam = experiment.beam
  crystal = experiment.crystal
  A = matrix.sqr(crystal.get_A())
  U = matrix.sqr(crystal.get_U())
  s0 = matrix.col(beam.get_s0())
  
  
  # The coordinates in the pixel
  i0 = i
  i1 = i+1
  ic = i+0.5
  j0 = j
  j1 = j+1
  jc = j+0.5

  # Index the pixel
  pixel_to_miller_index = PixelToMillerIndex(beam, detector, crystal)
  h_frac = pixel_to_miller_index.h(0, ic, jc)
  h = matrix.col((
    int(floor(h_frac[0]+0.5)),
    int(floor(h_frac[1]+0.5)),
    int(floor(h_frac[2]+0.5))))

  # The reciprocal lattice vector
  rlp = A * h

  # Construct the matrix for the wavelength spread
  L = matrix.sqr((
    0, 0, 0,
    0, 0, 0,
    0, 0, wavelength_spread))

  # Construct the matrix for the angular spread
  W = matrix.sqr((
    angular_spread[0], 0, 0,
    angular_spread[1], angular_spread[2], 0,
    0, 0, 0))

  # Construct the convolution. Both components scale with the length of the
  # reciprocal lattice vector so multiply this here.
  #LLWW = (L*L.transpose() + W*W.transpose()) * rlp.length()
  LW = (L + W)

  # The mosaicity matrix is defined in the reciprocal lattice frame.
  # Put the mosaicity matrix in the lab frame
  M = matrix.sqr((
    rlp_mosaicity[0], 0, 0,
    rlp_mosaicity[1], rlp_mosaicity[2], 0,
    rlp_mosaicity[3], rlp_mosaicity[4], rlp_mosaicity[5]))
  Sigma_M = U * (M * M.transpose()) * U.transpose()

  # THe coordinate system
  s2 = s0 + rlp
  e1 = s2.cross(s0).normalize()
  e2 = -e1.cross(rlp).normalize()
  e3 = rlp.normalize()

  # Put the LLWW matrix in the lab frame
  E = matrix.sqr((e1[0], e2[0], e3[0],
                  e1[1], e2[1], e3[1],
                  e1[2], e2[2], e3[2]))
  Sigma_LW = E * (rlp.length()*LW*LW.transpose()) * E.transpose()
  #Sigma_LW = E.transpose() * LLWW * E

  # The full covariance matrix in the lab frame
  assert E.is_r3_rotation_matrix()
  Sigma = Sigma_M + Sigma_LW
  # print LW
  # print Sigma_LW
  # print Sigma_M
  # print Sigma
  # exit(0)
  # print rlp
  # print s0
  # print L
  # print E
  # print Sigma
  #Sigma = Sigma_LW + Sigma_M
  # print Sigma
  # print Sigma.inverse()
  
  # Get beam vectors at each corner
  s0_length = s0.length()
  s00 = matrix.col(panel.get_pixel_lab_coord((i,j))).normalize() * s0_length
  s01 = matrix.col(panel.get_pixel_lab_coord((i+1,j))).normalize() * s0_length
  s10 = matrix.col(panel.get_pixel_lab_coord((i,j+1))).normalize() * s0_length
  s11 = matrix.col(panel.get_pixel_lab_coord((i+1,j+1))).normalize() * s0_length

  # Compute the pixel area
  A = s00
  B = s01
  C = s11
  D = s10
  area = ((A-B).cross(A-D)).length() / 2.0 + ((C-B).cross(C-D)).length() / 2.0

  I = 0
  for jj in range(1):
    for ii in range(1):
      iii = i + (ii+0.5) #/ 10
      jjj = j + (jj+0.5) #/ 10
      ss = matrix.col(panel.get_pixel_lab_coord((iii,jjj))).normalize() * s0_length
      I += normal_pdf(ss, s0+rlp, Sigma)

  I = area * I #/ 100
  return I


def simulate(experiment,
             intensity=1000,
             wavelength_spread=0.0,
             rlp_mosaicity=[0,0,0,0,0,0],
             angular_spread=[0,0,0]):
  '''
  Simulate diffraction image from a number of reflections

  '''
  from dials_scratch.jmp.stills import StillsSimulator
  from numpy.random import poisson

  # Get the image size and number of images
  xsize, ysize = experiment.detector[0].get_image_size()

  # Allocate the image
  image = flex.double(flex.grid(ysize, xsize))
  model = flex.double(flex.grid(ysize, xsize))

  for j in range(ysize//3):
    print j
    for i in range(xsize//3):
      try:
        model[j,i] = simulate_pixel(
          i, j,
          experiment,
          wavelength_spread = wavelength_spread,
          rlp_mosaicity    = rlp_mosaicity,
          angular_spread    = angular_spread)
        image[j,i] = poisson(model[j,i]*intensity)
      except Exception:
        pass

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
  pylab.imshow(data.as_numpy_array(), cmap=cm.Greys_r)
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

  intensity = 1000
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
