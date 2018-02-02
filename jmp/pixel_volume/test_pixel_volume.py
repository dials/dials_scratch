
from __future__ import division

def monte_carlo_pixel_volume(experiment, i, j):

  from dials.array_family import flex
  from scitbx import matrix
  from random import uniform
  from math import acos, pi, cos, sin, floor

  beam = experiment.beam
  goniometer = experiment.goniometer
  panel = experiment.detector[0]
  scan = experiment.scan

  s0_length = 1.0 / beam.get_wavelength()

  phi0, phi1 = scan.get_oscillation()

  DPHI = phi1 - phi0

  m2 = matrix.col(goniometer.get_rotation_axis_datum())
  s0 = matrix.col(beam.get_s0())

  N = 1000000

  s00 = matrix.col(panel.get_pixel_lab_coord((i,j))).normalize() * s0_length
  s01 = matrix.col(panel.get_pixel_lab_coord((i+1,j))).normalize() * s0_length
  s10 = matrix.col(panel.get_pixel_lab_coord((i,j+1))).normalize() * s0_length
  s11 = matrix.col(panel.get_pixel_lab_coord((i+1,j+1))).normalize() * s0_length
  scc = matrix.col(panel.get_pixel_lab_coord((i+0.5,j+0.5))).normalize() * s0_length

  phi0 = scc.angle(s00)
  phi1 = scc.angle(s01)
  phi2 = scc.angle(s11)
  phi3 = scc.angle(s10)

  phi_max = max(phi0, phi1, phi2, phi3)

  phi_max = phi_max * 1.1

  vmin = (cos(phi_max) + 1)/2.0

  count = 0
  for n in range(N):

    u = uniform(0,1)
    v = uniform(vmin,1)

    theta = 2*pi*u
    phi = acos(2*v - 1)

    a = matrix.col((0, 0, 1))
    b = scc.normalize()
    r = s0.length()

    x = matrix.col((
      r*sin(phi)*cos(theta),
      r*sin(phi)*sin(theta),
      r*cos(phi)))

    axis = a.cross(b)
    angle = acos(a.dot(b))

    if angle < 1e-7:
      s = x
    elif abs(angle - pi) < 1e-7:
      s = -x
    else:
      s = x.rotate_around_origin(axis=axis, angle=angle) * scc.length()

    px, py = panel.get_ray_intersection_px(s)

    px, py = int(floor(px)), int(floor(py))

    if px == i and py == j:
      count += 1


  box_area = 2*pi*s0.length()**2 * (1-cos(phi_max))

  V = box_area * count / float(N)

  return V


def pixel_polygon_area(experiment, i, j):

  from dials.array_family import flex
  from scitbx import matrix
  from random import uniform
  from math import acos, pi

  beam = experiment.beam
  goniometer = experiment.goniometer
  panel = experiment.detector[0]
  scan = experiment.scan

  s0_length = 1.0 / beam.get_wavelength()

  phi0, phi1 = scan.get_oscillation()

  DPHI = phi1 - phi0

  m2 = matrix.col(goniometer.get_rotation_axis_datum())
  s0 = matrix.col(beam.get_s0())

  s00 = matrix.col(panel.get_pixel_lab_coord((i,j))).normalize() * s0_length
  s01 = matrix.col(panel.get_pixel_lab_coord((i+1,j))).normalize() * s0_length
  s10 = matrix.col(panel.get_pixel_lab_coord((i,j+1))).normalize() * s0_length
  s11 = matrix.col(panel.get_pixel_lab_coord((i+1,j+1))).normalize() * s0_length
  scc = matrix.col(panel.get_pixel_lab_coord((i+0.5,j+0.5))).normalize() * s0_length

  s1 = scc

  A = s00
  B = s01
  C = s11
  D = s10

  A = ((A-B).cross(A-D)).length() / 2.0 + ((C-B).cross(C-D)).length() / 2.0

  return A

def pixel_volume(experiment, i, j):

  from dials.array_family import flex
  from scitbx import matrix
  from random import uniform
  from math import acos, pi, cos, sin

  beam = experiment.beam
  goniometer = experiment.goniometer
  panel = experiment.detector[0]
  scan = experiment.scan

  s0_length = 1.0 / beam.get_wavelength()

  phi0, phi1 = scan.get_oscillation()

  DPHI = phi1 - phi0

  m2 = matrix.col(goniometer.get_rotation_axis_datum())
  s0 = matrix.col(beam.get_s0())

  s00 = matrix.col(panel.get_pixel_lab_coord((i,j))).normalize() * s0_length
  s01 = matrix.col(panel.get_pixel_lab_coord((i+1,j))).normalize() * s0_length
  s10 = matrix.col(panel.get_pixel_lab_coord((i,j+1))).normalize() * s0_length
  s11 = matrix.col(panel.get_pixel_lab_coord((i+1,j+1))).normalize() * s0_length
  scc = matrix.col(panel.get_pixel_lab_coord((i+0.5,j+0.5))).normalize() * s0_length

  s1 = scc

  A = s00.normalize()
  B = s01.normalize()
  C = s11.normalize()
  D = s10.normalize()

  cosa = B.dot(C)
  cosb = A.dot(C)
  cosc = A.dot(B)

  sina = sin(acos(cosa))
  sinb = sin(acos(cosb))
  sinc = sin(acos(cosc))

  cosA = (cosa - cosb*cosc) / (sinb*sinc)
  cosB = (cosb - cosc*cosa) / (sinc*sina)
  cosC = (cosc - cosa*cosb) / (sina*sinb)

  Aa = acos(cosA)
  Ba = acos(cosB)
  Ca = acos(cosC)

  A1 = (Aa + Ba + Ca - pi)*s0_length**2

  cosa = D.dot(C)
  cosb = A.dot(C)
  cosc = A.dot(D)

  sina = sin(acos(cosa))
  sinb = sin(acos(cosb))
  sinc = sin(acos(cosc))

  cosA = (cosa - cosb*cosc) / (sinb*sinc)
  cosB = (cosb - cosc*cosa) / (sinc*sina)
  cosC = (cosc - cosa*cosb) / (sina*sinb)

  Aa = acos(cosA)
  Ba = acos(cosB)
  Ca = acos(cosC)

  A2 = (Aa + Ba + Ca - pi)*s0_length**2

  A = A1+A2
  #Linv = m2.dot(s1.cross(s0)) / s0_length**2

  #V = abs(A * DPHI * Linv * s0_length)

  return A



def monte_carlo_pixel_volume_all_pixels(experiment):

  from dials.array_family import flex
  from time import time
  from random import randint
  panel = experiment.detector[0]
  xsize, ysize = panel.get_image_size()
  volume = flex.double(flex.grid(ysize, xsize))
  for n in range(1000):
    i = randint(0, xsize-1)
    j = randint(0, ysize-1)
    st = time()
    V = monte_carlo_pixel_volume(experiment, i, j)
    V1 = pixel_polygon_area(experiment, i, j)
    V2 = pixel_volume(experiment, i, j)
    print i, j, V, V1, V2, time() - st
    volume[j,i] = V

  return volume



def pixel_volume_all_pixels(experiment):

  from dials.array_family import flex
  xsize, ysize = panel.get_image_size()
  volume = flex.double(flex.grid(ysize, xsize))

  for j in range(ysize):
    print j
    for i in range(xsize):

      volume[j,i] = pixel_volume(experiment, i, j)

  return volume


def display(data):

  class Formatter(object):
    def __init__(self, im):
        self.im = im
    def __call__(self, x, y):
        z = self.im.get_array()[int(y), int(x)]
        print z
        return 'x={:.01f}, y={:.01f}, z={:f}'.format(x, y, z)

  from matplotlib import pylab
  fig, ax = pylab.subplots()
  im = ax.imshow(data.as_numpy_array(), interpolation='none')
  ax.format_coord = Formatter(im)
  pylab.show()

if __name__ == '__main__':

  import sys
  from dxtbx.model.experiment_list import ExperimentListFactory

  experiments = ExperimentListFactory.from_json_file(sys.argv[1])

  #volume1 = pixel_volume_all_pixels(experiments[0])
  volume2 = monte_carlo_pixel_volume_all_pixels(experiments[0])

  import cPickle as pickle

  #pickle.dump(volume1, open("volume.pickle", "w"))
  # pickle.dump(volume2, open("monte_carlo_volume.pickle", "w"))

  # volume1 = pickle.load(open("volume.pickle"))

  # display(volume1)
