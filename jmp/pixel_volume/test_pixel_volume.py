
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

  N = 100000

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
      s = x.rotate_around_origin(axis=axis, angle=angle, deg=False) * scc.length()

    px, py = panel.get_ray_intersection_px(s)

    px, py = int(floor(px)), int(floor(py))

    if px == i and py == j:
      count += 1


  box_area = 2*pi*s0.length()**2 * (1-cos(phi_max))

  V = box_area * count / float(N)

  return V

def monte_carlo_pixel_volume2(experiment, i, j):

  from dials.array_family import flex
  from scitbx import matrix
  from random import uniform
  from math import acos, pi, cos, sin, floor, sqrt, atan2
  from dials.algorithms.spot_prediction import RotationAngles

  beam = experiment.beam
  goniometer = experiment.goniometer
  panel = experiment.detector[0]
  scan = experiment.scan

  s0_length = 1.0 / beam.get_wavelength()

  phi0, phi1 = scan.get_oscillation(deg=False)

  DPHI = phi1

  m2 = matrix.col(goniometer.get_rotation_axis_datum())
  s0 = matrix.col(beam.get_s0())

  N = 100000

  s00 = matrix.col(panel.get_pixel_lab_coord((i,j))).normalize() * s0_length
  s01 = matrix.col(panel.get_pixel_lab_coord((i+1,j))).normalize() * s0_length
  s10 = matrix.col(panel.get_pixel_lab_coord((i,j+1))).normalize() * s0_length
  s11 = matrix.col(panel.get_pixel_lab_coord((i+1,j+1))).normalize() * s0_length
  scc = matrix.col(panel.get_pixel_lab_coord((i+0.5,j+0.5))).normalize() * s0_length

  p0_00 = s00 - s0
  p0_01 = s01 - s0
  p0_10 = s10 - s0
  p0_11 = s11 - s0

  p1_00 = p0_00.rotate_around_origin(axis=m2, angle=-DPHI, deg=False)
  p1_01 = p0_01.rotate_around_origin(axis=m2, angle=-DPHI, deg=False)
  p1_10 = p0_10.rotate_around_origin(axis=m2, angle=-DPHI, deg=False)
  p1_11 = p0_11.rotate_around_origin(axis=m2, angle=-DPHI, deg=False)

  X0, Y0, Z0 = zip(*[p0_00, p0_01, p0_10, p0_11])
  X1, Y1, Z1 = zip(*[p1_00, p1_01, p1_10, p1_11])

  XX, YY, ZZ = X0 + X1, Y0 + Y1, Z0 + Z1
  Xmin, Xmax = min(XX), max(XX)
  Ymin, Ymax = min(YY), max(YY)
  Zmin, Zmax = min(ZZ), max(ZZ)

  Xrange = (Xmax - Xmin) * 1.1
  Yrange = (Ymax - Ymin) * 1.1
  Zrange = (Zmax - Zmin) * 1.1
  Xmid = (Xmax + Xmin) / 2.0
  Ymid = (Ymax + Ymin) / 2.0
  Zmid = (Zmax + Zmin) / 2.0

  Xmin = Xmid - Xrange / 2.0
  Ymin = Ymid - Yrange / 2.0
  Zmin = Zmid - Zrange / 2.0
  Xmax = Xmid + Xrange / 2.0
  Ymax = Ymid + Yrange / 2.0
  Zmax = Zmid + Zrange / 2.0

  count = 0
  min_angle = 10
  max_angle = -10
  X = []
  Y = []
  Z = []
  for n in range(N):
    x = uniform(Xmin, Xmax)
    y = uniform(Ymin, Ymax)
    z = uniform(Zmin, Zmax)

    r0 = matrix.col((x, y, z))

    try:
      angles = RotationAngles(s0, m2)(r0)
    except Exception:
      continue

    for angle in angles:
      #print angle, DPHI
      if angle >= 0 and angle <= DPHI:
        if angle > max_angle:
          max_angle = angle
        if angle < min_angle:
          min_angle = angle
        p = r0.rotate_around_origin(axis=m2, angle=angle, deg=False)
        s = s0 + p
        px, py = panel.get_ray_intersection_px(s)

        px, py = int(floor(px)), int(floor(py))
        #print px, py, i, j
        if px == i and py == j:
          count += 1
          X.append(x)
          Y.append(y)
          Z.append(z)

  box_volume = (Xmax - Xmin)*(Ymax-Ymin)*(Zmax-Zmin)
  V = box_volume * count / float(N)
  # print box_volume, V, count, min_angle * 180.0 / pi, max_angle*180.0 / pi, DPHI*180.0 /pi

  # from mpl_toolkits.mplot3d import Axes3D
  # from matplotlib import pylab
  # fig = pylab.figure()
  # ax = fig.add_subplot(111, projection='3d')
  # ax.scatter(X, Y, Z, color="blue")

  # # Draw box around volume
  # ax.plot3D([Xmin, Xmax], [Ymin, Ymin], [Zmin, Zmin], color="red")
  # ax.plot3D([Xmin, Xmax], [Ymax, Ymax], [Zmin, Zmin], color="red")
  # ax.plot3D([Xmin, Xmax], [Ymin, Ymin], [Zmax, Zmax], color="red")
  # ax.plot3D([Xmin, Xmax], [Ymax, Ymax], [Zmax, Zmax], color="red")
  # ax.plot3D([Xmin, Xmin], [Ymin, Ymax], [Zmin, Zmin], color="red")
  # ax.plot3D([Xmax, Xmax], [Ymin, Ymax], [Zmin, Zmin], color="red")
  # ax.plot3D([Xmin, Xmin], [Ymin, Ymax], [Zmax, Zmax], color="red")
  # ax.plot3D([Xmax, Xmax], [Ymin, Ymax], [Zmax, Zmax], color="red")
  # ax.plot3D([Xmin, Xmin], [Ymin, Ymin], [Zmin, Zmax], color="red")
  # ax.plot3D([Xmax, Xmax], [Ymin, Ymin], [Zmin, Zmax], color="red")
  # ax.plot3D([Xmin, Xmin], [Ymax, Ymax], [Zmin, Zmax], color="red")
  # ax.plot3D([Xmax, Xmax], [Ymax, Ymax], [Zmin, Zmax], color="red")

  # # Draw R0 and R0 rotated


  # pylab.show()

  # x00, y00, z00 = p0_00
  # x01, y01, z01 = p0_01
  # x10, y10, z10 = p0_10
  # x11, y11, z11 = p0_11

  # x00, y00, z00 = p0_00
  # x01, y01, z01 = p0_01
  # x10, y10, z10 = p0_10
  # x11, y11, z11 = p0_11


#
  #R00 = sqrt(z00**2 + y00**2)
  #R01 = sqrt(z01**2 + y01**2)
  #R10 = sqrt(z10**2 + y10**2)
  #R11 = sqrt(z11**2 + y11**2)

  #P00 = atan2(z00, y00)
  #P01 = atan2(z01, y01)
  #P10 = atan2(z10, y10)
  #P11 = atan2(z11, y11)

  #xmin, xmax = min(x00, x01, x10, x11), max(x00, x01, x10, x11)
  #Rmin, Rmax = min(R00, R01, R10, R11), max(R00, R01, R10, R11)
  #Pmin, Pmax = min(P00, P01, P10, P11), max(P00, P01, P10, P11)

  ## print xmin, xmax, Rmin, Rmax, Pmin, Pmax, DPHI

  #Pmax += DPHI
  #count = 0
  #min_angle = 10
  #max_angle = -10
  #for n in range(N):

  #  x = uniform(xmin, xmax)
  #  R = uniform(Rmin, Rmax)
  #  P = uniform(Pmin, Pmax)

  #  r0 = matrix.col((x, R*cos(P), R*sin(P)))
  #  # print r0.angle(scc - s0, deg=True)
  #  angles = RotationAngles(s0, m2)(r0)
  #  print P, angles
  #  for angle in angles:
  #    #print angle, DPHI
  #    if angle >= 0 and angle <= DPHI:
  #      if angle > max_angle:
  #        max_angle = angle
  #      if angle < min_angle:
  #        min_angle = angle
  #      p = r0.rotate_around_origin(axis=m2, angle=angle, deg=False)
  #      s = s0 + p
  #      px, py = panel.get_ray_intersection_px(s)

  #      px, py = int(floor(px)), int(floor(py))
  #      #print px, py, i, j
  #      if px == i and py == j:
  #        count += 1


  ##box_area = 2*pi*s0.length()**2 * (1-cos(phi_max))
  #box_volume = (Rmax**2 - Rmin**2) * (Pmax - Pmin)*(xmax-xmin)/2.0
  ##print box_volume, count, N
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

  phi0, phi1 = scan.get_oscillation(deg=False)

  DPHI = phi1

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

  Linv = m2.dot(s1.cross(s0)) / s0_length**2

  V = abs(A * DPHI * Linv * s0_length)

  return V

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

  phi0, phi1 = scan.get_oscillation(deg=False)

  DPHI = phi1

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
  Linv = m2.dot(s1.cross(s0)) / s0_length**2

  V = abs(A * DPHI * Linv * s0_length)

  return V



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
    V = monte_carlo_pixel_volume2(experiment, i, j)
    V1 = pixel_polygon_area(experiment, i, j)
    V2 = pixel_volume(experiment, i, j)
    print i, j, V, V1, V2, time() - st
    # exit(0)
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
