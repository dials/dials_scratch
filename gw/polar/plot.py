from __future__ import print_function
import matplotlib
matplotlib.use('Agg')

import numpy
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as axes3d

STEPS=100

def sph_harm_surf(l, m):
  from scitbx import math
  from scitbx.array_family import flex
  import math as pymath
  lfg = math.log_factorial_generator(2 * l + 1)
  nsssphe = math.nss_spherical_harmonics(l, 50000, lfg)

  theta = numpy.linspace(0, 2 * numpy.pi, 2*STEPS)
  phi = numpy.linspace(0, numpy.pi, STEPS)

  THETA, PHI = numpy.meshgrid(theta, phi)
  R = numpy.cos(PHI**2)
  C = numpy.empty(THETA.shape, dtype=str)

  sqrt2 = pymath.sqrt(2)

  for it, t in enumerate(theta):
    for ip, p in enumerate(phi):
      Ylm = nsssphe.spherical_harmonic(l, abs(m), p, t)
      if m < 0:
        r = sqrt2 * ((-1) ** m) * Ylm.imag
      elif m == 0:
        assert(Ylm.imag == 0.0)
        r = Ylm.real
      else:
        r = sqrt2 * ((-1) ** m) * Ylm.real
      R[ip, it] = pymath.fabs(r)
      if r < 0:
        C[ip, it] = 'y'
      else:
        C[ip, it] = 'b'

  X = R * numpy.sin(PHI) * numpy.cos(THETA)
  Y = R * numpy.sin(PHI) * numpy.sin(THETA)
  Z = R * numpy.cos(PHI)

  fig = plt.figure()
  ax = fig.add_subplot(1,1,1, projection='3d')
  plot = ax.plot_surface(
      X, Y, Z, rstride=1, cstride=1, facecolors=C,
      linewidth=0, antialiased=True, alpha=0.5)

  print('Saving %s...' % ('ylm%d%d.png' % (l, m)))
  plt.savefig('ylm%d%d.png' % (l, m))

import sys

lmax = int(sys.argv[1])

for l in range(1, lmax+1):
  for m in range(-l, l+1):
    sph_harm_surf(l, m)
