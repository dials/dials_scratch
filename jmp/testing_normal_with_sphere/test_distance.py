from __future__ import print_function

from scitbx import matrix
from math import exp, sqrt, pi

def normal_1d(x, mu, sigma_sq):
  A = (1.0/sqrt(2*pi*sigma_sq))
  B = (x - mu)**2 / sigma_sq
  return A * exp(-0.5 * B)

sigma = matrix.sqr((0.01, 0, 0, 0.01))
mu = matrix.col((0.9, 0))

mux, muy = mu
s11, s12, s21, s22 = sigma.inverse()

a = s11
b = s22
c = mux*(mux*s11 + muy*s12) + muy*(mux*s21 + muy*s22)

D = 0.5 * (a + b + 2*c)

print(D, normal_1d(D, 0, 1))

