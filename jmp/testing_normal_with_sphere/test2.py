

from math import sqrt, sin, cos, pi, exp

def normal_2d(x, mu, sigma):
  A = 1.0 / sqrt((2*pi)**2 * sigma.determinant())
  B = ((x - mu).transpose() * sigma.inverse() * (x - mu))[0]
  return A * exp(-0.5 * B)

def find_maximum(mu1, mu2, s1, s2):

  from scitbx import matrix
  mu = matrix.col((mu1, mu2))
  sigma = matrix.sqr((s1**2, 0,0,s2**2))

  t = []
  f = []

  a = -pi
  b = pi

  peak_theta = 0
  peak_value = 0

  for i in range(10000):
    theta = a + i*(b-a) / 10000
    x = matrix.col((cos(theta),sin(theta)))
    ft = normal_2d(x, mu, sigma)

    if ft > peak_value:
      peak_value = ft
      peak_theta = theta

    t.append(theta*180/pi)
    f.append(ft)

  x = cos(peak_theta)
  y = sin(peak_theta)

  print x, y


mu1 = 1.1*1 / sqrt(3)
mu2 = 1.1*1 / sqrt(3)
mu3 = 1.1*1 / sqrt(3)

s1 = 0.01
s2 = 0.02
s3 = 0.03

def func(x, mu1, mu2, mu3, s1, s2, s3):
  A = mu1**2 / (1 - x*s1**2)**2
  B = mu2**2 / (1 - x*s2**2)**2
  C = mu3**2 / (1 - x*s3**2)**2
  return A + B + C - 1


def df(x, mu1, mu2, mu3, s1, s2, s3):
  A = 2*s1**2 * mu1**2 / (1 - x*s1**2)**3
  B = 2*s2**2 * mu2**2 / (1 - x*s2**2)**3
  C = 2*s3**2 * mu3**2 / (1 - x*s3**2)**3
  return A + B + C

def d2f(x, mu1, mu2, mu3, s1, s2, s3):
  A = 6*s1**4 * mu1**2 / (1 - x*s1**2)**4
  B = 6*s2**4 * mu2**2 / (1 - x*s2**2)**4
  C = 6*s3**4 * mu3**2 / (1 - x*s3**2)**4
  return A + B + C


def compute_l0(x, mu, s):
  return (1-mu/x)/s**2

from scitbx import matrix
x = matrix.col((mu1, mu2, mu3)).normalize()

print "Max: ", 1.0 / s1**2
print "Max: ", 1.0 / s2**2
print "Max: ", 1.0 / s3**2

print compute_l0(x[0], mu1, s1)
print compute_l0(x[1], mu2, s2)
print compute_l0(x[2], mu3, s3)

l0 = compute_l0(x[0], mu1, s1)

for t in range(5):

  f = func(l0, mu1, mu2, mu3, s1, s2, s3)
  fp = df(l0, mu1, mu2, mu3, s1, s2, s3)
  fp2 = d2f(l0, mu1, mu2, mu3, s1, s2, s3)

  l0 = l0 - 2*f*fp / (2*fp**2 - f*fp2)

  x1 = mu1 / (1 - l0*s1**2)
  x2 = mu2 / (1 - l0*s2**2)
  x3 = mu3 / (1 - l0*s3**2)
  
  print l0, x1, x2, x3, x1**2 + x2**2 + x3**2

print tuple(x)

#find_maximum(mu1, mu2, s1, s2)

# X = []
# Y = []
# import numpy as np
# for x in np.arange(-100, 100, 0.01):

#   X.append(x)
#   Y.append(func(x, mu1, mu2, s1, s2))

# from matplotlib import pylab
# pylab.plot(X, Y)
# pylab.show()

