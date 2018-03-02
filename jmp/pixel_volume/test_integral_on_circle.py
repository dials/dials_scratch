
from __future__ import division


def normal(x, mu, sigma):
  '''
  Multi variate normal

  '''
  from math import sqrt, exp, pi
  N = len(mu)
  sigma_inv = sigma.inverse()
  A = (1.0 / (sqrt((2*pi)**N *sigma.determinant())))
  B = ((x-mu).transpose()*sigma_inv*(x-mu))[0]
  return A * exp(-0.5 * B)


def f(r, theta, mu, sigma):
  '''
  Multi variate normal at polar coords

  '''
  from math import cos, sin
  from scitbx import matrix
  x = r*cos(theta)
  y = r*sin(theta)
  return normal(matrix.col((x, y)), mu, sigma)

def g(r, theta, mu, sigma):
  '''
  The integral we want to evalute (the mean theta)

  '''
  from cmath import exp
  return r * exp(complex(0,theta)) * f(r, theta, mu, sigma)


def integrate(angle=0):
  from scitbx import matrix
  from math import pi, cos, sin

  # Spherical RLP at R = 10
  mu = matrix.col((-1, 0))
  sigma = matrix.sqr((0.25*0.25, 0, 0, 0.1*0.1))

  R = matrix.sqr((cos(angle), -sin(angle), sin(angle), cos(angle)))

  sigma = R * sigma * R.transpose()

  # Spherical RLP at R = 10
  # r = 10
  # theta = 10 * pi / 180.0
  # mu = matrix.col((r*cos(theta), r*sin(theta)))
  # sigma = matrix.sqr((0.5, 0, 0, 1.0))
  # print mu.angle(matrix.col((1, 0)), deg=True)
  N = 1000

  a = -pi
  b = pi

  r = 1

  I = g(r, a, mu, sigma) / 2.0 + g(r, b, mu, sigma) / 2.0
  for i in range(1, N):
    I += g(r, a + i*(b-a)/N, mu, sigma)
  I = I * (b-a) / N

  I2 = r*f(r, a, mu, sigma) / 2.0 + r*f(r, b, mu, sigma) / 2.0
  for i in range(1, N):
    I2 += r*f(r, a + i*(b-a)/N, mu, sigma)
  I2 = I2 * (b-a) / N

  return  I / I2



from math import pi, atan2
y = integrate()
print atan2(y.imag, y.real) * 180.0 / pi

exit(0)

#print (pi*10*10) * integrate() * 180 / pi
X = []
Y = []
for angle in range(0, 180):
  x = angle
  y = integrate(angle * pi / 180.0) * 180.0 / pi
  print x, y
  X.append(x)
  Y.append(y)

from matplotlib import pylab
pylab.scatter(X, Y)
pylab.show()
# print (pi*10*10)/180.0 * integrate() * 180.0 / pi


exit(0)


from matplotlib import pylab
from numpy.random import multivariate_normal
from math import sin, cos, pi
from scitbx import matrix

X = []
Y = []
X1 = []
Y1 = []
for i in range(1000):
  r = 10
  theta = i*2*pi / 1000.0
  X.append(r*cos(theta))
  Y.append(r*sin(theta))

  r = 10
  theta = -10 * pi / 180.0
  mu = matrix.col((r*cos(theta), r*sin(theta)))
  sigma = matrix.sqr((1, 0, 0, 1))
  p = multivariate_normal(mu, sigma.as_list_of_lists())

  X1.append(p[0])
  Y1.append(p[1])

pylab.plot(X, Y, color='red')
pylab.scatter(X1, Y1, color='blue')
pylab.show()
