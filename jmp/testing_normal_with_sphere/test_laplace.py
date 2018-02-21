from scitbx.linalg import eigensystem
from scitbx import matrix
from dials.array_family import flex
from math import pi, exp, sqrt, sin, atan2, cos, acos, erf
from random import uniform

def normal_1d(x, mu, sigma_sq):
  A = 1.0 / sqrt(2*pi * sigma_sq)
  B = ((x - mu)**2 / (sigma_sq))
  return A * exp(-0.5 * B)

def normal_2d(x, mu, sigma):
  A = 1.0 / sqrt((2*pi)**2 * sigma.determinant())
  B = ((x - mu).transpose() * sigma.inverse() * (x - mu))[0]
  return A * exp(-0.5 * B)

def f(x, mu1, mu2, s1, s2):
  A = mu1**2 / (1 - x*s1)**2
  B = mu2**2 / (1 - x*s2)**2
  return A + B - 1

def df(x, mu1, mu2, s1, s2):
  A = 2*s1 * mu1**2 / (1 - x*s1)**3
  B = 2*s2 * mu2**2 / (1 - x*s2)**3
  return A + B

def d2f(x, mu1, mu2, s1, s2):
  A = 6*s1 * mu1**2 / (1 - x*s1)**4
  B = 6*s2 * mu2**2 / (1 - x*s2)**4
  return A + B

def compute_l0(x, mu, s):
  return (1-mu[0]/x[0])/s[0]**2

def compute_mode(mu, sigma):

  # Compute the eigen decomposition of the covariance matrix
  eigen_decomposition = eigensystem.real_symmetric(sigma.as_flex_double_matrix())
  Q = matrix.sqr(eigen_decomposition.vectors())
  L = matrix.diag(eigen_decomposition.values())

  # To simplify calculation of the mode, we change basis
  mu_prime = Q.transpose() * mu

  # Compute the initial value of the lagrange multiplier
  l0 = compute_l0(mu_prime.normalize(), mu_prime, L)

  # Compute the value via Halley's method
  while True:
    mu1, mu2 = mu_prime
    s1, s2 = L[0], L[3]
    U = f(l0, mu1, mu2, s1, s2)
    Up = df(l0, mu1, mu2, s1, s2)
    Up2 = d2f(l0, mu1, mu2, s1, s2)
    l = l0 - 2*U*Up / (2*Up**2 - U*Up2)
    if abs(l - l0) < 1e-7:
      break
    l0 = l

  # Compute the mode
  x = matrix.col((
    mu_prime[0] / (1 - l * L[0]),
    mu_prime[1] / (1 - l * L[3])))

  # Now transform back into the original basis
  return Q * x

def compute_mean(mu, sigma):
  from matplotlib import pylab

  X = []
  Y = []
  XX = []
  YY = []

  peak_xx_value = 0
  peak_xx_coord = 0
  
  peak_yy_value = 0
  peak_yy_coord = 0

  I1 = matrix.col((0,0))
  I2 = 0
  for i in range(10000):
    theta = 0 + i*(pi/2)/10000
    x = cos(theta)
    y = sin(theta)
    v = matrix.col((x, y))
    g = normal_2d(v, mu, sigma)
    I1 += v * g
    I2 += g
    if v[0] * g > peak_xx_value:
      peak_xx_value = v[0] * g
      peak_xx_coord = x
    if v[1] * g > peak_yy_value:
      peak_yy_value = v[1] * g
      peak_yy_coord = x

    X.append(x)
    Y.append(y)
    XX.append(v[0] * g)
    YY.append(v[1] * g)

  print peak_xx_coord, peak_yy_coord

  peak = (0.7222307526374394, 0.6916521813344897)
  x0 = peak[0]
  A = sqrt(2*pi/abs(1/sigma[0]))*x0*exp(-0.5*(x0-mu[0])**2/sigma[0]**2)
  B = sqrt(2*pi/abs(1/sigma[0]))*exp(-0.5*(x0-mu[0])**2/sigma[0]**2)
  XX = A / B
  y0 = peak[1]
  A = sqrt(2*pi/abs(1/sigma[3]**2))*y0*exp(-0.5*(y0-mu[1])**2/sigma[3]**2)
  B = sqrt(2*pi/abs(1/sigma[3]**2))*exp(-0.5*(y0-mu[1])**2/sigma[3]**2)
  YY =  A / B

  peak = matrix.col((XX, YY)).normalize()
  print tuple(peak)

  peak = matrix.col((peak_xx_coord, peak_yy_coord)).normalize()
  print tuple(peak)

  pylab.plot(X, XX)
  pylab.plot(Y, YY)
  pylab.show()

  return (I1 / I2).normalize()

if __name__ == '__main__':

  mu = matrix.col((1, 1)).normalize() * 1.05
  sig = matrix.sqr((
    0.1, 0,
    0, 0.14))
  sigma = sig * sig.transpose()

  mean = compute_mean(mu, sigma)

  mode = compute_mode(mu, sigma)


  print tuple(mean)
  print tuple(mode)
