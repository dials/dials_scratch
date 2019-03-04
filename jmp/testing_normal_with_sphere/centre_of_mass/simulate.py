from __future__ import print_function
from scitbx import matrix
from math import pi, cos, sin, exp, sqrt, atan2, log, tan, acos
from random import uniform
import json
from time import time

def normal_3d(x, mu, sigma):
  A = 1.0 / sqrt((2*pi)**3 * sigma.determinant())
  B = ((x - mu).transpose() * sigma.inverse() * (x - mu))[0]
  return A * exp(-0.5 * B)

def compute_integrate(mu, sigma, N=1000, M=1000):

  I1 = matrix.col((0, 0, 0))
  I2 = 0
  T0 = 0
  T1 = pi/2
  P0 = 0
  P1 = pi/2
  for j in range(N):
    for i in range(M):
      theta = T0 + j * (T1-T0) / N
      phi = P0 + i * (P1-P0) / N

      sin_theta = sin(theta)
      cos_theta = cos(theta)
      sin_phi = sin(phi)
      cos_phi = cos(phi)

      x = matrix.col((
        sin_theta*cos_phi,
        sin_theta*sin_phi,
        cos_theta))

      f = normal_3d(x, mu, sigma)

      I1 += x * f * sin_theta
      I2 += f * sin_theta

  return (I1 / I2).normalize()

def generate_random_sigma():
  return matrix.sqr((
    uniform(0.0001, 0.005), 0, 0,
    0, uniform(0.0001, 0.005), 0,
    0, 0, uniform(0.0001, 0.005)))

def generate_random_mu():
  scale = uniform(0.95, 1.05)
  return matrix.col((1, 1, 1)).normalize() * scale


def simulate(num_spots,
             num_integral,
             filename):

  data = {
    'mu' : [],
    'sigma' : [],
    'x' : []
  }

  for i in range(num_spots):

    sigma = generate_random_sigma()
    mu = generate_random_mu()

    st = time()
    x = compute_integrate(mu, sigma, num_integral, num_integral)
    ft = time()

    print(i, mu.length(), tuple(sigma)[::4], x.angle(mu), int(ft-st))

    data['mu'].append(tuple(mu))
    data['sigma'].append(tuple(sigma))
    data['x'].append(tuple(x))

  json.dump(data, open(filename, "w"), indent=2)


if __name__ == "__main__":

  simulate(
    num_spots = 1000,
    num_integral = 1000,
    filename = "spots.txt")

  # mu = matrix.col((1, 1, 1)).normalize() * 0.95
  # sigma = matrix.sqr((
  #   0.01, 0, 0,
  #   0, 0.02, 0,
  #   0, 0, 0.03))

  # sest = compute_integrate(mu, sigma)

  # theta = acos(sest[2])
  # phi = atan2(sest[1], sest[0])
  # print theta, phi

  # print scal.angle(sest)
