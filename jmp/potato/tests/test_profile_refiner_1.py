from __future__ import division

import numpy.random
from scitbx import matrix
from math import sqrt, pi, sin, cos, log, exp, ceil
from dials_scratch.jmp.potato.util.generate_simple import generate_simple
from dials_scratch.jmp.potato.util.generate_simple import generate_simple_binned
from dials_scratch.jmp.potato.util.simplex import SimpleSimplex
from dials_scratch.jmp.potato.model import compute_change_of_basis_operation
from dials_scratch.jmp.potato.parameterisation import SimpleMosaicityParameterisation
from dials_scratch.jmp.potato.profile_refiner_target import MaximumLikelihoodTarget
from dials.array_family import flex


def log_likelihood(params, s0, s2_list, xbar_list, ctot_list, Sobs_list):
  '''
  The log likelihood given the data

  '''

  # Construct covariance
  M = matrix.sqr((
    params[0], 0, 0,
    params[1], params[2], 0,
    params[3], params[4], params[5]))
  sigma = M*M.transpose()

  # Compute the loglikelihood
  lnL = 0
  for i in range(len(s2_list)):

    # Set stuff
    s2 = s2_list[i]
    xbar = xbar_list[i]
    ctot = ctot_list[i]
    Sobs = Sobs_list[i]

    # Get the change of basis operation
    R = compute_change_of_basis_operation(s0, s2)

    # Compute rotated sigma
    S = R*sigma*R.transpose()
    S11 = matrix.sqr((
      S[0], S[1],
      S[3], S[4]))
    S12 = matrix.col((S[2], S[5]))
    S21 = matrix.col((S[6], S[7])).transpose()
    S22 = S[8]

    # Compute rotated mu
    mu = R*s2
    assert abs(1 - mu.normalize().dot(matrix.col((0, 0, 1)))) < 1e-7
    mu1 = matrix.col((mu[0], mu[1]))
    mu2 = mu[2]

    # Compute conditional mean and covariance
    mubar = mu1 + S12*(1/S22)*(s0.length()-mu2)
    Sbar = S11 - S12*(1/S22)*S21

    # Compute the log likelihood
    try:
      Sobs = matrix.sqr(Sobs)
      B1 = log(S22)
      B2 = (1/S22)*(s0.length()-mu2)**2
      A1 = log(Sbar.determinant())*ctot
      A2 = ((Sbar.inverse())*(Sobs)).trace()
      A3 = ((Sbar.inverse())*(ctot*(xbar-mubar)*(xbar-mubar).transpose())).trace()
      lnL += -0.5*(A1+A2+A3+(B1+B2))
    except Exception:
      raise
      lnL += -1e15

  # print tuple(sigma), lnL
  return lnL


class Target(object):

  def __init__(self,
               s0,
               s2_list,
               xbar_list,
               ctot_list,
               Sobs_list):
    self.s0 = s0
    self.s2_list = s2_list
    self.xbar_list = xbar_list
    self.ctot_list = ctot_list
    self.Sobs_list = Sobs_list

  def target(self, params):
    lnL = log_likelihood(
      params,
      self.s0,
      self.s2_list,
      self.xbar_list,
      self.ctot_list,
      self.Sobs_list)
    score = -lnL
    return score


def tst_ideal():

  numpy.random.seed(100)

  # The beam vector
  s0 = matrix.col((0, 0, 1))

  # The covariance matrix
  sigma = matrix.sqr((
    1e-6, 0, 0,
    0, 2e-6, 0,
    0, 0, 3e-6))

  # The number of reflections
  N = 1000

  # Generate a load of reflections
  s2_list, ctot_list, xbar_list, Sobs_list = generate_simple(s0, sigma, N = N)

  # Starting values for simplex
  values = flex.double((
    sqrt(1e-6),
    0, sqrt(1e-6),
    0, 0, sqrt(1e-6)))
  offset = flex.double(
    [sqrt(1e-7)  for v in values])

  # Do the simplex optimization
  optimizer = SimpleSimplex(
    values,
    offset,
    Target(
      s0,
      s2_list,
      xbar_list,
      ctot_list,
      Sobs_list),
    2000)
  params = optimizer.get_solution()

  # Create the covariance matrix
  M = matrix.sqr((
    params[0], 0, 0,
    params[1], params[2], 0,
    params[3], params[4], params[5]))
  sigma = M*M.transpose()

  print sigma

  expected = matrix.sqr((
    9.96410142906e-07, -2.21919325401e-09, 1.82450917959e-11,
    -2.21919325401e-09, 1.98250217906e-06, -1.13761980971e-09,
    1.82450917959e-11, -1.13761980971e-09, 2.99278336529e-06))
  assert all(1e6*abs(a-b) < 1e-7 for a, b in zip(sigma, expected))

  print 'OK'

def tst_binned():

  numpy.random.seed(100)

  # The beam vector
  s0 = matrix.col((0, 0, 1))

  # The covariance matrix
  sigma = matrix.sqr((
    1e-6, 0, 0,
    0, 2e-6, 0,
    0, 0, 3e-6))

  # The number of reflections
  N = 1000

  # Generate a load of reflections
  s2_list, ctot_list, xbar_list, Sobs_list = generate_simple_binned(s0, sigma, N = N)

  # Starting values for simplex
  values = flex.double((
    sqrt(1e-6),
    0, sqrt(1e-6),
    0, 0, sqrt(1e-6)))
  offset = flex.double(
    [sqrt(1e-7)  for v in values])

  # Do the simplex optimization
  optimizer = SimpleSimplex(
    values,
    offset,
    Target(
      s0,
      s2_list,
      xbar_list,
      ctot_list,
      Sobs_list),
    2000)
  params = optimizer.get_solution()

  # Create the covariance matrix
  M = matrix.sqr((
    params[0], 0, 0,
    params[1], params[2], 0,
    params[3], params[4], params[5]))
  sigma = M*M.transpose()

  print sigma

  expected = matrix.sqr((
    1.05670907537e-06, -1.78247842812e-09, -8.19788355586e-09,
    -1.78247842812e-09, 2.08085121753e-06, 6.30026722264e-09,
    -8.19788355586e-09, 6.30026722264e-09, 3.15006720619e-06))

  assert all(1e6*abs(a-b) < 1e-7 for a, b in zip(sigma, expected))

  print 'OK'


def tst_ml_target_class():


  class SimplexTarget(object):

    def __init__(self,
                 s0,
                 s2_list,
                 ctot_list,
                 xbar_list,
                 Sobs_list):
      self.s0 = s0
      self.s2_list = s2_list
      self.xbar_list = xbar_list
      self.ctot_list = ctot_list
      self.Sobs_list = Sobs_list

    def target(self, params):

      parameterisation = SimpleMosaicityParameterisation(params)

      t = MaximumLikelihoodTarget(
        parameterisation,
        self.s0,
        self.s2_list,
        self.ctot_list,
        self.xbar_list,
        self.Sobs_list)

      lnL = t.log_likelihood()

      # print tuple(parameterisation.sigma()), lnL

      return -lnL

  numpy.random.seed(100)

  # The beam vector
  s0 = matrix.col((0, 0, 1))

  # The covariance matrix
  sigma = matrix.sqr((
    1e-6, 0, 0,
    0, 2e-6, 0,
    0, 0, 3e-6))

  # The number of reflections
  N = 1000

  # Generate a load of reflections
  s2_list, ctot_list, xbar_list, Sobs_list = generate_simple(s0, sigma, N = N)

  # Starting values for simplex
  values = flex.double((
    sqrt(1e-6),
    0, sqrt(1e-6),
    0, 0, sqrt(1e-6)))
  offset = flex.double(
    [sqrt(1e-7)  for v in values])

  # Do the simplex optimization
  optimizer = SimpleSimplex(
    values,
    offset,
    SimplexTarget(
      s0,
      s2_list,
      ctot_list,
      xbar_list,
      Sobs_list),
    2000)
  params = optimizer.get_solution()

  # Create the covariance matrix
  M = matrix.sqr((
    params[0], 0, 0,
    params[1], params[2], 0,
    params[3], params[4], params[5]))
  sigma = M*M.transpose()

  # print sigma

  expected = matrix.sqr((
    9.96410142906e-07, -2.21919325401e-09, 1.82450917959e-11,
    -2.21919325401e-09, 1.98250217906e-06, -1.13761980971e-09,
    1.82450917959e-11, -1.13761980971e-09, 2.99278336529e-06))
  assert all(1e6*abs(a-b) < 1e-7 for a, b in zip(sigma, expected))

  print 'OK'


if __name__ == '__main__':
  #tst_ideal()
  #tst_binned()
  tst_ml_target_class()
