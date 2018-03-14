from __future__ import division
from dials.array_family import flex
from scitbx import matrix

class SimpleMosaicityParameterisation(object):
  '''
  A simple mosaicity parameterisation that uses 6 parameters to describe a
  multivariate normal reciprocal lattice profile. Sigma is enforced as positive
  definite by parameterising using the cholesky decomposition.

  M = | b1 0  0  |
      | b2 b3 0  |
      | b4 b5 b6 |

  S = M*M^T

  '''

  def __init__(self, params):
    '''
    Initialise with the parameters

    '''
    self.params = params

  def parameters(self):
    '''
    Return the parameters

    '''
    return self.params

  def sigma(self):
    '''
    Compute the covariance matrix of the MVN from the parameters

    '''
    M = matrix.sqr((
      self.params[0], 0, 0,
      self.params[1], self.params[2], 0,
      self.params[3], self.params[4], self.params[5]))
    return M*M.transpose()

  def first_derivatives(self):
    '''
    Compute the first derivatives of Sigma w.r.t the parameters

    '''
    b1, b2, b3, b4, b5, b6 = self.params

    dSdb1 = (
      2*b1,b2,b4,
      b2,0,0,
      b4,0,0)

    dSdb2 = (
      0,b1,0,
      b1,2*b2,b4,
      0,b4,0)

    dSdb3 = (
      0,0,0,
      0,2*b3,b5,
      0,b5,0)

    dSdb4 = (
      0,0,b1,
      0,0,b2,
      b1,b2,2*b4)

    dSdb5 = (
      0,0,0,
      0,0,b3,
      0,b3,2*b5)

    dSdb6 = (
      0,0,0,
      0,0,0,
      0,0,2*b6)

    return flex.mat3_double([dSdb1, dSdb2, dSdb3, dSdb4, dSdb5, dSdb6])

  def second_derivatives(self):
    '''
    Compute the second derivatives of Sigma w.r.t the parameters

    '''
    b1, b2, b3, b4, b5, b6 = self.params

    zero = (
      0, 0, 0,
      0, 0, 0,
      0, 0, 0)

    d11 = (
      2, 0, 0,
      0, 0, 0,
      0, 0, 0)
    d12 = (
      0, 1, 0,
      1, 0, 0,
      0, 0, 0)
    d13 = zero
    d14 = (
      0, 0, 1,
      0, 0, 0,
      1, 0, 0)
    d15 = zero
    d16 = zero

    d21 = d12
    d22 = (
      0, 0, 0,
      0, 2, 0,
      0, 0, 0)
    d23 = zero
    d24 = (
      0, 0, 0,
      0, 0, 1,
      0, 1, 0)
    d25 = zero
    d26 = zero

    d31 = zero
    d32 = zero
    d33 = (
      0, 0, 0,
      0, 2, 0,
      0, 0, 0)
    d34 = zero
    d35 = (
      0, 0, 0,
      0, 0, 1,
      0, 1, 0)
    d36 = zero

    d41 = d14
    d42 = d24
    d43 = zero
    d44 = (
      0, 0, 0,
      0, 0, 0,
      0, 0, 2)
    d45 = zero
    d46 = zero

    d51 = zero
    d52 = zero
    d53 = d35
    d54 = zero
    d55 = (
      0, 0, 0,
      0, 0, 0,
      0, 0, 2)
    d56 = zero

    d61 = zero
    d62 = zero
    d63 = zero
    d64 = zero
    d65 = zero
    d66 = (
      0, 0, 0,
      0, 0, 0,
      0, 0, 2)

    d2 = flex.mat3_double([
      d11, d12, d13, d14, d15, d16,
      d21, d22, d23, d24, d25, d26,
      d31, d32, d33, d34, d35, d36,
      d41, d42, d43, d44, d45, d46,
      d51, d52, d53, d54, d55, d56,
      d61, d62, d63, d64, d65, d66])
    d2.reshape(flex.grid(6,6))

    return d2
