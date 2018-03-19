from __future__ import division
import scitbx.linalg
from dials_scratch.jmp.potato.parameterisation import SimpleMosaicityParameterisation
from dials.array_family import flex
from scitbx import matrix


class SimpleMosaicityModel(object):
  '''
  A class to hold a simple 6 parameter mosaicity model

  '''

  def __init__(self, sigma):
    '''
    Initialise the with mosaicity matrix

    '''
    self._sigma = sigma

  def parameterisation(self):
    '''
    Return the parameterisation

    '''

    # Construct triangular matrix
    LL = flex.double()
    for j in range(3):
      for i in range(j+1):
        LL.append(self._sigma[j*3+i])

    # Do the cholesky decomposition
    ll = scitbx.linalg.l_l_transpose_cholesky_decomposition_in_place(LL)

    # Setup the parameters
    params = (
      LL[0],
      LL[1], LL[2],
      LL[3], LL[4], LL[5])

    # Return the parameterisation
    return SimpleMosaicityParameterisation(params)

  def compose(self, parameterisation):
    '''
    Compose the model from the parameterisation

    '''
    self._sigma = parameterisation.sigma()

  def sigma(self):
    '''
    Get the covariance matrix

    '''
    return self._sigma


def compute_change_of_basis_operation(s0, s2):
  '''
  Compute the change of basis operation that puts the s2 vector along the z axis

  '''
  e1 = s2.cross(s0).normalize()
  e2 = s2.cross(e1).normalize()
  e3 = s2.normalize()
  R = matrix.sqr(
    e1.elems +
    e2.elems +
    e3.elems)
  return R


# class ReflectionModel(object):

#   def __init__(self, model, s2):

#     # Save some stuff
#     self.model = model
#     self.s0 = s0
#     self.s2 = s2

#     # Compute the change of basis
#     self.R = compute_change_of_basis_operation(s0, s2)

#     # Rotate the covariance matrix
#     self.S = self.R*self.model.sigma()*self.R.transpose()

#     # Construct the marginal distribution
#     self._marginal = MarginalDistribution(self.S)

#     # Construct the conditional distribution
#     self._conditional = ConditionalDistribution(self.S, s0.length(), s2.length())

#   def marginal(self):
#     '''
#     Return the marginal

#     '''
#     return self._marginal

#   def conditional(self):
#     '''
#     Return the conditional

#     '''
#     return self._conditional

#   def s1(self):
#     '''
#     Compute the diffracted beam vector

#     '''
#     mubar = self._conditional.mean()
#     v = matrix.col((
#       mubar[0],
#       mubar[1],
#       self.s0.length())).normalize() * self.s0.length()
#     return self.R.transpose() * v

