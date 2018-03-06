from __future__ import division
from dials.array_family import flex
from scitbx import matrix
from math import log, exp
from dials_scratch.jmp.stills.potato.parameterisation import MosaicityParameterisation


class MarginalDistribution(object):
  '''
  A class to compute useful stuff about the marginal distribution

  '''

  def __init__(self, S, dS, d2S):

    # Compute the marginal variance
    self.S = S[8]

    # Compute the marginal derivatives
    self.dS = flex.double(d[8] for d in dS)

    # Compute the marginal second derivatives
    self.d2S = flex.double([d2[8] for d2 in d2S])
    self.d2S.reshape(d2S.accessor())

  def sigma(self):
    return self.S

  def first_derivatives(self):
    return self.dS

  def second_derivatives(self):
    return self.d2S

class ConditionalDistribution(object):

  def __init__(self, S, dS, d2S):

    # Partition the covariance matrix
    S11 = matrix.sqr((
      S[0], S[1],
      S[3], S[4]))
    S12 = matrix.col((S[2], S[5]))
    S21 = matrix.col((S[6], S[7])).transpose()
    S22 = S[8]

    # Compute the marginal covariance matrix
    self.S = S11 - S12*(1/S22)*S21
    #return
    def compute_dS(dS):

      dS11 = matrix.sqr((
        dS[0], dS[1],
        dS[3], dS[4]))
      dS12 = matrix.col((dS[2], dS[5]))
      dS21 = matrix.col((dS[6], dS[7])).transpose()
      dS22 = dS[8]

      S22_inv = 1 / S22

      A = dS11
      B = S12*S22_inv*dS22*S22_inv*S21
      C = S12*S22_inv*dS21
      D = dS12*S22_inv*S21
      return A + B - (C + D)

    self.dS = [compute_dS(d) for d in dS]
    return
    def compute_d2S(dSi, dSj, d2S):

      dSi12 = matrix.col((dSi[2], dSi[5]))
      dSi21 = matrix.col((dSi[6], dSi[7])).transpose()
      dSi22 = dSi[8]

      dSj12 = matrix.col((dSj[2], dSj[5]))
      dSj21 = matrix.col((dSj[6], dSj[7])).transpose()
      dSj22 = dSj[8]

      d2S11 = matrix.sqr((
        d2S[0], d2S[1],
        d2S[3], d2S[4]))
      d2S12 = matrix.col((d2S[2], d2S[5]))
      d2S21 = matrix.col((d2S[6], d2S[7])).transpose()
      d2S22 = d2S[8]

      S22_inv = 1 / S22

      A = d2S11
      B = dSj12*S22_inv*dSi22*S22_inv*S21
      C = S12*S22_inv*dSj22*S22_inv*dSi22*S22_inv*S21

      D = S12*S22_inv*d2S22*S22_inv*S21
      E = S12*S22_inv*dSi22*S22_inv*dSj22*S22_inv*S21
      F = S12*S22_inv*dSi22*S22_inv*dSj21

      G = dSj12*S22_inv*dSi21
      H = S12*S22_inv*dSj22*S22_inv*dSi21
      I = S12*S22_inv*d2S21

      J = d2S12*S22_inv*S21
      K = dSi12*S22_inv*(dSj22)*S22_inv*S21
      L = dSi12*S22_inv*dSj21

      return A+B-C+D-E+F-G+H-I-J+K-L

    self.d2S = [[
      compute_d2S(dS[i], dS[j], d2S[i,j])
      for j in range(d2S.all()[1])
    ] for i in range(d2S.all()[0])]


  def sigma(self):
    return self.S

  def first_derivatives(self):
    return self.dS

  def second_derivatives(self):
    return self.d2S


def rotate_mat3_double(R, A):
  accessor = A.accessor()
  RAR = flex.mat3_double([R*matrix.sqr(a)*R.transpose() for a in A])
  RAR.reshape(accessor)
  return RAR

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


class ReflectionProfileModel(object):

  def __init__(self,
               model,
               s0,
               s2,
               ctot,
               Sobs):

    self.s0 = s0
    self.s2 = s2
    self.R = compute_change_of_basis_operation(s0, s2)
    self.ctot = ctot
    self.Sobs = Sobs

    self.model = model

    # Rotate the covariance matrix
    S = self.R*self.model.sigma()*self.R.transpose()

    # Rotate the first derivative matrices
    dS = rotate_mat3_double(self.R, self.model.first_derivatives())

    # Rotate the second derivative matrices
    d2S = rotate_mat3_double(self.R, self.model.second_derivatives())

    # Construct the marginal distribution
    self._marginal = MarginalDistribution(S, dS, d2S)

    # Construct the conditional distribution
    self._conditional = ConditionalDistribution(S, dS, d2S)

    self.S = S
    self.dS = dS

  def marginal(self):
    return self._marginal

  def conditional(self):
    return self._conditional

  def log_likelihood(self):

    S22 = self.marginal().sigma()
    S22_inv = 1 / S22
    Sbar = self.conditional().sigma()
    Sbar_inv = Sbar.inverse()
    Sbar_det = Sbar.determinant()

    d = self.s0.length()-self.s2.length()
    A = log(S22)
    B = S22_inv*d**2
    C = log(Sbar_det)*self.ctot
    D = (Sbar_inv * self.Sobs).trace()
    return -0.5 * (A + B + C + D)


  def first_derivatives(self):

    # Get info about marginal distribution
    S22 = self.marginal().sigma()
    dS22 = self.marginal().first_derivatives()
    S22_inv = 1 / S22

    # Get info about conditional distribution
    Sbar = self.conditional().sigma()
    dSbar = self.conditional().first_derivatives()
    Sbar_inv = Sbar.inverse()

    # The distance from the ewald sphere
    d = self.s0.length()-self.s2.length()

    # Compute the derivative wrt parameter i
    dL = flex.double()
    for i in range(len(dS22)):

      I = matrix.sqr((
        1, 0,
        0, 1))

      U = S22_inv*dS22[i]*(1 - S22_inv*d**2)
      V = (Sbar_inv*dSbar[i]*(self.ctot*I - Sbar_inv*self.Sobs)).trace()

      dL.append(-0.5*(U+V))

    # Return the derivative of the log likelihood
    return dL

  def second_derivatives(self):

    # Get info about marginal distribution
    S22 = self.marginal().sigma()
    dS22 = self.marginal().first_derivatives()
    d2S22 = self.marginal().second_derivatives()
    S22_inv = 1 / S22

    # Get info about conditional distribution
    Sbar = self.conditional().sigma()
    dSbar = self.conditional().first_derivatives()
    d2Sbar = self.conditional().second_derivatives()
    Sbar_inv = Sbar.inverse()

    # The distance from the ewald sphere
    d = self.s0.length()-self.s2.length()

    # Compute the derivative wrt parameter i
    d2L = flex.double(d2S22.accessor())
    for j in range(d2S22.all()[0]):
      for i in range(d2S22.all()[1]):

        I = matrix.sqr((
          1, 0,
          0, 1))

        A1 = S22_inv*d2S22[j,i]*(1 - S22_inv*d**2)
        A2 = S22_inv*dS22[j]*S22_inv*dS22[i]*(1 - 2*S22_inv*d**2)
        B1 = Sbar_inv * d2Sbar[j][i]*(self.ctot*I - Sbar_inv*self.Sobs)
        B2 = Sbar_inv * dSbar[j] * Sbar_inv * dSbar[i] * (self.ctot*I - 2*Sbar_inv*self.Sobs)
        U = A1-A2
        V = (B1-B2).trace()
        d2L[j,i] = -0.5*(U+V)

    # Return the second derivative of the log likelihood
    return d2L

  def fisher_information(self):

    # Get info about marginal distribution
    S22 = self.marginal().sigma()
    dS22 = self.marginal().first_derivatives()
    S22_inv = 1 / S22

    # Get info about conditional distribution
    Sbar = self.conditional().sigma()
    dSbar = self.conditional().first_derivatives()
    Sbar_inv = Sbar.inverse()

    # The distance from the ewald sphere
    d = self.s0.length()-self.s2.length()

    I = flex.double(flex.grid(len(dS22), len(dS22)))
    for j in range(len(dS22)):
      for i in range(len(dS22)):

        U = S22_inv*dS22[j]*S22_inv*dS22[i]
        V = (Sbar_inv*dSbar[j]*Sbar_inv*dSbar[i]*self.ctot).trace()
        I[j,i] = 0.5*(U+V)

    return I


class ReflectionProfileModelList(object):
  '''
  A class to hold all the reflection profile models

  '''
  def __init__(self,
               parameterisation,
               s0,
               s2_list,
               ctot_list,
               Sobs_list):
    '''
    Initialise

    '''

    # Save the parameterisation
    self._parameterisation = parameterisation

    # Construct the models
    self._models = []
    for s2, ctot, Sobs in zip(s2_list, ctot_list, Sobs_list):
      self._models.append(
        ReflectionProfileModel(parameterisation, s0, s2, ctot, Sobs))

  def __len__(self):
    '''
    The number of models

    '''
    return len(self._models)

  def parameterisation(self):
    '''
    Return the parameterisation

    '''
    return self._parameterisation

  def log_likelihood(self):
    '''
    The joint log likelihood

    '''
    lnL = 0
    for i in range(len(self)):
      lnL += self._models[i].log_likelihood()
    return lnL

  def first_derivatives(self):
    '''
    The joint first derivatives

    '''
    dL = 0
    for i in range(len(self)):
      dL += self._models[i].first_derivatives()
    return dL

  def second_derivatives(self):
    '''
    The joint second derivatives

    '''
    d2L = 0
    for i in range(len(self)):
      d2L += self._models[i].second_derivatives()
    return d2L

  def fisher_information(self):
    '''
    The joint fisher information

    '''
    I = 0
    for i in range(len(self)):
      I += self._models[i].fisher_information()
    return I
