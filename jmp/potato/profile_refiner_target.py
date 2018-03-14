from __future__ import division
from scitbx import matrix
from dials.array_family import flex
from dials_scratch.jmp.potato.model import compute_change_of_basis_operation
from math import log

class MarginalDistribution(object):
  '''
  A class to compute useful stuff about the marginal distribution

  '''

  def __init__(self, S, dS, d2S=None):

    # Compute the marginal variance
    self.S = S[8]

    # Compute the marginal derivatives
    self.dS = flex.double(d[8] for d in dS)

    # Compute the marginal second derivatives
    if d2S is not None:
      self.d2S = flex.double([d2[8] for d2 in d2S])
      self.d2S.reshape(d2S.accessor())
    else:
      self.d2S = None

  def sigma(self):
    '''
    Return the marginal sigma

    '''
    return self.S

  def first_derivatives(self):
    '''
    Return the marginal first derivatives

    '''
    return self.dS

  def second_derivatives(self):
    '''
    Return the maginal second derivatives

    '''
    return self.d2S


class ConditionalDistribution(object):
  '''
  A class to compute useful stuff about the conditional distribution

  '''

  def __init__(self, s0, s2, S, dS, d2S=None):

    self._S = S
    self._dS = dS
    self._d2S = d2S

    # Partition the covariance matrix
    S11 = matrix.sqr((
      S[0], S[1],
      S[3], S[4]))
    S12 = matrix.col((S[2], S[5]))
    S21 = matrix.col((S[6], S[7])).transpose()
    S22 = S[8]

    # The partitioned mean vector
    mu1 = matrix.col((0, 0))
    mu2 = s2.length()
    a = s0.length()

    # The epsilon
    self.epsilon = a - mu2

    # Compute the conditional mean
    self.mubar = mu1 + S12*(1/S22)*self.epsilon

    # Compute the conditional covariance matrix
    self.Sbar = S11 - S12*(1/S22)*S21

    # Set to None and compute on demand
    self.dSbar = None
    self.dmbar = None
    self.d2Sbar = None
    self.d2mbar = None

  def mean(self):
    '''
    Return the conditional mean

    '''
    return self.mubar

  def sigma(self):
    '''
    Return the conditional sigma

    '''
    return self.Sbar

  def first_derivatives_of_sigma(self):
    '''
    Return the marginal first derivatives

    '''
    def compute_dSbar(S, dS):

      S12 = matrix.col((S[2], S[5]))
      S21 = matrix.col((S[6], S[7])).transpose()
      S22 = S[8]

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

    if self.dSbar is None:
      self.dSbar = [compute_dSbar(self._S, d) for d in self._dS]

    return self.dSbar

  def first_derivatives_of_mean(self):
    '''
    Return the marginal first derivatives

    '''
    def compute_dmbar(S, dS):

      S12 = matrix.col((S[2], S[5]))
      S21 = matrix.col((S[6], S[7])).transpose()
      S22 = S[8]

      dS11 = matrix.sqr((
        dS[0], dS[1],
        dS[3], dS[4]))
      dS12 = matrix.col((dS[2], dS[5]))
      dS21 = matrix.col((dS[6], dS[7])).transpose()
      dS22 = dS[8]

      S22_inv = 1 / S22

      A = dS12*S22_inv
      B = -S12*S22_inv*dS22*S22_inv
      return (A + B)*self.epsilon

    if self.dmbar is None:
      self.dmbar = [compute_dmbar(self._S, d) for d in self._dS]

    return self.dmbar

  def second_derivatives_of_mean(self):
    '''
    Return the maginal second derivatives

    '''
    def compute_d2mbar(S, dSi, dSj, d2S):

      S12 = matrix.col((S[2], S[5]))
      S21 = matrix.col((S[6], S[7])).transpose()
      S22 = S[8]

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

      A = d2S12*S22_inv
      B = dSi12*S22_inv*dSj22*S22_inv
      C = dSj12*S22_inv*dSi22*S22_inv

      D = S12*S22_inv*dSj22*S22_inv*dSi22*S22_inv
      E = S12*S22_inv*d2S22*S22_inv
      F = S12*S22_inv*dSi22*S22_inv*dSj22*S22_inv

      return (A-B-C+D-E+F)*self.epsilon

    if self.d2mbar is None:
      self.d2mbar = [[
        compute_d2mbar(self._S, self._dS[i], self._dS[j], self._d2S[i,j])
        for j in range(self._d2S.all()[1])
      ] for i in range(self._d2S.all()[0])]

    return self.d2mbar

  def second_derivatives_of_sigma(self):
    '''
    Return the maginal second derivatives

    '''
    def compute_d2S(S, dSi, dSj, d2S):

      S12 = matrix.col((S[2], S[5]))
      S21 = matrix.col((S[6], S[7])).transpose()
      S22 = S[8]

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

    if self.d2Sbar is None:
      self.d2Sbar = [[
        compute_d2S(self._S, self._dS[i], self._dS[j], self._d2S[i,j])
        for j in range(self._d2S.all()[1])
      ] for i in range(self._d2S.all()[0])]

    return self.d2Sbar



def rotate_mat3_double(R, A):
  '''
  Helper function to rotate a flex.mat3_double array of matrices

  '''
  accessor = A.accessor()
  RAR = flex.mat3_double([R*matrix.sqr(a)*R.transpose() for a in A])
  RAR.reshape(accessor)
  return RAR


class ReflectionData(object):

  def __init__(self, model, s0, s2, ctot, mobs, sobs):

    # Save stuff
    self.model = model
    self.s0 = s0
    self.s2 = s2
    self.ctot = ctot
    self.mobs = mobs
    self.sobs = sobs

    # Compute the change of basis
    self.R = compute_change_of_basis_operation(s0, s2)

    # Rotate the covariance matrix
    self.S = self.R*model.sigma()*self.R.transpose()

    # Rotate the first derivative matrices
    self.dS = rotate_mat3_double(self.R, model.first_derivatives())

    # Rotate the first derivative matrices
    self.d2S = rotate_mat3_double(self.R, model.second_derivatives())

    # Construct the marginal distribution
    self.marginal = MarginalDistribution(self.S, self.dS, self.d2S)

    # Construct the conditional distribution
    self.conditional = ConditionalDistribution(s0, s2, self.S, self.dS, self.d2S)

  def log_likelihood(self):
    '''
    Compute the log likelihood for the reflection

    '''

    # Get data
    s0 = self.s0
    s2 = self.s2
    ctot = self.ctot
    mobs = self.mobs
    Sobs = self.sobs

    # Get info about the marginal
    S22 = data.marginal.sigma()
    S22_inv = 1 / S22

    # Get info about the conditional
    Sbar = self.conditional.sigma()
    mubar = self.conditional.mean()
    Sbar_inv = Sbar.inverse()
    Sbar_det = Sbar.determinant()

    # Compute the marginal likelihood
    m_d = s0.length() - s2.length()
    m_lnL = log(S22) + S22_inv*m_d**2

    # Compute the conditional likelihood
    c_d = mobs - mubar
    c_lnL = log(Sbar_det)*ctot + (Sbar_inv * (Sobs + ctot*c_d*c_d.transpose())).trace()

    # Return the joint likelihood
    return -0.5 * (m_lnL + c_lnL)


  def first_derivatives(self):
    '''
    Compute the first derivatives

    '''
    # Get data
    s0 = self.s0
    s2 = self.s2
    ctot = self.ctot
    mobs = self.mobs
    Sobs = self.sobs

    # Get info about marginal distribution
    S22 = self.marginal.sigma()
    dS22 = self.marginal.first_derivatives()
    S22_inv = 1 / S22

    # Get info about conditional distribution
    Sbar = self.conditional.sigma()
    mubar = self.conditional.mean()
    dSbar = self.conditional.first_derivatives_of_sigma()
    dmbar = self.conditional.first_derivatives_of_mean()
    Sbar_inv = Sbar.inverse()

    # The distance from the ewald sphere
    m_d = s0.length() - s2.length()
    c_d = mobs - mubar

    # Compute the derivative wrt parameter i
    dL = flex.double()
    for i in range(len(dS22)):

      I = matrix.sqr((
        1, 0,
        0, 1))

      U = S22_inv*dS22[i]*(1 - S22_inv*m_d**2)
      V = (Sbar_inv*dSbar[i]*(ctot*I - Sbar_inv*(Sobs+ctot*c_d*c_d.transpose()))).trace()
      W = (-2*Sbar_inv*(ctot*c_d*dmbar[i].transpose())).trace()

      dL.append(-0.5*(U+V+W))

    # Return the derivative of the log likelihood
    return dL

  def second_derivatives(self):
    '''
    Compute the second derivatives

    '''
    # Get data
    s0 = self.s0
    s2 = self.s2
    ctot = self.ctot
    mobs = self.mobs
    Sobs = self.sobs

    # Get info about marginal distribution
    S22 = self.marginal.sigma()
    dS22 = self.marginal.first_derivatives()
    d2S22 = self.marginal.second_derivatives()
    S22_inv = 1 / S22

    # Get info about conditional distribution
    Sbar = self.conditional.sigma()
    dSbar = self.conditional.first_derivatives_of_sigma()
    d2Sbar = self.conditional.second_derivatives_of_sigma()
    mubar = self.conditional.mean()
    dmbar = self.conditional.first_derivatives_of_mean()
    d2mbar = self.conditional.second_derivatives_of_mean()
    Sbar_inv = Sbar.inverse()

    # The distance from the ewald sphere
    m_d = s0.length() - s2.length()
    c_d = mobs - mubar

    # Compute the derivative wrt parameter i j
    d2L = flex.double(d2S22.accessor())
    for j in range(d2S22.all()[0]):
      for i in range(d2S22.all()[1]):

        I = matrix.sqr((
          1, 0,
          0, 1))

        A1 = S22_inv*d2S22[j,i]*(1 - S22_inv*m_d**2)
        A2 = S22_inv*dS22[j]*S22_inv*dS22[i]*(1 - 2*S22_inv*m_d**2)
        B1 = Sbar_inv * d2Sbar[j][i]*(ctot*I - Sbar_inv*(Sobs+ctot*c_d*c_d.transpose()))
        B2 = Sbar_inv * dSbar[j] * Sbar_inv * dSbar[i] * (ctot*I - 2*Sbar_inv*(Sobs+ctot*c_d*c_d.transpose()))
        B3 = Sbar_inv * dSbar[i] * 2 * Sbar_inv * ctot * c_d*dmbar[j].transpose()
        B4 = Sbar_inv * dSbar[j] * 2 * Sbar_inv * ctot * c_d*dmbar[i].transpose()
        B5 = 2*Sbar_inv * ctot*dmbar[j]*dmbar[i].transpose()
        B6 = 2*Sbar_inv * ctot*c_d*d2mbar[j][i].transpose()
        U = A1-A2
        V = (B1-B2+B3+B4+B5-B6).trace()
        d2L[j,i] = -0.5*(U+V)

    # Return the second derivative of the log likelihood
    return d2L

  def fisher_information(self):
    '''
    Compute the fisher information

    '''

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

    # Compute the fisher information wrt parameter i j
    I = flex.double(flex.grid(len(dS22), len(dS22)))
    for j in range(len(dS22)):
      for i in range(len(dS22)):

        U = S22_inv*dS22[j]*S22_inv*dS22[i]
        V = (Sbar_inv*dSbar[j]*Sbar_inv*dSbar[i]*self.ctot).trace()
        W = c*(dmbar[i].transpose()*Sbar_inv*dmbar[j])[0]
        I[j,i] = 0.5*(U+V) + W

    return I


class MaximumLikelihoodTarget(object):

  def __init__(self,
               model,
               s0,
               s2_list,
               ctot_list,
               mobs_list,
               sobs_list):

    # Check input
    assert len(s2_list) == len(ctot_list)
    assert len(s2_list) == len(mobs_list)
    assert len(s2_list) == len(sobs_list)

    # Compute the change of basis for each reflection
    self.data = []
    for i in range(len(s2_list)):
      self.data.append(ReflectionData(
        model,
        s0,
        s2_list[i],
        ctot_list[i],
        mobs_list[i],
        matrix.sqr(sobs_list[i])))

  def log_likelihood(self):
    '''
    The joint log likelihood

    '''
    lnL = 0
    for i in range(len(self.data)):
      lnL += self.data[i].log_likelihood()
    return lnL

  def first_derivatives(self):
    '''
    The joint first derivatives

    '''
    dL = 0
    for i in range(len(self.data)):
      dL += self.data[i].first_derivatives()
    return dL

  def second_derivatives(self):
    '''
    The joint second derivatives

    '''
    d2L = 0
    for i in range(len(self.data)):
      d2L += self.data[i].second_derivatives()
    return d2L

  def fisher_information(self):
    '''
    The joint fisher information

    '''
    I = 0
    for i in range(len(self.data)):
      I += self.data[i].fisher_information()
    return I
