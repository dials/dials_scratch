from __future__ import division
from dials.array_family import flex
from scitbx import matrix

class MosaicityParameterisation(object):

  def __init__(self, params):
    self.params = params

  def parameters(self):
    return self.params

  def sigma(self):
    M = matrix.sqr((
      self.params[0], 0, 0,
      self.params[1], self.params[2], 0,
      self.params[3], self.params[4], self.params[5]))
    return M*M.transpose()

  def first_derivatives(self):
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


class ProfileModel(object):

  def __init__(self,
               parameterisation):
    self.parameterisation = parameterisation

  def parameters(self):
    return self.parameterisation.parameters()

  def sigma(self):
    return self.parameterisation.sigma()

  def first_derivatives(self):
    return self.parameterisation.first_derivatives()

  def second_derivatives(self):
    return self.parameterisation.second_derivatives()


class MarginalDistribution(object):
  
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

    def compute_d2S(dSi, dSj, d2S):
      
      dSi11 = matrix.sqr((
        dSi[0], dSi[1],
        dSi[3], dSi[4]))
      dSi12 = matrix.col((dSi[2], dSi[5]))
      dSi21 = matrix.col((dSi[6], dSi[7])).transpose()
      dSi22 = dSi[8]
      
      dSj11 = matrix.sqr((
        dSj[0], dSj[1],
        dSj[3], dSj[4]))
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

class ReflectionProfileModel(object):

  def __init__(self, 
               model,
               R,
               mu2,
               r,
               ctot,
               Sobs):

    self.R = R
    self.mu2 = mu2
    self.r = r
    self.ctot = ctot
    self.Sobs = Sobs

    self.model = model

    # Rotate the covariance matrix
    S = R*self.model.sigma()*R.transpose()
    
    # Rotate the first derivative matrices
    dS = rotate_mat3_double(R, self.model.first_derivatives())
   
    # Partition the rotated second derivative matrices
    d2S = rotate_mat3_double(R, self.model.second_derivatives())

    # Construct the marginal distribution
    self._marginal = MarginalDistribution(S, dS, d2S)
    
    # Construct the conditional distribution
    self._conditional = ConditionalDistribution(S, dS, d2S)

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

    d = self.r-self.mu2
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
   
    # Compute the derivative wrt parameter i
    dL = flex.double()
    for i in range(len(dS22)):
      
      I = matrix.sqr((
        1, 0,
        0, 1))

      U = S22_inv*dS22[i]*(1 - S22_inv*(self.r - self.mu2)**2)
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
   
    # Compute the derivative wrt parameter i
    d2L = flex.double(d2S22.accessor())
    for j in range(d2S22.all()[0]):
      for i in range(d2S22.all()[1]):

        I = matrix.sqr((
          1, 0,
          0, 1))

        d = self.r - self.mu2
        A1 = S22_inv*d2S22[j,i]*(1 - S22_inv*d**2)
        A2 = S22_inv*dS22[j]*S22_inv*dS22[i]*(1 - 2*S22_inv*d**2)
        B1 = Sbar_inv * d2Sbar[j][i]*(self.ctot*I - Sbar_inv*self.Sobs)
        B2 = Sbar_inv * dSbar[j] * Sbar_inv * dSbar[i] * (self.ctot*I - 2*Sbar_inv*self.Sobs)
        U = A1-A2
        V = (B1-B2).trace()
        d2L[j,i] = -0.5*(U+V)
   
    # Return the second derivative of the log likelihood
    return d2L

if __name__ == '__main__':

  p = MosaicityParameterisation((1, 0.1, 2, 0.2, 0.3, 3))

  print p.sigma()
  print list(p.first_derivatives())
  print list(p.second_derivatives())
