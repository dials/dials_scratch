from __future__ import division
from dials.array_family import flex
from scitbx import matrix
from scitbx import linalg
from math import log, exp

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
    # M = matrix.sqr((
    #   exp(self.params[0]), 0, 0,
    #   self.params[1], exp(self.params[2]), 0,
    #   self.params[3], self.params[4], exp(self.params[5])))
    return M*M.transpose()

  def first_derivatives(self):
    b1, b2, b3, b4, b5, b6 = self.params

    # dSdb1 = (
    #   2*exp(2*b1),b2*exp(b1),b4*exp(b1),
    #   b2*exp(b1),0,0,
    #   b4*exp(b1),0,0)

    # dSdb2 = (
    #   0,exp(b1),0,
    #   exp(b1),2*b2,b4,
    #   0,b4,0)

    # dSdb3 = (
    #   0,0,0,
    #   0,2*exp(2*b3),b5*exp(b3),
    #   0,b5*exp(b3),0)

    # dSdb4 = (
    #   0,0,exp(b1),
    #   0,0,b2,
    #   exp(b1),b2,2*b4)

    # dSdb5 = (
    #   0,0,0,
    #   0,0,exp(b3),
    #   0,exp(b3),2*b5)

    # dSdb6 = (
    #   0,0,0,
    #   0,0,0,
    #   0,0,2*exp(2*b6))

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

    b1, b2, b3, b4, b5, b6 = self.params

    zero = (
      0, 0, 0,
      0, 0, 0,
      0, 0, 0)

    # d11 = (
    #   4*exp(2*b1), b2*exp(b1), b4*exp(b1),
    #   b2*exp(b1), 0, 0,
    #   b4*exp(b1), 0, 0)
    # d12 = (
    #   0, exp(b1), 0,
    #   exp(b1), 0, 0,
    #   0, 0, 0)
    # d13 = zero
    # d14 = (
    #   0, 0, exp(b1),
    #   0, 0, 0,
    #   exp(b1), 0, 0)
    # d15 = zero
    # d16 = zero

    # d21 = d12
    # d22 = (
    #   0, 0, 0,
    #   0, 2, 0,
    #   0, 0, 0)
    # d23 = zero
    # d24 = (
    #   0, 0, 0,
    #   0, 0, 1,
    #   0, 1, 0)
    # d25 = zero
    # d26 = zero

    # d31 = zero
    # d32 = zero
    # d33 = (
    #   0, 0, 0,
    #   0, 4*exp(b3), b5*exp(b3),
    #   0, b5*exp(b3), 0)
    # d34 = zero
    # d35 = (
    #   0, 0, 0,
    #   0, 0, exp(b3),
    #   0, exp(b3), 0)
    # d36 = zero

    # d41 = d14
    # d42 = d24
    # d43 = zero
    # d44 = (
    #   0, 0, 0,
    #   0, 0, 0,
    #   0, 0, 2)
    # d45 = zero
    # d46 = zero

    # d51 = zero
    # d52 = zero
    # d53 = d35
    # d54 = zero
    # d55 = (
    #   0, 0, 0,
    #   0, 0, 0,
    #   0, 0, 2)
    # d56 = zero

    # d61 = zero
    # d62 = zero
    # d63 = zero
    # d64 = zero
    # d65 = zero
    # d66 = (
    #   0, 0, 0,
    #   0, 0, 0,
    #   0, 0, 4*exp(2*b6))

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
    #return
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

        # S12 = matrix.col((self.S[2], self.S[5]))
        # dSi12 = matrix.col((self.dS[i][2],self.dS[i][5]))
        # dSj12 = matrix.col((self.dS[j][2],self.dS[j][5]))
        # dSi22 = self.dS[i][8]
        # dSj22 = self.dS[j][8]

        # dmbari = dSi12*S22_inv*d - S12*S22_inv*dSi22*S22_inv*d
        # dmbarj = dSi12*S22_inv*d - S12*S22_inv*dSj22*S22_inv*d
        #W = (dmbari.transpose()*Sbar_inv*dmbarj)[0]

        U = S22_inv*dS22[j]*S22_inv*dS22[i]
        V = (Sbar_inv*dSbar[j]*Sbar_inv*dSbar[i]*self.ctot).trace()
        I[j,i] = 0.5*(U+V)

    return I

  
def line_search(func, x, p, tau=0.5, delta=1.0, tolerance=1e-7):
  '''
  Perform a line search
  :param func The function to minimize
  :param x The initial position
  :param p The direction to search
  :param tau: The backtracking parameter
  :param delta: The initial step
  :param tolerance: The algorithm tolerance
  :return: The amount to move in the given direction

  '''
  fa = func(x)
  min_delta = min(tolerance, tolerance / p.length())
  while delta > min_delta:
    fb = func(x + delta*p)
    if fb >= fa:
      return delta
    delta *= tau
  return 0


def gradient_descent(f, df, x0, max_iter=1000, tolerance=1e-10):
  '''
  Find the minimum using gradient descent and a line search
  :param f The function to minimize
  :param df The function to compute derivatives 
  :param x0 The initial position
  :param max_iter: The maximum number of iterations
  :param tolerance: The algorithm tolerance
  :return: The amount to move in the given direction

  '''
  delta = 0.5
  for it in range(max_iter):
    p = -matrix.col(df(x0))
    delta = line_search(f, x0, p, delta=min(1.0, delta*2), tolerance=tolerance)
    assert delta > 0
    x = x0 + delta*p
    assert f(x) <= f(x0)
    if (x - x0).length() < tolerance:
      break
    x0 = x
  return x


class FisherScoringMaximumLikelihoodBase(object):
  '''
  A class to solve maximum likelihood equations using fisher scoring

  '''

  def __init__(self, x0, max_iter=1000, tolerance=1e-7):
    '''
    Configure the algorithm
    
    :param x0: The initial parameter estimates
    :param max_iter: The maximum number of iterations
    :param tolerance: The parameter tolerance

    '''
    self.x0 = x0
    self.max_iter = max_iter
    self.tolerance = tolerance

  def solve(self):
    '''
    Find the maximum likelihood estimate

    '''
    x0 = self.x0
   
    # Loop through the maximum number of iterations
    for it in range(self.max_iter):
      
      # Compute the derivative and fisher information at x0
      S, I = self.score_and_fisher_information(x0)
     
      # Solve the update equation to get direction
      p = matrix.col(self.solve_update_equation(S, I))

      # Perform a line search to ensure that each step results in an increase the
      # in log likelihood. In the rare case where the update does not result in an
      # increase in the likelihood (only observed for absurdly small samples 
      # (e.g. 2 reflections) do an iteration of gradient descent 
      delta = self.line_search(x0, p) 
      if delta > 0:
        x = x0 + delta*p
      else:
        assert False
        x = self.gradient_search(x0) 

      # Call an update
      self.callback(x)

      # Break the loop if the parameters change less than the tolerance
      if (x - x0).length() < self.tolerance:
        break

      # Update the parameter
      x0 = x

    # Save the parameters
    self.num_iter = it+1
    self.parameters = x

  def solve_update_equation(self, D, I):
    '''
    Solve the update equation using cholesky decomposition
    :param D: The first derivatives
    :param I: The fisher information
    :return: The parameter delta

    '''
    
    # Construct triangular matrix
    LL = flex.double()
    for j in range(6):
      for i in range(j+1):
        LL.append(I[j*6+i])

    # Perform the decomposition
    ll = linalg.l_l_transpose_cholesky_decomposition_in_place(LL)
    p = flex.double(D)
    return ll.solve(p)

  def line_search(self,  x, p, tau=0.5, delta=1.0, tolerance=1e-7):
    '''
    Perform a line search
    :param x The initial position
    :param p The direction to search
    :return: The amount to move in the given direction

    '''
    def f(x):
      return -self.log_likelihood(x)

    return line_search(f, x, -p, tolerance=self.tolerance)

  def gradient_descent(self, x0):
    '''
    Find the minimum using gradient descent and a line search
    :param x0 The initial position
    :return: The amount to move in the given direction

    '''
    def f(x):
      return -self.log_likelihood(x)

    def df(x):
      return -self.score(x)

    return gradient_descent(f, df, x0, max_iter=1, tolerance=self.tolerance)


class FisherScoringMaximumLikelihood(FisherScoringMaximumLikelihoodBase):
  '''
  A class to solve the maximum likelihood equations

  '''
  def __init__(self, 
               x0, 
               s0, 
               s2_list, 
               ctot_list, 
               Sobs_list,
               max_iter=1000,
               tolerance=1e-7):
    '''
    Initialise the algorithm:

    '''
    
    # Initialise the super class
    super(FisherScoringMaximumLikelihood, self).__init__(
      x0,
      max_iter=max_iter,
      tolerance=tolerance)

    # Save some stuff
    self.s0 = s0
    self.s2_list = s2_list
    self.ctot_list = ctot_list
    self.Sobs_list = Sobs_list

  def log_likelihood(self, x):
    '''
    :param x: The parameter estimate
    :return: The log likelihood at x

    '''
    parameterisation = MosaicityParameterisation(x)
    profile_model = ProfileModel(parameterisation)
    lnL = 0
    for i in range(len(self.s2_list)):
      s2 = self.s2_list[i]
      ctot = self.ctot_list[i]
      Sobs = self.Sobs_list[i]
      r = ReflectionProfileModel(profile_model, self.s0, s2, ctot, Sobs)
      lnL += r.log_likelihood()
    return lnL

  def score(self, x):
    '''
    :param x: The parameter estimate
    :return: The score at x

    '''
    parameterisation = MosaicityParameterisation(x)
    profile_model = ProfileModel(parameterisation)
    S = 0
    for i in range(len(self.s2_list)):
      s2 = self.s2_list[i]
      ctot = self.ctot_list[i]
      Sobs = self.Sobs_list[i]
      r = ReflectionProfileModel(profile_model, self.s0, s2, ctot, Sobs)
      S += r.first_derivatives()
    return S

  def score_and_fisher_information(self, x):
    '''
    :param x: The parameter estimate
    :return: The score and fisher information at x

    '''
    parameterisation = MosaicityParameterisation(x)
    profile_model = ProfileModel(parameterisation)
    S = 0
    I = 0
    for i in range(len(self.s2_list)):
      s2 = self.s2_list[i]
      ctot = self.ctot_list[i]
      Sobs = self.Sobs_list[i]
      r = ReflectionProfileModel(profile_model, self.s0, s2, ctot, Sobs)
      S += r.first_derivatives()
      I += r.fisher_information()
    return S, I

  def callback(self, x):
    '''
    Handle and update in parameter values

    '''
    parameterisation = MosaicityParameterisation(x)
    profile_model = ProfileModel(parameterisation)
    sigma = profile_model.sigma()
    lnL = self.log_likelihood(x)
    format_string = "( %.3g, %.3g, %.3g, %.3g, %.3g, %.3g, %.3g, %.3g, %.3g ): L = %f" 
    print format_string % (tuple(sigma) + (lnL,))


class ProfileRefiner(object):

  def __init__(self, s0, s2_list, ctot_list, Sobs_list):
    self.s0 = s0
    self.s2_list = s2_list
    self.ctot_list = ctot_list
    self.Sobs_list = Sobs_list

  def refine(self):

    # Set the initial parameter values
    x0 = matrix.col((1, 0, 1, 0, 0, 1))

    # Initialise the algorithm
    ml = FisherScoringMaximumLikelihood(
      x0, 
      self.s0, 
      self.s2_list, 
      self.ctot_list,
      self.Sobs_list)
    
    # Solve the maximum likelihood equations
    ml.solve()

    # Get the parameters
    self.parameters = ml.parameters
    # ml.log_likelihood()
    

def estimate_parameters(s0, s2_list, ctot_list, Sobs_list, verbose=False):

  refiner = ProfileRefiner(
    s0,
    s2_list,
    ctot_list,
    Sobs_list)

  refiner.refine()

  return refiner.parameters

  # def f(x):
  #   parameterisation = MosaicityParameterisation(x)
  #   profile_model = ProfileModel(parameterisation)
  #   L = 0
  #   for i in range(len(s2_list)):
  #     s2 = s2_list[i]
  #     ctot = ctot_list[i]
  #     Sobs = Sobs_list[i]
  #     r = ReflectionProfileModel(profile_model, s0, s2, ctot, Sobs)
  #     L += r.log_likelihood()
  #   return L

  # def df(x):
  #   parameterisation = MosaicityParameterisation(x)
  #   profile_model = ProfileModel(parameterisation)
  #   dL = 0
  #   for i in range(len(s2_list)):
  #     s2 = s2_list[i]
  #     ctot = ctot_list[i]
  #     Sobs = Sobs_list[i]
  #     r = ReflectionProfileModel(profile_model, s0, s2, ctot, Sobs)
  #     dL += r.first_derivatives()
  #   return dL
  
  # def d2f(x):
  #   parameterisation = MosaicityParameterisation(x)
  #   profile_model = ProfileModel(parameterisation)
  #   d2L = 0
  #   for i in range(len(s2_list)):
  #     s2 = s2_list[i]
  #     ctot = ctot_list[i]
  #     Sobs = Sobs_list[i]
  #     r = ReflectionProfileModel(profile_model, s0, s2, ctot, Sobs)
  #     d2L += r.fisher_information()
  #   return d2L

  # return fisher_scoring_maximum_likelihood(f, df, d2f, x0, tolerance=1e-10, max_iter=1000)











def gradient_descent_old(func, dL, params):
  delta = line_search(func, params, matrix.col(-dL))
  p = params+delta*matrix.col(-dL)
  print delta, func(params), func(p)
  assert func(params) >= func(p)
  return p


def estimate_parameters2(s0, s2_list, ctot_list, Sobs_list, verbose=False):

  params_old = matrix.col((1, 0, 1, 0, 0, 1))
  while True:

    parameterisation = MosaicityParameterisation(params_old)

    profile_model = ProfileModel(parameterisation)

    reflection_models = []
    for i in range(len(s2_list)):
      s2 = s2_list[i]
      ctot = ctot_list[i]
      Sobs = Sobs_list[i]
      reflection_models.append(ReflectionProfileModel(profile_model, s0, s2, ctot, Sobs))

    L = 0
    dL = flex.double(6)
    I = flex.double(flex.grid(6,6))
    for r in reflection_models:
      L += r.log_likelihood()
      dL += r.first_derivatives()
      I += r.fisher_information()

    dL = matrix.col(list(dL))
    I = matrix.sqr(list(I))

    from scitbx import linalg
    LL = flex.double()
    for j in range(6):
      for i in range(j+1):
        LL.append(I[j*6+i])

    ll = linalg.l_l_transpose_cholesky_decomposition_in_place(LL)
    S = flex.double(-dL)
    S = ll.solve(S)

    def log_likelihood(params):
      parameterisation = MosaicityParameterisation(params)
      profile_model = ProfileModel(parameterisation)
      L = 0
      for i in range(len(s2_list)):
        s2 = s2_list[i]
        ctot = ctot_list[i]
        Sobs = Sobs_list[i]
        r = ReflectionProfileModel(profile_model, s0, s2, ctot, Sobs)
        L += r.log_likelihood()
      return L

    def negative_log_likelihood(params):
      return -log_likelihood(params)

    TINY = 1e-7
    delta = 1
    use_submitted = False
    while True:
      L_test = log_likelihood(params_old + delta*matrix.col(S))
      if L_test >= (L-TINY):
        break
      else:
        delta /= 2
      if delta < 1.0 / (2**10):
        #params = console(locals())
        params = gradient_descent_old(negative_log_likelihood, dL, params_old)
        use_submitted = True
        break 
        #raise RuntimeError("Cant find small step")
    
    if use_submitted == False:
      params = params_old + delta*matrix.col(S)
    M = matrix.sqr((
      params[0], 0, 0,
      params[1], params[2], 0,
      params[3], params[4], params[5]))
    sigma = M*M.transpose()
    
    if verbose:
      print "( %.3g, %.3g, %.3g, %.3g, %.3g, %.3g, %.3g, %.3g, %.3g, ): L = %f"  % (tuple(sigma) + (L,))

    if all(abs(a-b)<1e-7 for a, b in zip(params, params_old)):
      break
    params_old = params

  return params

if __name__ == '__main__':

  p = MosaicityParameterisation((1, 0.1, 2, 0.2, 0.3, 3))

  print p.sigma()
  print list(p.first_derivatives())
  print list(p.second_derivatives())
