from __future__ import division
from dials.array_family import flex
from scitbx import matrix
from scitbx import linalg
from dials_scratch.jmp.stills.potato.parameterisation import MosaicityParameterisation
from dials_scratch.jmp.stills.potato.profile_model import ReflectionProfileModelList


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
    if fb <= fa:
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

  def solve_update_equation(self, S, I):
    '''
    Solve the update equation using cholesky decomposition
    :param S: The score
    :param I: The fisher information
    :return: The parameter delta

    '''

    # Construct triangular matrix
    LL = flex.double()
    for j in range(len(S)):
      for i in range(j+1):
        LL.append(I[j*len(S)+i])

    # Perform the decomposition
    ll = linalg.l_l_transpose_cholesky_decomposition_in_place(LL)
    p = flex.double(S)
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

    return line_search(f, x, p, tolerance=self.tolerance)

  def gradient_search(self, x0):
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

    # Store the parameter history
    self.history = []

  def log_likelihood(self, x):
    '''
    :param x: The parameter estimate
    :return: The log likelihood at x

    '''
    parameterisation = MosaicityParameterisation(x)
    model = ReflectionProfileModelList(
      parameterisation,
      self.s0,
      self.s2_list,
      self.ctot_list,
      self.Sobs_list)
    return model.log_likelihood()

  def score(self, x):
    '''
    :param x: The parameter estimate
    :return: The score at x

    '''
    parameterisation = MosaicityParameterisation(x)
    model = ReflectionProfileModelList(
      parameterisation,
      self.s0,
      self.s2_list,
      self.ctot_list,
      self.Sobs_list)
    return model.first_derivatives()

  def score_and_fisher_information(self, x):
    '''
    :param x: The parameter estimate
    :return: The score and fisher information at x

    '''
    parameterisation = MosaicityParameterisation(x)
    model = ReflectionProfileModelList(
      parameterisation,
      self.s0,
      self.s2_list,
      self.ctot_list,
      self.Sobs_list)
    S = model.first_derivatives()
    I = model.fisher_information()
    return S, I

  def callback(self, x):
    '''
    Handle and update in parameter values

    '''
    parameterisation = MosaicityParameterisation(x)
    sigma = parameterisation.sigma()
    lnL = self.log_likelihood(x)
    format_string = "( %.3g, %.3g, %.3g, %.3g, %.3g, %.3g, %.3g, %.3g, %.3g ): L = %f"
    print format_string % (tuple(sigma) + (lnL,))
    self.history.append(x)


class ProfileRefiner(object):
  '''
  Top level profile refiner class that handles book keeping etc

  '''

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
