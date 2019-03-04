from __future__ import division
from __future__ import print_function
from dials.array_family import flex
from scitbx import matrix
from scitbx import linalg
from scitbx.linalg import eigensystem
from dials_scratch.jmp.stills.potato.parameterisation import MosaicityParameterisation
from dials_scratch.jmp.stills.potato.profile_model import ReflectionProfileModelList
from dials.algorithms.profile_model.gaussian_rs import CoordinateSystem2d
from math import sqrt, pi


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
        print("NOO")
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

  def line_search(self, x, p, tau=0.5, delta=1.0, tolerance=1e-7):
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

    # Print initial
    self.callback(x0)

  def log_likelihood(self, x):
    '''
    :param x: The parameter estimate
    :return: The log likelihood at x

    '''
    return self.model(x).log_likelihood()

  def score(self, x):
    '''
    :param x: The parameter estimate
    :return: The score at x

    '''
    return self.model(x).first_derivatives()

  def score_and_fisher_information(self, x):
    '''
    :param x: The parameter estimate
    :return: The score and fisher information at x

    '''
    model = self.model(x)
    S = model.first_derivatives()
    I = model.fisher_information()
    return S, I

  def model(self, x):
    '''
    :param x: The parameter estimate
    :return: The model

    '''
    parameterisation = MosaicityParameterisation(x)
    model = ReflectionProfileModelList(
      parameterisation,
      self.s0,
      self.s2_list,
      self.ctot_list,
      self.Sobs_list)
    return model

  def callback(self, x):
    '''
    Handle and update in parameter values

    '''
    parameterisation = MosaicityParameterisation(x)
    sigma = parameterisation.sigma()
    lnL = self.log_likelihood(x)
    format_string = "( %.3g, %.3g, %.3g, %.3g, %.3g, %.3g, %.3g, %.3g, %.3g ): L = %f"
    print(format_string % (tuple(sigma) + (lnL,)))
    self.history.append(x)


class ProfileRefiner(object):
  '''
  High level profile refiner class that handles book keeping etc

  '''

  def __init__(self, data, initial_parameters=None):
    '''
    Set the data and initial parameters

    '''
    self.s0 = data.s0
    self.s2_list = data.s2_list
    self.ctot_list = data.ctot_list
    self.Sobs_list = data.Sobs_list
    self.initial_parameters = initial_parameters

  def refine(self):
    '''
    Perform the profile refinement

    '''

    # Set the initial parameter values
    if self.initial_parameters is None:
      x0 = matrix.col((1, 0, 1, 0, 0, 1))
    else:
      x0 = matrix.col((self.initial_parameters))

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

    # Get the covariance matrix
    sigma = ml.model(self.parameters).parameterisation().sigma()

    # Print the eigen values and vectors
    print_eigen_values_and_vectors(sigma)


class ProfileRefinerData(object):
  '''
  A class for holding the data needed for the profile refinement

  '''
  def __init__(self, s0, s2_list, ctot_list, xbar_list, Sobs_list):
    '''
    Init the data

    '''
    self.s0 = s0
    self.s2_list = s2_list
    self.ctot_list = ctot_list
    self.xbar_list = xbar_list
    self.Sobs_list = Sobs_list

  @classmethod
  def from_reflections(self, experiment, reflections):
    '''
    Generate the required data from the reflections

    '''

    # Get the beam vector
    s0 = matrix.col(experiment.beam.get_s0())

    # Get the reciprocal lattice vector
    s2_list = reflections['s2']

    # Initialise the list of observed intensities and covariances
    ctot_list = flex.double(len(s2_list))
    Sobs_list = flex.double(flex.grid(len(s2_list), 4))

    print("Computing observed covariance for %d reflections" % len(reflections))
    s0_length = s0.length()
    assert len(experiment.detector) == 1
    panel = experiment.detector[0]
    sbox = reflections['shoebox']
    for r in range(len(reflections)):

      # Create the coordinate system
      cs = CoordinateSystem2d(s0, s2_list[r])

      # Get data and compute total counts
      data = sbox[r].data
      mask = sbox[r].mask
      bgrd = sbox[r].background

      # Get array of vectors
      i0 = sbox[r].bbox[0]
      j0 = sbox[r].bbox[2]
      assert data.all()[0] == 1
      X = flex.vec2_double(flex.grid(data.all()[1], data.all()[2]))
      ctot = 0
      C = flex.double(X.accessor())
      for j in range(data.all()[1]):
        for i in range(data.all()[2]):
          c = data[0,j,i] - bgrd[0,j,i]
          if mask[0,j,i] == 5 and c > 0:
            ctot += c
            ii = i + i0
            jj = j + j0
            s = panel.get_pixel_lab_coord((ii+0.5,jj+0.5))
            s = matrix.col(s).normalize() * s0_length
            X[j,i] = cs.from_beam_vector(s)
            C[j,i] = c

      # Compute the mean vector
      xbar = matrix.col((0,0))
      for j in range(X.all()[0]):
        for i in range(X.all()[1]):
          x = matrix.col(X[j,i])
          xbar += C[j,i] * x
      xbar /= ctot

      # Compute the covariance matrix
      Sobs = matrix.sqr((0, 0, 0, 0))
      for j in range(X.all()[0]):
        for i in range(X.all()[1]):
          x = matrix.col(X[j,i])
          Sobs += (x-xbar)*(x-xbar).transpose()*C[j,i]

      # Add to the lists
      ctot_list[r] = ctot
      Sobs_list[r,0] = Sobs[0]
      Sobs_list[r,1] = Sobs[1]
      Sobs_list[r,2] = Sobs[2]
      Sobs_list[r,3] = Sobs[3]

    # Print some information
    print("I_min = %.2f, I_max = %.2f" % (flex.min(ctot_list),
                                          flex.max(ctot_list)))

    # Print the mean covariance
    Smean = matrix.sqr((0,0,0,0))
    for r in range(Sobs_list.all()[0]):
      Smean += matrix.sqr(tuple(Sobs_list[r:r+1,:]))
    Smean /= Sobs_list.all()[0]
    print("")
    print("Mean observed covariance:")
    print_matrix(Smean)

    # Return the profile refiner data
    return ProfileRefinerData(s0, s2_list, ctot_list, Sobs_list)


def print_eigen_values_and_vectors(A):
  '''
  Print the eigen values and vectors of a matrix

  '''

  # Compute the eigen decomposition of the covariance matrix
  eigen_decomposition = eigensystem.real_symmetric(A.as_flex_double_matrix())
  Q = matrix.sqr(eigen_decomposition.vectors())
  L = matrix.diag(eigen_decomposition.values())

  # Print the matrix eigen values
  print("")
  print("Eigen Values:")
  print("")
  print_matrix(L, indent=2)
  print("")

  print("Eigen Vectors:")
  print("")
  print_matrix(Q, indent=2)
  print("")

  print("Mosaicity in degrees equivalent units")
  print("M1: %.5f degrees" % (sqrt(L[0])*180.0/pi))
  print("M2: %.5f degrees" % (sqrt(L[4])*180.0/pi))
  print("M3: %.5f degrees" % (sqrt(L[8])*180.0/pi))

def print_matrix(A, fmt='%.3g', indent=0):
  '''
  Pretty print matrix

  '''
  t = [fmt % a for a in A]
  l = [len(tt) for tt in t]
  max_l = max(l)
  fmt = "%" + ("%d" % (max_l+1)) + "s"
  prefix = ' ' * indent
  lines = []
  for j in range(A.n[0]):
    line = ''
    for i in range(A.n[1]):
      line += fmt % t[i+j*A.n[1]]
    lines.append("%s|%s|" % (prefix, line))
  print('\n'.join(lines))
