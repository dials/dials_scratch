from numpy.random import normal, uniform, multivariate_normal
from math import exp, pi, sqrt, log, sin, cos
from dials.array_family import flex
from scitbx import matrix
from dials_scratch.jmp.stills.potato.profile_model import compute_change_of_basis_operation

def generate_simple(s0, sigma, N = 100):
  '''
  Generate a list of normally distributed observations

  '''
  s2_list = []
  ctot_list = []
  xbar_list = []
  Sobs_list = []

  # Loop through the list
  for k in range(N):

    # Compute position in reciprocal space of the centre of the rlp
    s2_direction = matrix.col((
      uniform(0, 1),
      uniform(0, 1),
      uniform(0, 1))).normalize()

    # Rotate the covariance matrix
    R = compute_change_of_basis_operation(s0, s2_direction)
    sigmap = R*sigma*R.transpose()
    s2_magnitude = normal(s0.length(), sqrt(sigmap[8]))
    s2 = s2_direction * s2_magnitude

    # Rotate to get mu
    mu = R*s2
    assert abs(mu.normalize().dot(matrix.col((0, 0, 1))) - 1) < 1e-7

    # Partition the matrix
    sigma11 = matrix.sqr((
      sigmap[0], sigmap[1],
      sigmap[3], sigmap[4]))
    sigma12 = matrix.col((sigmap[2], sigmap[5]))
    sigma21 = matrix.col((sigmap[6], sigmap[7])).transpose()
    sigma22 = sigmap[8]

    # Compute the conditional distribution
    mu1 = matrix.col((mu[0], mu[1]))
    mu2 = mu[2]
    z = s0.length()
    sigma_bar = sigma11 - sigma12*(1/sigma22)*sigma21
    mu_bar = mu1 + sigma12*(z-mu2)/sigma22

    # Compute the scale factor and intensity
    scale = exp(-0.5 * (z-mu2)**2 / sigma22)
    I = uniform(50, 1000)
    if I <= 1:
      continue

    # Simulate some observations
    points = multivariate_normal(mu_bar, sigma_bar.as_list_of_lists(), int(I))
    
    # Compute the observed mean for each observation
    ctot = 0
    xbar = matrix.col((0, 0))
    for x in points:
      xbar += matrix.col(x)
  
    ctot = len(points)

    xbar /= ctot
    
    # Compute the observed covariance for each observation
    Sobs = matrix.sqr((0, 0, 0, 0))
    for x in points:
      x = matrix.col(x)
      Sobs += (x-xbar)*(x-xbar).transpose()

    s2_list.append(s2)
    ctot_list.append(ctot)
    xbar_list.append(xbar)
    Sobs_list.append(list(Sobs))

  return s2_list, ctot_list, xbar_list, Sobs_list


def generate_simple_binned(s0, sigma, N = 100):
  '''
  Generate a list of normally distributed observations

  '''
  s2_list = []
  ctot_list = []
  xbar_list = []
  Sobs_list = []

  # Loop through the list
  for k in range(N):

    # Compute position in reciprocal space of the centre of the rlp
    s2_direction = matrix.col((
      uniform(0, 1),
      uniform(0, 1),
      uniform(0, 1))).normalize()

    # Rotate the covariance matrix
    R = compute_change_of_basis_operation(s0, s2_direction)
    sigmap = R*sigma*R.transpose()
    s2_magnitude = normal(s0.length(), sqrt(sigmap[8]))
    s2 = s2_direction * s2_magnitude

    # Rotate to get mu
    mu = R*s2
    assert abs(mu.normalize().dot(matrix.col((0, 0, 1))) - 1) < 1e-7

    # Partition the matrix
    sigma11 = matrix.sqr((
      sigmap[0], sigmap[1],
      sigmap[3], sigmap[4]))
    sigma12 = matrix.col((sigmap[2], sigmap[5]))
    sigma21 = matrix.col((sigmap[6], sigmap[7])).transpose()
    sigma22 = sigmap[8]

    # Compute the conditional distribution
    mu1 = matrix.col((mu[0], mu[1]))
    mu2 = mu[2]
    z = s0.length()
    sigma_bar = sigma11 - sigma12*(1/sigma22)*sigma21
    mu_bar = mu1 + sigma12*(z-mu2)/sigma22

    # Compute the scale factor and intensity
    scale = exp(-0.5 * (z-mu2)**2 / sigma22)
    I = uniform(50, 1000)
    if I <= 1:
      continue

    # Simulate some observations
    a = 7.5
    b = 15.0 / (12*sqrt(sigma_bar[0]))
    c = 15.0 / (12*sqrt(sigma_bar[3]))
    D = flex.double(flex.grid(15,15))
    points = multivariate_normal(mu_bar, sigma_bar.as_list_of_lists(), int(I))
    for x, y in points:
      i = int(a + b*x)
      j = int(a + c*y)
      if j >= 0 and j < D.all()[0] and i >= 0 and i < D.all()[1]:
        D[j,i] += 1
    
    # Compute the observed mean for each observation
    ctot = 0
    xbar = matrix.col((0, 0))
    for j in range(D.all()[0]):
      for i in range(D.all()[1]):
        ctot += D[j,i]
        x = matrix.col((
          (i+0.5-a)/b,
          (j+0.5-a)/c))
        xbar += x*D[j,i]

    if ctot <= 0:
      continue

    xbar /= ctot

    # Compute the observed covariance for each observation
    Sobs = matrix.sqr((0, 0, 0, 0))
    for j in range(D.all()[0]):
      for i in range(D.all()[1]):
        x = matrix.col((
          (i+0.5-a)/b,
          (j+0.5-a)/c))
        Sobs += (x-xbar)*(x-xbar).transpose()*D[j,i]

    s2_list.append(s2)
    ctot_list.append(ctot)
    xbar_list.append(xbar)
    Sobs_list.append(list(Sobs))

  return s2_list, ctot_list, xbar_list, Sobs_list


def generate_from_reflections(s0, sigma, reflections):
  '''
  Generate a list of normally distributed observations

  '''
  s2_list = []
  ctot_list = []
  xbar_list = []
  Sobs_list = []

  # Loop through the list
  for k in range(len(reflections)):

    # Compute position in reciprocal space of the centre of the rlp
    s2 = matrix.col(reflections['s2'][k])

    # Rotate the covariance matrix
    R = compute_change_of_basis_operation(s0, s2)
    sigmap = R*sigma*R.transpose()

    # Rotate to get mu
    mu = R*s2
    assert abs(mu.normalize().dot(matrix.col((0, 0, 1))) - 1) < 1e-7

    # Partition the matrix
    sigma11 = matrix.sqr((
      sigmap[0], sigmap[1],
      sigmap[3], sigmap[4]))
    sigma12 = matrix.col((sigmap[2], sigmap[5]))
    sigma21 = matrix.col((sigmap[6], sigmap[7])).transpose()
    sigma22 = sigmap[8]

    # Compute the conditional distribution
    mu1 = matrix.col((mu[0], mu[1]))
    mu2 = mu[2]
    z = s0.length()
    sigma_bar = sigma11 - sigma12*(1/sigma22)*sigma21
    mu_bar = mu1 + sigma12*(z-mu2)/sigma22

    # Perform rejection sampling to get a normally distributed set of
    # reflections
    P = exp(-0.5*(z-mu2)**2 / sigma22)
    R = uniform(0, 1.0)
    if P < R:
      continue

    # Compute the scale factor and intensity
    scale = exp(-0.5 * (z-mu2)**2 / sigma22)
    I = uniform(50, 1000)
    if I <= 1:
      continue

    # Simulate some observations
    points = multivariate_normal(mu_bar, sigma_bar.as_list_of_lists(), int(I))
    
    # Compute the observed mean for each observation
    ctot = 0
    xbar = matrix.col((0, 0))
    for x in points:
      xbar += matrix.col(x)
  
    ctot = len(points)

    xbar /= ctot
    
    # Compute the observed covariance for each observation
    Sobs = matrix.sqr((0, 0, 0, 0))
    for x in points:
      x = matrix.col(x)
      Sobs += (x-xbar)*(x-xbar).transpose()

    s2_list.append(s2)
    ctot_list.append(ctot)
    xbar_list.append(xbar)
    Sobs_list.append(list(Sobs))

  return s2_list, ctot_list, xbar_list, Sobs_list


def generate_from_reflections_binned(s0, sigma, reflections):
  '''
  Generate a list of normally distributed observations

  '''
  s2_list = []
  ctot_list = []
  xbar_list = []
  Sobs_list = []

  # Loop through the list
  for k in range(len(reflections)):

    # Compute position in reciprocal space of the centre of the rlp
    s2 = matrix.col(reflections['s2'][k])

    # Rotate the covariance matrix
    R = compute_change_of_basis_operation(s0, s2)
    sigmap = R*sigma*R.transpose()

    # Rotate to get mu
    mu = R*s2
    assert abs(mu.normalize().dot(matrix.col((0, 0, 1))) - 1) < 1e-7

    # Partition the matrix
    sigma11 = matrix.sqr((
      sigmap[0], sigmap[1],
      sigmap[3], sigmap[4]))
    sigma12 = matrix.col((sigmap[2], sigmap[5]))
    sigma21 = matrix.col((sigmap[6], sigmap[7])).transpose()
    sigma22 = sigmap[8]

    # Compute the conditional distribution
    mu1 = matrix.col((mu[0], mu[1]))
    mu2 = mu[2]
    z = s0.length()
    sigma_bar = sigma11 - sigma12*(1/sigma22)*sigma21
    mu_bar = mu1 + sigma12*(z-mu2)/sigma22

    # Perform rejection sampling to get a normally distributed set of
    # reflections
    P = exp(-0.5*(z-mu2)**2 / sigma22)
    R = uniform(0, 1.0)
    if P < R:
      continue

    # Compute the scale factor and intensity
    scale = exp(-0.5 * (z-mu2)**2 / sigma22)
    I = uniform(50, 1000)
    if I <= 1:
      continue

    # Simulate some observations
    a = 7.5
    b = 15.0 / (12*sqrt(sigma_bar[0]))
    c = 15.0 / (12*sqrt(sigma_bar[3]))
    D = flex.double(flex.grid(15,15))
    points = multivariate_normal(mu_bar, sigma_bar.as_list_of_lists(), int(I))
    for x, y in points:
      i = int(a + b*x)
      j = int(a + c*y)
      if j >= 0 and j < D.all()[0] and i >= 0 and i < D.all()[1]:
        D[j,i] += 1
    
    # Compute the observed mean for each observation
    ctot = 0
    xbar = matrix.col((0, 0))
    for j in range(D.all()[0]):
      for i in range(D.all()[1]):
        ctot += D[j,i]
        x = matrix.col((
          (i+0.5-a)/b,
          (j+0.5-a)/c))
        xbar += x*D[j,i]

    if ctot <= 0:
      continue

    xbar /= ctot

    # Compute the observed covariance for each observation
    Sobs = matrix.sqr((0, 0, 0, 0))
    for j in range(D.all()[0]):
      for i in range(D.all()[1]):
        x = matrix.col((
          (i+0.5-a)/b,
          (j+0.5-a)/c))
        Sobs += (x-xbar)*(x-xbar).transpose()*D[j,i]

    s2_list.append(s2)
    ctot_list.append(ctot)
    xbar_list.append(xbar)
    Sobs_list.append(list(Sobs))

  return s2_list, ctot_list, xbar_list, Sobs_list
