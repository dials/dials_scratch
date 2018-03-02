
from scitbx import matrix
from math import pi, cos, sin, exp, sqrt, atan2, log, tan

def normal_2d(x, mu, sigma):
  A = 1.0 / sqrt((2*pi)**2 * sigma.determinant())
  B = ((x - mu).transpose() * sigma.inverse() * (x - mu))[0]
  return A * exp(-0.5 * B)

def compute_peak_f(mu, sigma):

  p = 0
  w = 0

  N = 10000
  for i in range(N):
    theta = 0 + i*2*pi / N
    x = cos(theta)
    y = sin(theta)
    v = matrix.col((x, y))
    f = normal_2d(v, mu, sigma)

    if f > p:
      p = f
      w = theta

  return w

def compute_peak_g(mu, sigma):

  p = 0
  w = 0

  N = 10000
  for i in range(N):
    theta = 0 + i*2*pi / N
    x = cos(theta)
    y = sin(theta)
    v = matrix.col((x, y))
    f = theta*normal_2d(v, mu, sigma)

    if f > p:
      p = f
      w = theta

  return w

def compute_peak_g1(mu, sigma):

  p = 0
  w = 0

  N = 10000
  for i in range(N):
    theta = 0 + i*2*pi / N
    x = cos(theta)
    y = sin(theta)
    v = matrix.col((x, y))
    f = x*normal_2d(v, mu, sigma)

    if f > p:
      p = f
      w = x

  return w


def compute_peak_g1(mu, sigma):

  p = None
  w = 0
  N = 10000
  for i in range(N):
    theta = 0.001 + i*(pi/2-0.001) / N
    x = cos(theta)
    y = sin(theta)
    v = matrix.col((x, y))
    f = 0.5 * ((v - mu).transpose() * sigma.inverse() * (v - mu))[0] - log(x)
    #f = x*normal_2d(v, mu, sigma)

    if p is None or f < p:
      p = f
      w = x

  return w

def compute_peak_g2(mu, sigma):

  p = 0
  w = 0

  N = 10000
  for i in range(N):
    theta = 0 + i*2*pi / N
    x = cos(theta)
    y = sin(theta)
    v = matrix.col((x, y))
    f = y*normal_2d(v, mu, sigma)

    if f > p:
      p = f
      w = y

  return w

def compute_mean(mu, sigma):

  I1 = matrix.col((0,0))
  I2 = 0

  N = 10000
  for i in range(N):
    theta = 0 + i*2*pi / N
    x = cos(theta)
    y = sin(theta)
    v = matrix.col((x, y))
    f = normal_2d(v, mu, sigma)
    I1 += f * v
    I2 += f

  return (I1 / I2).normalize()


def func(theta, mu, sigma):
  x = matrix.col((cos(theta), sin(theta)))
  x1, x2 = x
  mu1, mu2 = mu
  s1, _, _, s2 = sigma
  return 0.5 * (x1 - mu1)**2 / s1 + 0.5 * (x2 - mu2)**2 / s2 - log(x[0])
  #return 0.5 * ((x - mu).transpose()*sigma.inverse()*(x-mu))[0] - log(x[0])

def compute_peak_g1_it(mu, sigma):

  #theta0 = atan2(sqrt(1 - 0.700506727462**2), 0.700506727462)
  theta0 = atan2(mu[1], mu[0])
  m1, m2 = mu
  s1, _, _, s2 = sigma

  while True:
    A = sin(theta0)*cos(theta0)*(1/s2 - 1/s1)
    B = sin(theta0)*(m1 / s1)
    C = cos(theta0)*(m2 / s2)
    D = tan(theta0)

    E = cos(2*theta0)*(1/s2 - 1/s1)
    F = cos(theta0)*(m1 / s1)
    G = sin(theta0)*(m2 / s2)
    H = 1/cos(theta0)**2

    df = A + B - C + D
    d2f = E + F + G + H

    theta = theta0 - df / d2f
    # print theta, theta0, df, func(theta, mu, sigma), func(theta0, mu, sigma)
    if abs(theta - theta0) < 1e-7:
      break
    theta0 = theta

  # print "---"
  x = matrix.col((cos(theta), sin(theta)))
  # x0 = matrix.col((mu[0], mu[1], 0))

  # while True:
  #   m1, m2 = mu
  #   s1, _, _, s2 = sigma
  #   x1, x2, l = x0

  #   A = (x1 - m1) / s1**2 - 1 / x1 + 2*l*x1
  #   B = (x2 - m2) / s2**2 + 2*l*x2
  #   C = (1 / s1**2 + 1 / x1**2 + 2*l)
  #   D = (1 / s2**2 + 2*l)
  #   E = (x1**2 + x2**2 - 1)
  #   dLdx1 = 2*A*C + 2*E*(2*x1)
  #   dLdx2 = 2*B*D + 2*E*(2*x2)
  #   dLdx3 = 2*A*(2*x1) + 2*B*(2*x2)

  #   G = matrix.col((dLdx1, dLdx2, dLdx3))

  #   x = x0 -0.00000001 * G
  #   print tuple(G), tuple(x), matrix.col((x[0], x[1])).length()


  #   if (x - x0).length() < 1e-7:
  #     break
  #   x0 = x

  # while True:

  #   m1, m2 = mu
  #   s1, _, _, s2 = sigma
  #   x1, x2, l = x0

  #   H = matrix.sqr((
  #     1/s1**2 + 1/x1**2 + 2*l, 0, 2*x1,
  #     0, 1/s2**2 + 2*l, 2*x2,
  #     2*x1, 2*x2, 0))

  #   G = matrix.col((
  #     (x1-m1)/s1**2 - 1/x1 + 2*l*x1,
  #     (x2-m2)/s2**2 + 2*l*x2,
  #     x1**2 + x2**2 -1))

  #   x = x0 - H.inverse() * G
  #   print tuple(x)

  #   if (x - x0).length() < 1e-7:
  #     break
  #   x0 = x

  #x = matrix.col((x[0], x[1]))
  return x

def compute_peak_g2_it(mu, sigma):

  theta0 = atan2(mu[1], mu[0])
  m1, m2 = mu
  s1, _, _, s2 = sigma

  while True:
    A = sin(theta0)*cos(theta0)*(1/s2 - 1/s1)
    B = sin(theta0)*(m1 / s1)
    C = cos(theta0)*(m2 / s2)
    D = 1/tan(theta0)

    E = cos(2*theta0)*(1/s2 - 1/s1)
    F = cos(theta0)*(m1 / s1)
    G = sin(theta0)*(m2 / s2)
    H = 1/sin(theta0)**2

    df = A + B - C - D
    d2f = E + F + G + H

    theta = theta0 - df / d2f
    print theta
    # print theta, theta0, df, func(theta, mu, sigma), func(theta0, mu, sigma)
    if abs(theta - theta0) < 1e-7:
      break
    theta0 = theta

  # print "---"
  x = matrix.col((cos(theta), sin(theta)))

  # x0 = matrix.col((mu[0], mu[1], 0))

  # while True:

  #   m1, m2 = mu
  #   s1, _, _, s2 = sigma
  #   x1, x2, l = x0

  #   H = matrix.sqr((
  #     1/s1**2 + 2*l, 0, 2*x1,
  #     0, 1/s2**2 +1/x2**2 + 2*l, 2*x2,
  #     2*x1, 2*x2, 0))

  #   G = matrix.col((
  #     (x1-m1)/s1**2  + 2*l*x1,
  #     (x2-m2)/s2**2 -1/x1+ 2*l*x2,
  #     x1**2 + x2**2 -1))

  #   x = x0 - H.inverse() * G
  #   print tuple(x)

  #   if (x - x0).length() < 1e-7:
  #     break
  #   x0 = x

  # x = matrix.col((x[0], x[1]))
  # print x.length()
  return x

# def f(x, mu1, mu2, s1, s2):
#   A = mu1**2 / (1 - x*s1)**2
#   B = mu2**2 / (1 - x*s2)**2
#   return A + B - 1

# def df(x, mu1, mu2, s1, s2):
#   A = 2*s1 * mu1**2 / (1 - x*s1)**3
#   B = 2*s2 * mu2**2 / (1 - x*s2)**3
#   return A + B

# def d2f(x, mu1, mu2, s1, s2):
#   A = 6*s1 * mu1**2 / (1 - x*s1)**4
#   B = 6*s2 * mu2**2 / (1 - x*s2)**4
#   return A + B

def f(x, mu1, mu2, s1, s2):
  A = mu1**2 * (1 - x*s2)**2
  B = mu2**2 * (1 - x*s1)**2
  C = (1 - x*s1)**2 * (1 -  x*s2)**2
  return A * B - C

def df(x, mu1, mu2, s1, s2):
  A = mu1**2 * s2 * (1 - x * s2)
  B = mu2**2 * s1 * (1 - x * s1)
  C = s1 * (1 - x * s2)**2 * (1 - x * s1)
  D = s2 * (1 - x * s1)**2 * (1 - x * s2)
  return -2 * (A + B - C - D)

def d2f(x, mu1, mu2, s1, s2):
  A = mu1**2 * s2**2
  B = mu2**2 * s1**2
  C = s1**2 * (1 - x * s2)**2
  D = s2**2 * (1 - x * s1)**2
  E = 4 * s1 * s2 * (1 - x*s1)*(1 - x*s2)
  return 2 * (A + B - C - D - E)

def compute_l0(x, mu, s):
  return (1-mu[0]/x[0])/s[0]**2

def compute_peak_f_it(mu, sigma):

  theta0 = atan2(mu[1], mu[0])
  m1, m2 = mu
  s1, _, _, s2 = sigma

  while True:
    A = sin(theta0)*cos(theta0)*(1/s2 - 1/s1)
    B = sin(theta0)*(m1 / s1)
    C = cos(theta0)*(m2 / s2)

    E = cos(2*theta0)*(1/s2 - 1/s1)
    F = cos(theta0)*(m1 / s1)
    G = sin(theta0)*(m2 / s2)

    df = A + B - C
    d2f = E + F + G

    theta = theta0 - df / d2f
    print "P: ", theta
    # print theta, theta0, df, func(theta, mu, sigma), func(theta0, mu, sigma)
    if abs(theta - theta0) < 1e-7:
      break
    theta0 = theta

  x = matrix.col((cos(theta), sin(theta)))

  # # Compute the initial value of the lagrange multiplier
  # l0 = compute_l0(mu.normalize(), mu, sigma)

  # # Compute the value via Halley's method
  # while True:
  # #for i in range(20):
  #   mu1, mu2 = mu
  #   s1, s2 = sigma[0], sigma[3]
  #   U = f(l0, mu1, mu2, s1, s2)
  #   Up = df(l0, mu1, mu2, s1, s2)
  #   Up2 = d2f(l0, mu1, mu2, s1, s2)
  #   l = l0 - 2*U*Up / (2*Up**2 - U*Up2)
  #   print l, df(l, mu1, mu2, s1, s2)

  #   if abs(l - l0) < 1e-7:
  #     break
  #   l0 = l

  # # Compute the mode
  # x = matrix.col((
  #   mu[0] / (1 - l * sigma[0]),
  #   mu[1] / (1 - l * sigma[3])))

  # print x
  # exit(0)

  # Now transform back into the original basis
  return x


def compute_mean_estimate(peak_f, peak_fx, peak_fy, mu, sigma):
  sigma_inv = sigma.inverse()

  H_f = -matrix.sqr((
    sigma_inv[0], 0.5*(sigma_inv[1] + sigma_inv[2]),
    0.5*(sigma_inv[1] + sigma_inv[2]), sigma_inv[3]))

  H_g1 = -matrix.sqr((
    sigma_inv[0] + 1/peak_fx[0]**2, 0.5*(sigma_inv[1] + sigma_inv[2]),
    0.5*(sigma_inv[1] + sigma_inv[2]), sigma_inv[3]))

  H_g2 = -matrix.sqr((
    sigma_inv[0], 0.5*(sigma_inv[1] + sigma_inv[2]),
    0.5*(sigma_inv[1] + sigma_inv[2]), sigma_inv[3] + 1/peak_fy[1]**2))

  A = (2*pi) * (1/sqrt((-H_f).determinant())) * exp(
    -0.5*((peak_f-mu).transpose() * sigma_inv * (peak_f-mu))[0])

  B = (2*pi) * (1/sqrt((-H_g1).determinant())) * peak_fx[0] * exp(
    -0.5*((peak_fx-mu).transpose() * sigma_inv * (peak_fx-mu))[0])

  C = (2*pi) * (1/sqrt((-H_g2).determinant())) * peak_fy[1] * exp(
    -0.5*((peak_fy-mu).transpose() * sigma_inv * (peak_fy-mu))[0])

  print A, B, C
  print tuple(peak_f)
  xc = B / A
  yc = C / A

  x = matrix.col((xc, yc)).normalize()

  return x



mu = matrix.col((1, 1)).normalize() * 0.95
sigma = matrix.sqr((0.03, 0,
                    0, 0.01))
sigma_inv = sigma.inverse()

mean = compute_mean(mu, sigma)
peak_f = compute_peak_f(mu, sigma)
peak_g = compute_peak_g(mu, sigma)
peak_x = compute_peak_g1(mu, sigma)
peak_y = compute_peak_g2(mu, sigma)

peak_f_it = compute_peak_f_it(mu, sigma)
peak_g1_it = compute_peak_g1_it(mu, sigma)
peak_g2_it = compute_peak_g2_it(mu, sigma)


x_f_search = matrix.col((cos(peak_f), sin(peak_f)))

x_g_search = matrix.col((cos(peak_g), sin(peak_g)))
x_g1_search = matrix.col((peak_x, sqrt(1 - peak_x**2)))
x_g2_search = matrix.col((sqrt(1 - peak_y**2), peak_y))


print "---"
print peak_x, peak_g1_it
print peak_y, peak_g2_it
print "---"

x_f = peak_f_it
x_g1 = peak_g1_it
x_g2 = peak_g2_it

x_search1 = compute_mean_estimate(x_f_search, x_g1_search, x_g2_search, mu, sigma)
x_search2 = compute_mean_estimate(x_f_search, x_g1_search, x_g2_search, mu, sigma)

x = compute_mean_estimate(x_f, x_g1, x_g2, mu, sigma)


print "Mean: ", tuple(mean)
print "Peak f(search): ", tuple(x_f_search), mean.angle(x_f_search)
print "Peak g1(search): ", tuple(x_g1_search), mean.angle(x_g1_search)
print "Peak g2(search): ", tuple(x_g2_search), mean.angle(x_g2_search)
print "Peak f: ", tuple(x_f), mean.angle(x_f)
print "Peak g1: ", tuple(x_g1), mean.angle(x_g1)
print "Peak g2: ", tuple(x_g2), mean.angle(x_g2)
print "Estimate (search 1): ", tuple(x_search1), mean.angle(x_search1)
print "Estimate (search 2): ", tuple(x_search2), mean.angle(x_search2)
print "Esimate: ", tuple(x), mean.angle(x)
print "Mu: ", tuple(mu), mean.angle(mu)
