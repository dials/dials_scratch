
from scitbx import matrix
from math import pi, cos, sin, exp, sqrt, atan2, log, tan, acos

def normal_3d(x, mu, sigma):
  A = 1.0 / sqrt((2*pi)**3 * sigma.determinant())
  B = ((x - mu).transpose() * sigma.inverse() * (x - mu))[0]
  return A * exp(-0.5 * B)

def func_f(theta, phi, mu, s_sq):
  sin_theta = sin(theta)
  cos_theta = cos(theta)
  sin_phi = sin(phi)
  cos_phi = cos(phi)
  A = (sin_theta * cos_phi - mu[0])**2 / s_sq[0]
  B = (sin_theta * sin_phi - mu[1])**2 / s_sq[1]
  C = (cos_theta           - mu[2])**2 / s_sq[2]
  return 0.5 * (A + B + C)

def func_f_gradient(theta, phi, mu, s_sq):
    sin_theta = sin(theta)
    cos_theta = cos(theta)
    sin_phi = sin(phi)
    cos_phi = cos(phi)
    A = (sin_theta * cos_theta * cos_phi**2 - mu[0]*cos_theta*cos_phi) / s_sq[0]
    B = (sin_theta * cos_theta * sin_phi**2 - mu[1]*cos_theta*sin_phi) / s_sq[1]
    C = (cos_theta * sin_theta - mu[2]*sin_theta) / s_sq[2]
    D = (sin_theta**2 * sin_phi * cos_phi - mu[0]*sin_theta*sin_phi) / s_sq[0]
    E = (sin_theta**2 * sin_phi * cos_phi - mu[1]*sin_theta*cos_phi) / s_sq[1]
    return matrix.col((
      A + B - C,
      -D + E))

def func_f_hessian(theta, phi, mu, s_sq):
  sin_theta = sin(theta)
  cos_theta = cos(theta)
  sin_phi = sin(phi)
  cos_phi = cos(phi)
  cos_2theta = cos(2*theta)
  sin_2theta = sin(2*theta)
  sin_2phi = sin(2*phi)
  cos_2phi = cos(2*phi)
  A = (cos_2theta * cos_theta**2 + mu[0]*sin_theta*cos_phi) / s_sq[0]
  B = (cos_2theta * sin_theta**2 + mu[1]*sin_theta*sin_phi) / s_sq[1]
  C = (cos_2theta - mu[2]*cos_theta) / s_sq[2]
  D = (0.5 * sin_2theta*sin_2phi - mu[0]*cos_theta*sin_phi) / s_sq[0]
  E = (0.5 * sin_2theta*sin_2phi - mu[1]*cos_theta*cos_phi) / s_sq[1]
  F = (sin_theta**2 * cos_2phi - mu[0]*sin_theta*cos_phi) / s_sq[0]
  G = (sin_theta**2 * cos_2phi + mu[1]*sin_theta*sin_phi) / s_sq[1]
  H1 = A + B - C
  H2 = -D + E
  H3 = -F + G
  return matrix.sqr((
    H1, H2,
    H2, H3))

def compute_peak_f(mu, sigma, max_iter=100):

  def gradient(theta, phi, mu, s_sq):
    return func_f_gradient(theta, phi, mu, s_sq)

  def hessian(theta, phi, mu, s_sq):
    return func_f_hessian(theta, phi, mu, s_sq)

  phi = atan2(mu[1], mu[0])
  theta = acos(mu[2] / mu.length())

  x0 = matrix.col((theta, phi))

  s_sq = matrix.col((sigma[0], sigma[4], sigma[8]))

  for it in range(max_iter):

    G = gradient(x0[0], x0[1], mu, s_sq)
    H = hessian(x0[0], x0[1], mu, s_sq)
    x = x0 - H.inverse() * G
    if (x - x0).length() < 1e-7:
      break

    x0 = x

  assert(it < max_iter-1)

  return x

def compute_peak_fx(mu, sigma, max_iter=100):

  def gradient(theta, phi, mu, s_sq):
    G = func_f_gradient(theta, phi, mu, s_sq)
    return G + matrix.col((
      - 1 / tan(theta),
      tan(phi)))

  def hessian(theta, phi, mu, s_sq):
    H = func_f_hessian(theta, phi, mu, s_sq)
    return H + matrix.sqr((
      1 / sin(theta)**2, 0,
      0, 1 / cos(phi)**2))

  phi = atan2(mu[1], mu[0])
  theta = acos(mu[2] / mu.length())

  x0 = matrix.col((theta, phi))

  s_sq = matrix.col((sigma[0], sigma[4], sigma[8]))

  for it in range(max_iter):

    G = gradient(x0[0], x0[1], mu, s_sq)
    H = hessian(x0[0], x0[1], mu, s_sq)
    x = x0 - H.inverse() * G
    if (x - x0).length() < 1e-7:
      break

    x0 = x

  assert(it < max_iter-1)

  return x

def compute_peak_fy(mu, sigma, max_iter=100):

  def gradient(theta, phi, mu, s_sq):
    G = func_f_gradient(theta, phi, mu, s_sq)
    return G + matrix.col((
      - 1 / tan(theta),
      - 1 / tan(phi)))

  def hessian(theta, phi, mu, s_sq):
    H = func_f_hessian(theta, phi, mu, s_sq)
    return H + matrix.sqr((
      1 / sin(theta)**2, 0,
      0, 1 / sin(phi)**2))

  phi = atan2(mu[1], mu[0])
  theta = acos(mu[2] / mu.length())

  x0 = matrix.col((theta, phi))

  s_sq = matrix.col((sigma[0], sigma[4], sigma[8]))

  for it in range(max_iter):

    G = gradient(x0[0], x0[1], mu, s_sq)
    H = hessian(x0[0], x0[1], mu, s_sq)
    x = x0 - H.inverse() * G
    if (x - x0).length() < 1e-7:
      break

    x0 = x

  assert(it < max_iter-1)

  return x

def compute_peak_fz(mu, sigma, max_iter=100):

  def gradient(theta, phi, mu, s_sq):
    G = func_f_gradient(theta, phi, mu, s_sq)
    return G + matrix.col((
      tan(theta),
      0))

  def hessian(theta, phi, mu, s_sq):
    H = func_f_hessian(theta, phi, mu, s_sq)
    return H + matrix.sqr((
      1 / cos(theta)**2, 0,
      0, 0))

  phi = atan2(mu[1], mu[0])
  theta = acos(mu[2] / mu.length())

  x0 = matrix.col((theta, phi))

  s_sq = matrix.col((sigma[0], sigma[4], sigma[8]))

  for it in range(max_iter):

    G = gradient(x0[0], x0[1], mu, s_sq)
    H = hessian(x0[0], x0[1], mu, s_sq)
    x = x0 - H.inverse() * G
    if (x - x0).length() < 1e-7:
      break

    x0 = x

  assert(it < max_iter-1)



  return x


def compute_mean_estimate(peak_f, peak_fx, peak_fy, peak_fz, mu, sigma):

  sigma_inv = sigma.inverse()

  s1, s2, s3 = sigma[0], sigma[4], sigma[8]

  H_f = matrix.sqr((
    -1/s1, 0, 0,
    0, -1/s2, 0,
    0, 0, -1/s3))

  H_fx = H_f + matrix.sqr((
    -1/peak_fx[0]**2, 0, 0,
    0, 0, 0,
    0, 0, 0))

  H_fy = H_f + matrix.sqr((
    0, 0, 0,
    0, -1/peak_fy[1]**2, 0,
    0, 0, 0))

  H_fz = H_f + matrix.sqr((
    0, 0, 0,
    0, 0, 0,
    0, 0, -1/peak_fz[2]**2))

  E_f = exp((-0.5 * (peak_f-mu).transpose() * sigma_inv * (peak_f-mu))[0])
  E_fx = peak_fx[0] * exp((-0.5 * (peak_fx-mu).transpose() * sigma_inv * (peak_fx-mu))[0])
  E_fy = peak_fy[1] * exp((-0.5 * (peak_fy-mu).transpose() * sigma_inv * (peak_fy-mu))[0])
  E_fz = peak_fz[2] * exp((-0.5 * (peak_fz-mu).transpose() * sigma_inv * (peak_fz-mu))[0])

  A = (1 / sqrt((-H_f).determinant())) * E_f
  B = (1 / sqrt((-H_fx).determinant())) * E_fz
  C = (1 / sqrt((-H_fy).determinant())) * E_fy
  D = (1 / sqrt((-H_fz).determinant())) * E_fz

  xc = B / A
  yc = C / A
  zc = D / A

  v = matrix.col((xc, yc, zc)).normalize()

  return v


def compute_laplace(mu, sigma):


  theta_f, phi_f = compute_peak_f(mu, sigma)
  theta_fx, phi_fx = compute_peak_fx(mu, sigma)
  theta_fy, phi_fy = compute_peak_fy(mu, sigma)
  theta_fz, phi_fz = compute_peak_fz(mu, sigma)

  # print theta_f, phi_f
  # print theta_fx, phi_fx
  # print theta_fy, phi_fy
  # print theta_fz, phi_fz

  peak_f = matrix.col((
    sin(theta_f)*cos(phi_f),
    sin(theta_f)*sin(phi_f),
    cos(theta_f)))

  peak_fx = matrix.col((
    sin(theta_fx)*cos(phi_fx),
    sin(theta_fx)*sin(phi_fx),
    cos(theta_fx)))

  peak_fy = matrix.col((
    sin(theta_fy)*cos(phi_fy),
    sin(theta_fy)*sin(phi_fy),
    cos(theta_fy)))

  peak_fz = matrix.col((
    sin(theta_fz)*cos(phi_fz),
    sin(theta_fz)*sin(phi_fz),
    cos(theta_fz)))

  x = compute_mean_estimate(peak_f, peak_fx, peak_fy, peak_fz, mu, sigma)

  theta = acos(x[2])
  phi = atan2(x[1], x[0])

  return x
  # mu_phi = atan2(mu[1], mu[0])
  # mu_theta = acos(mu[2] / mu.length())
  # s_sq = matrix.col((sigma[0], sigma[4], sigma[8]))
  # print mu_theta, mu_phi, func_f(mu_theta, mu_phi, mu, s_sq)
  # print theta, phi, func_f(theta, phi, mu, s_sq)



def compute_mean_plane(mu, sigma):
  z = matrix.col((0, 0, 1))
  axis = z.cross(mu.normalize())
  angle = acos(z.dot(mu.normalize()))
  R = axis.axis_and_angle_as_r3_rotation_matrix(angle=angle, deg=False)

  sigma_1 = R.transpose() * sigma * R
  mu_1 = R.transpose() * mu

  sigma_11 = matrix.sqr((
    sigma_1[0], sigma_1[1],
    sigma_1[3], sigma_1[4]))
  sigma_12 = matrix.col((sigma_1[2], sigma_1[5]))
  sigma_21 = matrix.col((sigma_1[3], sigma_1[7])).transpose()
  sigma_22 = sigma_1[8]

  mu1 = matrix.col((mu_1[0], mu_1[1]))
  mu2 = mu_1[2]

  z = 1

  mu_new_1 = mu1 + sigma_12 * (1/sigma_22) * (z - mu2)
  v = matrix.col((
    mu_new_1[0],
    mu_new_1[1],
    z)).normalize()

  x_new = R * v

  return x_new

def compute_mean_peak(mu, sigma):
  return compute_laplace(mu, sigma)


if __name__ == '__main__':
  import json

  data = json.load(open("spots.txt"))

  mu = data['mu']
  sigma = data['sigma']
  mean_numerical = data['x']
  mean_plane = []
  mean_peak = []
  angle_plane = []
  angle_peak = []

  for i in range(len(mu)):
    x = matrix.col(mean_numerical[i])
    m = matrix.col(mu[i])
    s = matrix.sqr(sigma[i])
    x2 = compute_mean_plane(m, s)
    try:
      x3 = compute_mean_peak(m , s)
    except Exception:
      x3 = None
    mean_plane.append(x2)
    mean_peak.append(x3)

    a2 = x.angle(x2)
    angle_plane.append(a2)
    if x3 is None:
      a3 = None
    else:
      a3 = x.angle(x3)
      angle_peak.append(a3)

    print i, a2, a3

  from matplotlib import pylab
  pylab.hist(angle_plane, bins=50)
  pylab.show()
  pylab.hist(angle_peak, bins=50)
  pylab.show()
