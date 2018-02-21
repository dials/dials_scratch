
from scitbx.linalg import eigensystem
from scitbx import matrix
from dials.array_family import flex
from math import pi, exp, sqrt, sin, atan2, cos, acos, erf
from random import uniform

def normal_1d(x, mu, sigma_sq):
  A = 1.0 / sqrt(2*pi * sigma_sq)
  B = ((x - mu)**2 / (sigma_sq))
  return A * exp(-0.5 * B)

def normal_3d(x, mu, sigma):
  A = 1.0 / sqrt((2*pi)**3 * sigma.determinant())
  B = ((x - mu).transpose() * sigma.inverse() * (x - mu))[0]
  return A * exp(-0.5 * B)

def f(x, mu1, mu2, mu3, s1, s2, s3):
  A = mu1**2 / (1 - x*s1)**2
  B = mu2**2 / (1 - x*s2)**2
  C = mu3**2 / (1 - x*s3)**2
  return A + B + C - 1

def df(x, mu1, mu2, mu3, s1, s2, s3):
  A = 2*s1 * mu1**2 / (1 - x*s1)**3
  B = 2*s2 * mu2**2 / (1 - x*s2)**3
  C = 2*s3 * mu3**2 / (1 - x*s3)**3
  return A + B + C

def d2f(x, mu1, mu2, mu3, s1, s2, s3):
  A = 6*s1 * mu1**2 / (1 - x*s1)**4
  B = 6*s2 * mu2**2 / (1 - x*s2)**4
  C = 6*s3 * mu3**2 / (1 - x*s3)**4
  return A + B + C

def compute_l0(x, mu, s):
  return (1-mu[0]/x[0])/s[0]**2

def compute_mode(mu, sigma):

  # Compute the eigen decomposition of the covariance matrix
  eigen_decomposition = eigensystem.real_symmetric(sigma.as_flex_double_matrix())
  Q = matrix.sqr(eigen_decomposition.vectors())
  L = matrix.diag(eigen_decomposition.values())

  # To simplify calculation of the mode, we change basis
  mu_prime = Q.transpose() * mu

  # Compute the initial value of the lagrange multiplier
  l0 = compute_l0(mu_prime.normalize(), mu_prime, L)

  # Compute the value via Halley's method
  while True:
    mu1, mu2, mu3 = mu_prime
    s1, s2, s3 = L[0], L[4], L[8]
    U = f(l0, mu1, mu2, mu3, s1, s2, s3)
    Up = df(l0, mu1, mu2, mu3, s1, s2, s3)
    Up2 = d2f(l0, mu1, mu2, mu3, s1, s2, s3)
    l = l0 - 2*U*Up / (2*Up**2 - U*Up2)
    if abs(l - l0) < 1e-7:
      break
    l0 = l

  # Compute the mode
  x = matrix.col((
    mu_prime[0] / (1 - l * L[0]),
    mu_prime[1] / (1 - l * L[4]),
    mu_prime[2] / (1 - l * L[8])))

  # Now transform back into the original basis
  return Q * x

def compute_plane_average(mu, sigma):

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

  def func(z):
    return z - mu2
  
  def compute_z():
    mu3 = (mu2 + 1) / 2.0
    a = 0
    b = 1
    s2 = sigma_22 / 16
    s = sqrt(s2)
    alpha = (a - mu3) / s
    beta = (b - mu3) / s
    def p(x):
      return normal_1d(x, 0, 1)
    def P(x):
      return 0.5 * (1 + erf(x/sqrt(2)))
    Z = P(beta) - P(alpha)
    return mu3 + (p(alpha) - p(beta)) * s / Z
  
  def compute_z():
    mu3 = (mu2+1) / 2.0
    a = 0
    b = 1
    s = sqrt(sigma_22) / 2
    alpha = (a - mu3) / s
    beta = (b - mu3) / s
    def p(x):
      return normal_1d(x, 0, 1)
    def P(x):
      return 0.5 * (1 + erf(x/sqrt(2)))
    Z = P(beta) - P(alpha)
    return mu3 + (p(alpha) - p(beta)) * s / Z
  
  # def compute_z():
  #   sigma_new = sigma_11 - sigma_12*(1/sigma_22)*sigma_21
  #   sigma_new_inv = sigma_new.inverse()
  #   sigma_inv = 0.5 * (sigma_new_inv[0] + sigma_new_inv[3])
  #   A = 1/sigma_22 + 0.5 * sigma_inv
  #   B = -2*(mu2 / sigma_22 + 0.5 * sigma_inv)
  #   H = -B / (2*A)
  #   S = 1/A
  #   mu3 = H
  #   a = 0
  #   b = 1
  #   s = sqrt(S)
  #   print s, sqrt(sigma_22), mu2, mu3
  #   alpha = (a - mu3) / s
  #   beta = (b - mu3) / s
  #   def p(x):
  #     return normal_1d(x, 0, 1)
  #   def P(x):
  #     return 0.5 * (1 + erf(x/sqrt(2)))
  #   Z = P(beta) - P(alpha)
  #   return mu3 + (p(alpha) - p(beta)) * s / Z

  def compute_z2():
    b = mu2/sigma_22
    e = exp(b)
    return (e*(b-1)+1) / b*(e-1)

  #a = 1
  #b = 0.0
  #n = 100
  #I = a + b
  #I2 = normal_1d(a, mu2, sigma_22) + normal_1d(b, mu2, sigma_22)
  #for i in range(1, n):
  #  z = a + (b - a)*i / n
  #  I += a * normal_1d(z, mu2, sigma_22)
  #  I2 += normal_1d(z, mu2, sigma_22)
  ##I *= abs(b - a) / n

  # z = I/I2
  # print I / I2
  z = compute_z()
  z2 = compute_z2()

  print z, mu2
  print z2

  mu_new_1 = mu1 + sigma_12 * (1/sigma_22) * (z - mu2)
  v = matrix.col((
    mu_new_1[0],
    mu_new_1[1],
    z)).normalize()
  # a = 1
  # b = 0.9
  
  # w1 = normal_1d(a, mu2, sigma_22)
  # w2 = normal_1d(b, mu2, sigma_22)
  # w1, w2 = w1 / (w1 + w2), w2 / (w1 + w2)

  
  # mu_new_1 = mu1 + sigma_12 * (1/sigma_22) * (a - mu2)
  # mu_new_2 = mu1 + sigma_12 * (1/sigma_22) * (b - mu2)

  # print w1, w2, w1+w2

  # v = matrix.col((
  #   w1 * mu_new_1[0] + w2 * mu_new_2[0],
  #   w1 * mu_new_1[0] + w2 * mu_new_2[1],
  #   w1 * a + w2 * b)).normalize()

  #v = matrix.col((mu_new[0], mu_new[1], 1)).normalize()

  x_new = R * v
  
  return x_new


def compute_peak_and_mean(mu, sigma):

  I1 = matrix.col((0, 0, 0))
  I2 = 0
  peak_value = 0
  peak_coord = None
  N = 1200

  # theta_0 = (45-30) * pi / 180
  # theta_1 = (45+30) * pi / 180
  # phi_0 = (45-30) * pi / 180
  # phi_1 = (45+30) * pi / 180
  
  theta_0 = 0
  theta_1 = pi / 2
  phi_0 =  0
  phi_1 =  pi / 2
  
  A = 1.0 / sqrt((2*pi)**3 * sigma.determinant())
  sigma_inverse = sigma.inverse()

  for j in range(N):
    for i in range(N):
      theta = theta_0 + j * (theta_1 - theta_0) / N
      phi = phi_0 + i * (phi_1 - phi_0) / N
      x = sin(theta)*cos(phi)
      y = sin(theta)*sin(phi)
      z = cos(theta)

      v = matrix.col((x, y, z))
      B = ((v - mu).transpose() * sigma_inverse * (v - mu))[0]
      g = A * exp(-0.5 * B)

      if g > peak_value:
        peak_value = g
        peak_coord = v
      
      sin_theta = sin(theta)
      I1 += g * v * sin_theta
      I2 += g * sin_theta
  # for i in range(200000):
  #   x = uniform(mu[0]-1, mu[0]+1)
  #   y = uniform(mu[1]-1, mu[1]+1)
  #   z = uniform(mu[2]-1, mu[2]+1)
  #   x = matrix.col((x, y, z)).normalize()
  #   g = normal_3d(x, mu, sigma) 
  #   theta = atan2(x[1],x[0])
  #   if g > peak_value:
  #     peak_value = g
  #     peak_coord = x
  #   I1 += g * x * sin(theta) 
  #   I2 += g * sin(theta)
  return peak_coord, (I1 / I2).normalize()



def plot(mu, sigma):
  from mpl_toolkits.mplot3d import Axes3D
  import matplotlib.pyplot as plt
  import numpy as np

  fig = plt.figure(figsize=plt.figaspect(1))  # Square figure
  ax = fig.add_subplot(111, projection='3d')

  def draw_ewald_sphere():
  
    u, v = np.mgrid[0:2*np.pi:100j, 0:np.pi:100j]
    x = np.cos(u)*np.sin(v)
    y = np.sin(u)*np.sin(v)
    z = np.cos(v)
    ax.plot_surface(x, y, z, color="r", rstride=4, cstride=4, alpha=0.2)
  
  def draw_rlp():

    # Compute the eigen decomposition of the covariance matrix
    eigen_decomposition = eigensystem.real_symmetric(sigma.as_flex_double_matrix())
    Q = matrix.sqr(eigen_decomposition.vectors())
    L = eigen_decomposition.values()
    radii = 3* L

    # now carry on with EOL's answer
    u = np.linspace(0.0, 2.0 * np.pi, 100)
    v = np.linspace(0.0, np.pi, 100)
    x = radii[0] * np.outer(np.cos(u), np.sin(v))
    y = radii[1] * np.outer(np.sin(u), np.sin(v))
    z = radii[2] * np.outer(np.ones_like(u), np.cos(v))
    for i in range(len(x)):
      for j in range(len(x)):
        [x[i,j],y[i,j],z[i,j]] = np.dot([x[i,j],y[i,j],z[i,j]],
                                        Q.as_list_of_lists()) + mu

    # plot
    ax.plot_wireframe(x, y, z,  rstride=4, cstride=4, color='b', alpha=0.2)
  ax.set_xlim3d(0.4, 0.8)
  ax.set_ylim3d(0.4, 0.8)
  ax.set_zlim3d(0.4, 0.8)
  ax.set_xlim(0.4, 0.8)
  ax.set_ylim(0.4, 0.8)
  ax.set_zlim(0.4, 0.8)
  draw_ewald_sphere()
  draw_rlp()
  plt.show()
  

if __name__ == '__main__':

  mu = matrix.col((1, 1, 1)).normalize() * 1.05
  sig = matrix.sqr((
    0.1, 0, 0,
    0, 0.14, 0,
    0, 0, 0.17))
  sigma = sig * sig.transpose()

  mode = compute_mode(mu, sigma)

  mest = compute_plane_average(mu, sigma)
  
  peak, mean = compute_peak_and_mean(mu, sigma)

  print "Mode: ", tuple(mode), mode.length()
  print "Peak: ", tuple(peak), peak.length()
  print "Mean: ", tuple(mean), mean.length()
  print "Mest: ", tuple(mest), mest.length()
  print "Vect: ", tuple(mu.normalize())

  print "Mode - Mean angle", mean.angle(mode) * 180.0 / pi
  print "Mest - Mean angle", mean.angle(mest) * 180.0 / pi
  print "Vect - Mean angle", mean.angle(mu) * 180.0 / pi


  plot(mu, sigma)
