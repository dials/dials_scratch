def generate_gaussian(xc, sx, width, scale=1000, random=False):
  from dials.array_family import flex
  from numpy.random import poisson
  from math import exp, erf, sqrt, floor
  X = flex.double(list(range(width)))

  if sx > 0:
    ax = 0.5 * flex.double(map(erf, (X - xc) / (sqrt(2.0)*sx)))
    bx = 0.5 * flex.double(map(erf, (X+1 - xc) / (sqrt(2.0)*sx)))
    xcomp = bx - ax
  else:
    x = int(floor(xc))
    xcomp = flex.double(width, 0)
    xcomp[x] = 1

  data = scale * xcomp
  if random == True:
    data = flex.double(map(poisson, data))
  return data

def generate_rect(xc, sx, width, scale=1000, random=False):
  from dials.array_family import flex
  from numpy.random import poisson
  from math import exp, erf, sqrt
  if sx > 0:
    data = flex.double(width)
    x0 = xc - sx/2
    x1 = xc + sx/2
    m = 1.0 / (x1 - x0)
    c = -x0 / (x1 - x0)
    def integral(i):
      if i < x0:
        return 0
      elif i > x1:
        return 1.0
      else:
        return m * i + c
    for i in range(width):
      a = integral(i)
      b = integral(i+1)
      mean = 1000*(b - a)
      mean = max([mean,0])
      #mean = scale * exp(-(i-centre_x)**2 / (2.0 * sigma_x**2))
      data[i] = poisson(mean)
  else:
    from math import floor
    x = int(floor(xc))
    data = flex.double(width, 0)
    data[x] = scale
  return data

def compute_centroid(profile):
  '''
  Compute the centroid as normal

  '''
  from dials.array_family import flex
  width = len(profile)
  X = flex.double([list(range(width))]) + 0.5
  I = profile
  xc = flex.sum(X * profile) / flex.sum(profile)

  #W = profile / flex.sum(profile)
  #xv = flex.sum(W * (X*X)) - xc*xc
  W = profile
  xv = flex.sum(W * (X - xc)**2) / (flex.sum(W)-1)
  xv1 = xv
  xv2 = xv / flex.sum(profile)

  xv3 = flex.sum(W * ((X - xc)**2 + 1/12.0)) / flex.sum(W)

  xv3 = (1.0/12.0)*(1.0/flex.sum(profile)**2)*flex.sum(profile**2)

  return xc, xv, xv2, xv2 + xv3

def generate_and_compute_centroid(xc, sx, width, random=True):

  profile = generate_gaussian(xc, sx, width, random=random)


  return compute_centroid(profile)


def compute_bias(sigma):
  if sigma <= 0:
    return 1.0/12.0
  def func2(C, S, N=1000):
    from math import sqrt, erf
    sum1 = 0
    sum2 = N * erf((N+1-C)/(sqrt(2)*S))
    sum3 = N * erf((-N-C)/(sqrt(2)*S))
    for i in range(-N, N):
      x = (i+1-C)/(sqrt(2)*S)
      sum1 += (-1)*erf(x)
    return 0.5 * (sum1 + sum2 + sum3)
  C = []
  V = []
  for i in range(100):
    c = 0.01*i
    v = (func2(c, sigma) - c + 0.5)**2
    C.append(c)
    V.append(v)

  b = 1.0
  a = 0.0
  n = len(V)-2
  sum1 = sum(V[1:-1])
  return ((b - a) / n)*(V[0]/2.0 + sum1 + V[-1]/2.0)

def compute_variance(sigma, num, width, random=True):
  from math import floor, sqrt
  from random import uniform
  xt_list = []
  xc_list = []
  xv_list = []
  xv0_list = []
  xv2_list = []

  for i in range(num):
    x = uniform(-2, 3)

    xt = int(floor(width / 2)) + float(x)

    xc, xv, xv2, xv3 = generate_and_compute_centroid(
      xt, sigma, width, random=random)

    xt_list.append(xt)
    xc_list.append(xc)
    xv0_list.append(xv)
    xv_list.append(xv2)
    xv2_list.append(xv3)

  mu_obs = sum(xc_list) / len(xc_list)
  var_obs = sum((xc - xt)**2 for xc, xt in zip(xc_list, xt_list)) / len(xt_list)
  var_obs2 = sum((xc - mu_obs)**2 for xc in xc_list) / len(xc_list)

  sigma_cal = sqrt((sum(xv0_list)/len(xv0_list)))

  #var_cal = sum(compute_bias(sqrt(v)) for v in xv_list) / len(xv_list)
  var_cal = 1.0/12.0 + sum(xv_list) / len(xc_list)
  var_cal2 = compute_bias(sigma) + sum(xv_list) / len(xv_list)
  var_cal3 = sum(xv2_list) / len(xv2_list)
  #var_cal2 = compute_bias(sigma_cal) + sum(xv_list) / len(xv_list)
  #var_cal = sum(xv_list) / len(xv_list)
  print sigma, sigma_cal, var_obs, var_cal, var_cal2, var_cal3

  return var_obs, var_cal, var_cal2, var_cal3


if __name__ == '__main__':

  random = True
  num = 1000
  width = 30

  sigma_list = []
  var_obs_list = []
  var_cal1_list = []
  var_cal2_list = []
  var_cal3_list = []

  for sigma in [ 0.1 * i for i in range(50) ]:

    var_obs, var_cal1, var_cal2, var_cal3 = compute_variance(sigma, num, width, random=random)

    sigma_list.append(sigma)
    var_obs_list.append(var_obs)
    var_cal1_list.append(var_cal1)
    var_cal2_list.append(var_cal2)
    var_cal3_list.append(var_cal3)

  from matplotlib import pylab
  figure = pylab.figure(figsize=(10,10/1.6), dpi=600)
  pylab.plot(sigma_list, var_obs_list, color='blue', label='obs')
  pylab.plot(sigma_list, var_cal1_list, color='red', label='cal(original)')
  pylab.plot(sigma_list, var_cal2_list, color='green', label='cal(with bias)')
  pylab.plot(sigma_list, var_cal3_list, color='yellow', label='cal(david)')
  pylab.xlabel("Sigma")
  pylab.ylabel("Var X")
  pylab.legend()
  pylab.savefig("variance_vs_sigma.png")
  #pylab.show()
