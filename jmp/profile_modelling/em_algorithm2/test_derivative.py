

def gaussian(A, B, N, mu, sigma):
  from math import sqrt, erf
  from numpy.random import poisson
  G = []
  for a, b in zip(A,B):
    G.append(
      poisson(
        N*0.5*(
          erf((b-mu)/(sqrt(2)*sigma)) -
          erf((a-mu)/(sqrt(2)*sigma)))))
  return G


def compute_derivatives(A, B, N, mu, sigma, USE):
  from math import erf, sqrt, exp, pi, log
  erf_A = erf((A[0] - mu) / (sqrt(2) * sigma))
  erf_B = erf((B[-1] - mu) / (sqrt(2) * sigma))
  exp_A = exp(-(A[0] - mu)**2 / (2 * sigma**2))
  exp_B = exp(-(B[-1] - mu)**2 / (2 * sigma**2))

  ntot = sum(N)
  vtot = ntot / (0.5 * (erf_B - erf_A))
  dvtot = 2.0*ntot*(sqrt(2/pi)/sigma**2)*((B[-1]-mu)*exp_B - (A[0]-mu)*exp_A) / (erf_B - erf_A)**2

  restart = True
  while restart:
    sum_lnvi = 0
    dsum_lnvi = 0
    restart = False
    for j, (a, b, n) in enumerate(zip(A, B, N)):
      if USE[j]:
        erf_a = erf((a - mu) / (sqrt(2) * sigma))
        erf_b = erf((b - mu) / (sqrt(2) * sigma))
        exp_a = exp(-(a - mu)**2 / (2 * sigma**2))
        exp_b = exp(-(b - mu)**2 / (2 * sigma**2))
        p = 0.5 * (erf_b - erf_a)
        dp = -(1.0 / (sqrt(2*pi) *sigma**2))*((b-mu)*exp_b - (a-mu)*exp_a)
        vi = vtot * p
        if vi < 1e-100:
          # USE[j] = False
          # restart = True
          # break
          dlnvi = 0.0
          sum_lnvi += n * (-200)
        else:
          dvi = vtot * dp + p * dvtot
          dlnvi = dvi / vi
          dsum_lnvi += n * dlnvi
          sum_lnvi += n*log(vi)

  L = sum_lnvi - vtot
  dL = dsum_lnvi - dvtot

  return vtot, dvtot, sum_lnvi, dsum_lnvi, L, dL

def estimate(A, B, N, mu, sigma):
  USE = [True] * len(A)

  for i in range(100):
    v,dv,s,ds,l,dl = compute_derivatives(A, B, N, mu, sigma, USE)

    sigma += 0.01 * dl
    print dl, sigma

if __name__ == '__main__':
  mu = 50
  sigma = 5.0
  A = list(range(100))
  B = list(range(1,101))
  N = gaussian(A, B, 1000, mu, sigma)

  estimate(A, B, N, mu, 5.0)

  USE = [True] * len(A)

  X = []
  Y = []
  DY = []

  min_sigma = 4.0
  max_sigma = 3.0
  num = 1000
  for j in range(num):

    sigma = min_sigma + j * (max_sigma - min_sigma) / (num - 1)

    X.append(sigma)


    ntot = sum(N)

    vtot, dvtot, sum_lnvi, dsum_lnvi, L, dL = compute_derivatives(A, B, N, mu, sigma,
                                                           USE)

    # Y.append(vtot)
    # DY.append(dvtot)
    # Y.append(sum_lnvi)
    # DY.append(dsum_lnvi)
    Y.append(L)
    DY.append(dL)

  from numpy import gradient

  DY2 = gradient(Y, gradient(X))

  from matplotlib import pylab
  # pylab.plot(X, Y)
  pylab.plot(X, DY)
  pylab.plot(X, DY2)
  pylab.show()
