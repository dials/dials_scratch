

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
    sum_p = 0
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
        sum_p += p

  ptot = 0.5 * (erf_B - erf_A)
  nrem = vtot - sum(N)
  prem = 1.0 - ptot
  vi_rem = vtot * prem
  sum_lnvi += nrem * log(vi_rem)
  #L = sum_lnvi - sum(N)*log(vtot * ptot)
  L = sum_lnvi - vtot*log(vtot)
  dL = dsum_lnvi - dvtot

  return vtot, dvtot, sum_lnvi, dsum_lnvi, L, dL

def compute_all_derivatives(A, B, N, I, mu, sigma, USE):
  i0 = 0
  DL = 0
  L = 0
  for j in range(len(I)):

    i1 = i0 + I[j]
    AA = A[i0:i1]
    BB = B[i0:i1]
    NN = N[i0:i1]
    UU = USE[i0:i1]

    v,dv,s,ds,l,dl = compute_derivatives(AA, BB, NN, mu, sigma, UU)
    L += l
    DL += dl
    USE[i0:i1] = UU

    i0 = i1
  return L, DL

def estimate(A, B, N, I, mu, sigma):
  from scipy.optimize import minimize

  USE = [True] * len(A)

  def func(sigma):
    L, DL = compute_all_derivatives(A, B, N, I, mu, sigma, USE)


  # for i in range(100):

  #   DL = func(sigma)


  #   sigma += 0.001 * DL

  # print minimize(func, sigma)

  # print DL, sigma

if __name__ == '__main__':
  mu = 0
  sigma = 5.0

  A = []
  B = []
  N = []
  I = []

  for j in range(1):
    from random import uniform, randint
    XMIN = mu-1*sigma
    XMAX = mu+1*sigma
    NUM = 10
    COUNTS = 10000
    # XMIN = uniform(mu-3*sigma, mu-0.1*sigma)
    # XMAX = uniform(mu+3*sigma, mu+0.1*sigma)
    # NUM = randint(1, 10)
    # COUNTS = uniform(100, 10000)

    AA = [XMIN + i*(XMAX-XMIN)/(NUM) for i in range(NUM)]
    BB = [XMIN + (i+1)*(XMAX-XMIN)/(NUM) for i in range(NUM)]

    NN = gaussian(AA, BB, COUNTS, mu, sigma)
    A.extend(AA)
    B.extend(BB)
    N.extend(NN)
    I.append(len(AA))



  # estimate(A, B, N, I, mu, 5.0)


  USE = [True] * len(A)

  X = []
  Y = []
  DY = []

  min_sigma = 2.0
  max_sigma = 8.0
  num = 1000
  for j in range(num):

    sigma = min_sigma + j * (max_sigma - min_sigma) / (num - 1)

    X.append(sigma)


    ntot = sum(N)

    L, dL = compute_all_derivatives(A, B, N, I, mu, sigma, USE)

    # Y.append(vtot)
    # DY.append(dvtot)
    # Y.append(sum_lnvi)
    # DY.append(dsum_lnvi)
    Y.append(L)
    DY.append(dL)

  from numpy import gradient

  DY2 = gradient(Y, gradient(X))

  from matplotlib import pylab

  D = [abs(dy) for dy in DY]
  minimum = X[D.index(min(D))]
  print minimum

  pylab.plot(N)
  pylab.ylim((0, max(N)))
  pylab.show()

  pylab.plot(X, Y)
  pylab.axvline(minimum)
  # pylab.plot(X, DY)
  # pylab.plot(X, DY2)
  pylab.show()
