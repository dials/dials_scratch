

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
    dsum_tot = 0
    dsum2 = 0
    dsum_tot2 = 0
    sum_p = 0
    restart = False
    for j, (a, b, n) in enumerate(zip(A, B, N)):
      if USE[j]:
        erf_a = erf((a - mu) / (sqrt(2) * sigma))
        erf_b = erf((b - mu) / (sqrt(2) * sigma))
        exp_a = exp(-(a - mu)**2 / (2 * sigma**2))
        exp_b = exp(-(b - mu)**2 / (2 * sigma**2))
        Z = 0.5 * (erf_B - erf_A)
        kt = sum(N)

        zi = erf_b - erf_a

        if zi < 1e-10:
          zi = 1e-10

        A1 = sqrt(2/pi)*(((b-mu)*exp_b - (a-mu)*exp_a)**2 / sigma**4) / zi**2
        A2 = (2.0*((b-mu)*exp_b - (a-mu)*exp_a)/(sigma**3)) / zi
        A3 = (((b-mu)**3 * exp_b - (a-mu)**3 * exp_a) / sigma**5) / zi

        sum_lnvi += n * (log(kt) + log(zi) - log(erf_B - erf_A))
        dsum_lnvi += n* sqrt(2/pi)*(1/sigma**2)*((b-mu)*exp_b - (a-mu)*exp_a)/zi
        dsum2 += n * sqrt(2/pi) * (A1 - A2 + A3)

  zi = erf_B - erf_A
  A1 = sqrt(2/pi)*(((B[-1]-mu)*exp_B - (A[0]-mu)*exp_A)**2 / sigma**4) / zi**2
  A2 = (2.0*((B[-1]-mu)*exp_B - (A[0]-mu)*exp_A)/(sigma**3)) / zi
  A3 = (((B[-1]-mu)**3 * exp_B - (A[0]-mu)**3 * exp_A) / sigma**5) / zi

  dsum_tot2 = ntot * sqrt(2/pi) * (A1 - A2 + A3)
  dsum_tot = ntot* sqrt(2/pi)*(1/sigma**2)*((B[-1]-mu)*exp_B - (A[0]-mu)*exp_A)/zi
  L = sum_lnvi - kt
  dL = dsum_lnvi - dsum_tot
  d2L = dsum2 - dsum_tot2

  return vtot, dvtot, sum_lnvi, dsum_lnvi, -L, dL, d2L

def compute_all_derivatives(A, B, N, I, mu, sigma, USE):
  i0 = 0
  DL = 0
  D2L = 0
  L = 0
  for j in range(len(I)):

    i1 = i0 + I[j]
    AA = A[i0:i1]
    BB = B[i0:i1]
    NN = N[i0:i1]
    UU = USE[i0:i1]

    v,dv,s,ds,l,dl,d2l = compute_derivatives(AA, BB, NN, mu, sigma, UU)
    L += l
    DL += dl
    D2L += d2l
    USE[i0:i1] = UU

    i0 = i1
  return L, DL, D2L

def estimate(A, B, N, I, mu, sigma):
  from scipy.optimize import minimize

  USE = [True] * len(A)

  def func(sigma):
    L, DL, D2L = compute_all_derivatives(A, B, N, I, mu, sigma, USE)
    return L, DL, D2L


  for i in range(100):

    L, DL, D2L = func(sigma)

    delta = DL / abs(D2L)

    print sigma, DL, D2L, DL/D2L, delta, sigma -delta
    sigma -= delta

    if sigma < 0.0001:
       sigma = 0.0001



  # print minimize(func, sigma)

  # print DL, sigma

if __name__ == '__main__':
  mu = 0
  sigma = 5.0

  A = []
  B = []
  N = []
  I = []

  for j in range(100):
    from random import uniform, randint
    # XMIN = mu-4*sigma
    # XMAX = mu+4*sigma
    # NUM = 10
    # COUNTS = 10000
    XMIN = uniform(mu-3*sigma, mu-0.1*sigma)
    XMAX = uniform(mu+3*sigma, mu+0.1*sigma)
    NUM = randint(1, 5)
    COUNTS = uniform(100, 10000)

    AA = [XMIN + i*(XMAX-XMIN)/(NUM) for i in range(NUM)]
    BB = [XMIN + (i+1)*(XMAX-XMIN)/(NUM) for i in range(NUM)]

    NN = gaussian(AA, BB, COUNTS, mu, sigma)
    A.extend(AA)
    B.extend(BB)
    N.extend(NN)
    I.append(len(AA))





  USE = [True] * len(A)

  X = []
  Y = []
  DY = []
  D2Y = []

  min_sigma = 2.0
  max_sigma = 20.0
  num = 1000
  for j in range(num):

    sigma = min_sigma + j * (max_sigma - min_sigma) / (num - 1)

    X.append(sigma)


    ntot = sum(N)

    L, dL, d2L = compute_all_derivatives(A, B, N, I, mu, sigma, USE)

    # Y.append(vtot)
    # DY.append(dvtot)
    # Y.append(sum_lnvi)
    # DY.append(dsum_lnvi)
    Y.append(L)
    DY.append(dL)
    D2Y.append(d2L)

  print USE.count(False)

  from numpy import gradient

  DY2 = gradient(Y, gradient(X))
  D2Y2 = gradient(DY, gradient(X))

  from matplotlib import pylab

  D = [abs(dy) for dy in DY]
  minimum = X[D.index(min(D))]
  print minimum

  # pylab.plot(N)
  # pylab.ylim((0, max(N)))
  # pylab.show()

  pylab.plot(X, Y, color='black')
  pylab.axvline(minimum)
  pylab.plot(X, DY, color='blue')
  # pylab.plot(X, DY2, color='orange')
  pylab.plot(X, D2Y, color='red')
  # pylab.plot(X, D2Y2, color='purple')
  pylab.show()

  estimate(A, B, N, I, mu, 0.1)
