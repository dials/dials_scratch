


def compute_derivatives(A, B, N, mu, sigma, USE):
  from math import erf, sqrt, exp, pi, log
  erf_A = erf((A[0] - mu) / (sqrt(2) * sigma))
  erf_B = erf((B[-1] - mu) / (sqrt(2) * sigma))
  exp_A = exp(-(A[0] - mu)**2 / (2 * sigma**2))
  exp_B = exp(-(B[-1] - mu)**2 / (2 * sigma**2))
  Z = erf_B - erf_A
  if Z < 1e-10:
    print "ARG"
    return 0, 0, 0

  ntot = sum(N)

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
        kt = sum(N)

        zi = erf_b - erf_a

        if zi < 1e-10:
          print "NOOO"
          zi = 1e-10

        A1 = sqrt(2/pi)*(((b-mu)*exp_b - (a-mu)*exp_a)**2 / sigma**4) / zi**2
        A2 = (2.0*((b-mu)*exp_b - (a-mu)*exp_a)/(sigma**3)) / zi
        A3 = (((b-mu)**3 * exp_b - (a-mu)**3 * exp_a) / sigma**5) / zi

        sum_lnvi += n * (log(kt) + log(zi) - log(erf_B - erf_A))
        dsum_lnvi += n* sqrt(2/pi)*(1/sigma**2)*((b-mu)*exp_b - (a-mu)*exp_a)/zi
        dsum2 += n * sqrt(2/pi) * (A1 - A2 + A3)

  A1 = sqrt(2/pi)*(((B[-1]-mu)*exp_B - (A[0]-mu)*exp_A)**2 / sigma**4) / Z**2
  A2 = (2.0*((B[-1]-mu)*exp_B - (A[0]-mu)*exp_A)/(sigma**3)) / Z
  A3 = (((B[-1]-mu)**3 * exp_B - (A[0]-mu)**3 * exp_A) / sigma**5) / Z

  dsum_tot2 = ntot * sqrt(2/pi) * (A1 - A2 + A3)
  dsum_tot = ntot* sqrt(2/pi)*(1/sigma**2)*((B[-1]-mu)*exp_B - (A[0]-mu)*exp_A)/Z
  L = sum_lnvi - kt
  dL = dsum_lnvi - dsum_tot
  d2L = dsum2 - dsum_tot2

  return -L, dL, d2L

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

    l,dl,d2l = compute_derivatives(AA, BB, NN, mu, sigma, UU)
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

    delta = DL / abs(D2L) # abs ensures that we always go down hill

    from math import pi
    print sigma, (sigma - delta)
    sigma -= delta

    if sigma < 0.0001:
       sigma = 0.0001
