


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
    restart = False
    sum1 = 0

    for j, (a, b, n) in enumerate(zip(A, B, N)):
      if USE[j]:
        erf_a = erf((a - mu) / (sqrt(2) * sigma))
        erf_b = erf((b - mu) / (sqrt(2) * sigma))
        exp_a = exp(-(a - mu)**2 / (2 * sigma**2))
        exp_b = exp(-(b - mu)**2 / (2 * sigma**2))
        zi = erf_b - erf_a
        if zi < 1e-10:
          zi = 1e-10
        sum1 += n * log(0.5*zi)
    L = (sum1 - ntot * log(0.5*Z))
    dL= 0
    d2L = 0

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



