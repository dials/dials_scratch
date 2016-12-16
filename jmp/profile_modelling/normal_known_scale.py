


def compute_derivatives(A, B, N, K, mu, sigma, USE):
  from math import erf, sqrt, exp, pi, log
  erf_A = erf((A[0] - mu) / (sqrt(2) * sigma))
  erf_B = erf((B[-1] - mu) / (sqrt(2) * sigma))
  exp_A = exp(-(A[0] - mu)**2 / (2 * sigma**2))
  exp_B = exp(-(B[-1] - mu)**2 / (2 * sigma**2))
  #Z = 0.5*(erf_B - erf_A)
  # if Z < 1e-100:
  #   for j in range(len(USE)):
  #     USE[j] = False
  #   print "ARG"
  #   return 0, 0, 0

  # ntot = sum(N)

  def f(x):
    return (1.0/(sqrt(2*pi)*sigma)) * exp(-(x-mu)**2 / (2*sigma**2))

  def integral(a, b):
    return (b-a)*(f(a) + f(b)) / 2.0

  from scipy.integrate import romberg
  restart = True
  while restart:
    restart = False
    sum1 = 0
    ntot = 0
    Z = 0
    use_count = 0
    bad_count = 0
    for j, (a, b, n) in enumerate(zip(A, B, N)):
      if USE[j]:
        erf_a = erf((a - mu) / (sqrt(2) * sigma))
        erf_b = erf((b - mu) / (sqrt(2) * sigma))
        # exp_a = exp(-(a - mu)**2 / (2 * sigma**2))
        # exp_b = exp(-(b - mu)**2 / (2 * sigma**2))
        zi = 0.5*(erf_b - erf_a)
       # zi = romberg(f, a, b)
        if zi < 1e-10:
          zi = 1e-10
          bad_count += 1
          # USE[j] = False
          # restart =True
          # print "RESTART"
          # break
        Z += zi
        sum1 += n * log(K*zi)
        use_count += 1
        ntot += n
        assert n > 0
    if use_count > 0:
      L = sum1 - K * Z
      #L = sum1 - K * log(Z)
      dL= 0
      d2L = 0
    else:
      L = 0
      dL = 0
      d2L = 0

  return -L, dL, d2L, bad_count

def compute_all_derivatives(A, B, N, I, mu, sigma, USE):
  i0 = 0
  DL = 0
  D2L = 0
  L = 0
  bad_count = 0
  for j in range(len(I)):

    i1 = i0 + I[j]
    AA = A[i0:i1]
    BB = B[i0:i1]
    NN = N[i0:i1]
    UU = USE[i0:i1]
    KK = sum(NN)

    l,dl,d2l,bc = compute_derivatives(AA, BB, NN, KK, mu, sigma, UU)
    L += l
    DL += dl
    D2L += d2l
    USE[i0:i1] = UU
    bad_count += bc

    i0 = i1
  from math import pi
  print "Num: ", sigma * 180 / pi, bad_count
  return L, DL, D2L

def estimate(A, B, N, I, mu, a, b):
  from scipy.optimize import minimize
  #print "ESTIMATE"

  USE = [True] * len(A)

  def f(sigma):
    L, DL, D2L = compute_all_derivatives(A, B, N, I, mu, sigma, USE)
    return L

  from math import sqrt, pi
  tol = 0.0001 * pi / 180.0
  gr = (sqrt(5) + 1)/2.0
  a_start = a
  b_start = b
  restart = True
  while restart == True:
    use_count= USE.count(True)
    restart = False
    a = a_start
    b = b_start
    c = b - (b - a) / gr
    d = a + (b - a) / gr
    while abs(c - d) > tol:
      if f(c) < f(d):
        b = d
      else:
        a = c
      if USE.count(True) != use_count:
       # print "RES"
        restart = True
        break

      # we recompute both c and d here to avoid loss of precision which may lead to incorrect results or infinite loop
      c = b - (b - a) / gr
      d = a + (b - a) / gr

      #print c*180/pi, d*180/pi, f(c), f(d)
    #print c, d

  return (b + a) / 2
  # for i in range(100):

  #   L, DL, D2L = func(sigma)

  #   delta = DL / abs(D2L) # abs ensures that we always go down hill

  #   from math import pi
  #   print sigma, (sigma - delta)
  #   sigma -= delta

  #   if sigma < 0.0001:
  #      sigma = 0.0001
