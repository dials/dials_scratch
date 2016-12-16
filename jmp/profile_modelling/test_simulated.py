
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


if __name__ == '__main__':
  mu = 0
  sigma = 5.0

  A = []
  B = []
  N = []
  I = []

  for j in range(1000):
    from random import uniform, randint
    # XMIN = mu-4*sigma
    # XMAX = mu+4*sigma
    # NUM = 10
    # COUNTS = 10000
    XMIN = uniform(mu-4*sigma, mu-0.1*sigma)
    XMAX = uniform(mu+4*sigma, mu+0.1*sigma)
    NUM = randint(1, 3)
    COUNTS = uniform(100, 10000)

    AA = [XMIN + i*(XMAX-XMIN)/(NUM) for i in range(NUM)]
    BB = [XMIN + (i+1)*(XMAX-XMIN)/(NUM) for i in range(NUM)]

    NN = gaussian(AA, BB, COUNTS, mu, sigma)
    A.extend(AA)
    B.extend(BB)
    N.extend(NN)
    I.append(len(AA))



  from truncated_normal import compute_all_derivatives, estimate
  from normal import compute_all_derivatives, estimate

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

  minimum = X[Y.index(min(Y))]
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

  estimate(A, B, N, I, mu, 0.1, 10)
