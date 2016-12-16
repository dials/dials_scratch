
def compute_centroid(X, Y, fixed_mean=None):
  from math import sqrt
  if fixed_mean is not None:
    mean = fixed_mean
  else:
    mean = sum(xx * yy for xx, yy in zip(X, Y)) / sum(Y)
  var = sum(yy * (xx - mean)**2 for xx, yy in zip(X, Y)) / (sum(Y) - 1)
  return mean, sqrt(var)

def compute_mean_centroid(A, B, N, I):
  from math import sqrt
  MM, _ = compute_centroid((B+A)/2.0, N)
  VV = 0
  i0 = 0
  for num in I:
    i1 = i0 + num
    AA = A[i0:i1]
    BB = B[i0:i1]
    NN = N[i0:i1]
    _, sigma = compute_centroid((BB+AA)/2.0, NN, fixed_mean=0)
    VV += sigma**2
  return MM, sqrt(VV / len(I))

if __name__ == '__main__':


  def gaussian(x, mu, sigma):
    from math import pi, sqrt, exp
    return [(1.0 / (sqrt(2.0*pi) * sigma)) * exp(-(xx - mu)**2 / (2*sigma**2)) for xx in x]


  import pickle
  import sys

  A, B, N, I = pickle.load(open(sys.argv[1]))

  mean0, sigma0 = compute_centroid((B+A)/2.0, N)

  print mean0, sigma0

  mean1, sigma1 = compute_mean_centroid(A, B, N, I)

  print mean1, sigma1

  from dials.algorithms.statistics import BinnedGMMSingle1DFixedMean
  from dials.algorithms.statistics import BinnedGMMSingle1D
  from dials.array_family import flex

  result = BinnedGMMSingle1DFixedMean(
    A, B, N,
    mean0,
    sigma0,
    1e-12,
    10000)

  mean = result.mu()
  sigma = result.sigma()
  print mean, sigma
  print result.num_iter()

  # xd = X[1] - X[0]
  # area = sum(yy * xd for yy in Y)
  # Y = [yy / area for yy in Y]


  # from matplotlib import pylab
  # pylab.plot(X, Y)
  # pylab.plot(X, gaussian(X, mean0, sigma0))
  # pylab.plot(X, gaussian(X, mean, sigma))
  # pylab.show()
