
def compute_centroid(X, Y):
  from math import sqrt
  mean = sum(xx * yy for xx, yy in zip(X, Y)) / sum(Y)
  var = sum(yy * (xx - mean)**2 for xx, yy in zip(X, Y)) / (sum(Y) - 1)
  return mean, sqrt(var)


if __name__ == '__main__':

  
  def gaussian(x, mu, sigma):
    from math import pi, sqrt, exp
    return [(1.0 / (sqrt(2.0*pi) * sigma)) * exp(-(xx - mu)**2 / (2*sigma**2)) for xx in x]


  import pickle
  import sys

  X, Y = pickle.load(open(sys.argv[1]))
  
  mean0, sigma0 = compute_centroid(X, Y)
  
  print mean0, sigma0

  from dials.algorithms.statistics import BinnedGMMSingle1DFixedMean 
  from dials.algorithms.statistics import BinnedGMMSingle1D
  from dials.array_family import flex
  a_list = flex.double(X)
  b_list = flex.double(X) + X[1] - X[0]
  n_list = flex.double(Y)

  result = BinnedGMMSingle1DFixedMean(
    a_list, 
    b_list, 
    n_list, 
    mean0, 
    sigma0, 
    1e-7,
    100)

  mean = result.mu()
  sigma = result.sigma()
  print mean, sigma
  print result.num_iter()

  xd = X[1] - X[0]
  area = sum(yy * xd for yy in Y)
  Y = [yy / area for yy in Y]


  from matplotlib import pylab
  pylab.plot(X, Y)
  pylab.plot(X, gaussian(X, mean0, sigma0))
  pylab.plot(X, gaussian(X, mean, sigma))
  pylab.show()
