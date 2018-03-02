
from matplotlib import pylab
from math import sqrt, exp
X = []
Y = []
s = 0.02
for i in range(1001):
  x = 0.9+i *0.001 * 0.1
  w = exp(-0.5 * (1-x**2) / s**2)
  # w = 1/sqrt(1 - x**2)
  # w = w * s
  X.append(x)
  Y.append(w)

pylab.plot(X, Y)
pylab.show()
