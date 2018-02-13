

def normal(x, mu, sigma):
  from math import sqrt, pi, exp
  return (1.0 / (sqrt(2.0*pi) * sigma)) * exp(-0.5*(x - mu)**2 / sigma**2)



from math import sqrt


X = []
Y1 = []
Y2 = []
Y3 = []
Y4 = []

mu1 = 100
mu2 = 200
sig1 = 100
sig2 = 100
X0 = -500
X1 = 500

for t in range(1000):

  x = X0 + t

  y1 = normal(x, mu1, sig1)
  y2 = normal(x, mu2, sig2)
  y3 = y1 * y2
  y4 = normal(x, 0, sqrt(sig1**2+sig2**2))

  X.append(x)
  Y1.append(y1)
  Y2.append(y2)
  Y3.append(y3)
  Y4.append(y4)


print sum(Y3), Y4[600]

from matplotlib import pylab
pylab.plot(X, Y1, color='blue')
pylab.plot(X, Y2, color='red')
pylab.plot(X, Y3, color='orange')
pylab.plot(X, Y4, color='black')
pylab.show()
