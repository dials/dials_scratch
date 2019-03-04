from __future__ import print_function
from math import sqrt, log

def f(x, mu, s):
  A = mu
  B = (1 - s**2)
  C = -(1+mu)
  D = -s**2
  return A*x**3 + B*x**2 + C*x + D

def df(x, mu, s):
  A = mu
  B = (1 - s**2)
  C = -(1+mu)
  return 3*A*x**2 + 2*B*x + C

def d2f(x, mu, s):
  A = mu
  B = (1 - s**2)
  return 6*A*x + 2*B

def compute_peak(mu, s):
  l0 = 0
  while True:
    U = f(l0, mu, s)
    Up = df(l0, mu, s)
    Up2 = d2f(l0, mu, s)
    l = l0 - 2*U*Up / (2*Up**2 - U*Up2)
    if abs(l - l0) < 1e-7:
      break
    l0 = l

  x = ((1 + 2*l)+s**2)/mu
  return x


def min_index(x):
  m = x[0]
  i = 0
  for j in range(1, len(x)):
    if x[j] < m:
      i = j
      m = x[j]
  return i

m = 1
s = 0.1

X = []
Y = []

for i in range(200000):
  x = 0.00001 + i * 0.00001
  y = 0.5*(x - m)**2 - s**2 * log(x)

  X.append(x)
  Y.append(y)

i = min_index(Y)
print(X[i])

# print compute_peak(m, s)

# from matplotlib import pylab
# pylab.plot(X, Y)
# pylab.show()
