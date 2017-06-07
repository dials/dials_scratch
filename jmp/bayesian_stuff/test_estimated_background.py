from __future__ import division

from numpy.random import seed
seed(100)



def gen_shapes(N):
  from math import exp
  mid = N / 2
  sigma = (N / 2)/5
  bg = [1] * N
  fg = [0] * N
  mk = [False] * N
  for i in range(N):
    fg[i] = exp(-(i-mid)**2 / (2*sigma**2))
    mk[i] = abs(i - mid) < 3*sigma
  bg = [b / sum(bg) for b in bg]
  fg = [f / sum(fg) for f in fg]
  return bg, fg, mk

def simulate(B, S, bg, fg):
  from numpy.random import poisson
  data = [poisson(B * b + S * s) for b, s in zip(bg, fg)]
  return data

def estimate(data, background, mk):
  from math import factorial
  N = mk.count(True)
  M = mk.count(False)
  Cb = 0
  Cf = 0
  for i in range(len(data)):
    if mk[i]:
      Cf += data[i]
    else:
      Cb += data[i]
  s1 = 0
  for j in range(Cf+1):
    s1 += (1 + M/N)**j * factorial(Cf + Cb - j) / factorial(Cf - j)
  s2 = 0
  for i in range(Cf+1):
    Ci = ((1 + M/N)**i * factorial(Cf + Cb - i) / factorial(Cf - i))# / s1
    s2 += Ci * (i+1)
  return s2 / (N * s1)


# def plot_weights(data, background, mk):
#   from math import exp, factorial
#   N = mk.count(True)
#   M = mk.count(False)
#   Cb = 0
#   Cf = 0

#   X = []
#   Y = []

#   for i in range(len(data)):
#     if mk[i]:
#       Cf += data[i]
#     else:
#       Cb += data[i]
#   for j in range(Cf+1):
#     a = (1 + M/N)**(j-Cf)
#     b = factorial(Cf + Cb - j) / factorial(Cf - j)
#     jj = a
#     jj = a*b
#     X.append(j)
#     Y.append(jj)

#   from matplotlib import pylab
#   pylab.plot(X,Y)
#   pylab.show()

def Ps(data, background, mk, s):
  from math import exp, factorial
  N = mk.count(True)
  M = mk.count(False)
  Cb = 0
  Cf = 0

  for i in range(len(data)):
    if mk[i]:
      Cf += data[i]
    else:
      Cb += data[i]
  s1 = 0
  for j in range(Cf+1):
    jj = (1 + M/N)**j * factorial(Cf + Cb - j) / factorial(Cf - j)
    s1+= jj

  s2 = 0
  for i in range(Cf+1):
    Ci = ((1 + M/N)**i * factorial(Cf + Cb - i) / factorial(Cf - i)) / s1
    s2 += Ci * N*(s*N)**i * exp(-s*N) / factorial(i)
  return s2

# def Ps2(data, background, mk, s):
#   from math import exp, factorial
#   N = mk.count(True)
#   M = mk.count(False)
#   Cb = 0
#   Cf = 0
#   for i in range(len(data)):
#     if mk[i]:
#       Cf += data[i]
#     else:
#       Cb += data[i]
#   s1 = 0
#   for k in range(Cf+1):
#     a = factorial(Cf) / (factorial(k)*factorial(Cf-k))
#     b = s**k * exp(-s*N)
#     c = (N+M)**(-(Cf+Cb-k)-1) * factorial(Cf+Cb-k)
#     s1 += a * b * c
#   s2 = 0
#   for i in range(Cf+1):
#     a = factorial(Cf) / (factorial(k)*factorial(Cf-k))
#     b = N**(-k-1) * factorial(k)
#     c = (N+M)**(-(Cf+Cb-k)-1) * factorial(Cf+Cb-k)
#     s2 += a * b * c
#   return (s1 / s2)


bg, fg, mk = gen_shapes(20)

B = 100.0
S = 1.0

data = simulate(B, S, bg, fg)
background = [B * bg[i] for i in range(len(data))]


# plot_weights(data, background, mk)
# exit(0)

Nb = mk.count(False)
Nf = mk.count(True)
Cb = 0
Cf = 0
for i in range(len(data)):
  if mk[i]:
    Cf += data[i]
  else:
    Cb += data[i]

print Cb, Cf, Nb, Nf

X = []
Y1 = []
Y2 = []





X = []
Y1 = []

print Cf, Cb, (Cf - Cb) / Nb

from math import log
for i in range(100):
  s = -1.6 + i * 0.1
  p1 = Ps(data, background, mk, s)
  X.append(s)
  Y1.append(p1)


print estimate(data, background, mk) * mk.count(True)

from matplotlib import pylab
pylab.plot(X, Y1)
pylab.show()
