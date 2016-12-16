import math

def func(a, b, mu, sigma):
  from math import erf, exp, sqrt, pi, log, log10, floor
  e1 = erf((b-mu)/((sqrt(2.0))*sigma));
  e2 = erf((a-mu)/((sqrt(2.0))*sigma));
  e3 = exp(-(a-mu)*(a-mu) / ((2) * sigma*sigma)) / ((sqrt(2.0*pi))*sigma);
  e4 = exp(-(b-mu)*(b-mu) / ((2) * sigma*sigma)) / ((sqrt(2.0*pi))*sigma);
  exp0_ = (0.5) * (e1 - e2);
  exp1_ = (0.5)*mu*(e1 - e2) + (sigma*sigma)*(e3 - e4);
  exp2_ = (sigma*sigma / (2.0))*(e1 - e2) + sigma*sigma * ((a-mu)*e3 - (b-mu)*e4);

  # mag = floor(log10(abs(exp1_)))
  # F = 10**abs(mag)
  # print F
  # exp0_ = 0.5 * (F*e1 - F*e2)
  # exp1_ = (0.5)*mu*(F*e1 - F*e2) + (sigma*sigma)*(F*e3 - F*e4);
  # exp2_ = (sigma*sigma / (2.0))*(F*e1 - F*e2) + sigma*sigma * (F*(a-mu)*e3 - F*(b-mu)*e4);
  # print e1, e2, exp0_, exp1_, exp2_


  if exp0_ < 1e-16:
    exp0_ = 1e-16


  exp3 = exp1_ / exp0_
  exp4 = exp2_ / exp0_
  
  # exp3 = mu+ 1.0 + 2.0*sigma*sigma*(e3 - e4) / (e1 - e2 + (b-a)*1e-10)
  #exp3 = mu+exp(abs(log(2.0*sigma*sigma*(e3 - e4))- log((e1 - e2))))
  # exp3 = mu + 2.0*sigma*sigma*(e3 - e4) / (e1 - e2 + 1e-9) 
  # exp4 = sigma*sigma + 2.0*sigma*sigma*((a-mu)*e3 - (b-mu)*e4) / (e1 - e2 + 1e-9)

  # exp3 = mu + 2.0*sigma*sigma*(e3 - e4) / (e1 - e2 + 1e-9) 
  # exp4 = sigma*sigma + 2.0*sigma*sigma*((a-mu)*e3 - (b-mu)*e4) / (e1 - e2 + 1e-9)
  # if (exp0_ > 1e-9):
  #   exp3 = mu + 2.0*sigma*sigma*(e3 - e4) / (e1 - e2) 
  #   exp4 = sigma*sigma + 2.0*sigma*sigma*((a-mu)*e3 - (b-mu)*e4) / (e1 - e2)
  # else:
  #   d = exp(-(b-a)**2 / 2)
  #   print d
  #   exp3 = d*(b+a)/2.0
  #   exp4 = d*((b+a)/2.0)**2#sigma*sigma

  return exp0_, exp1_, exp2_, exp3, exp4


from scipy.integrate import romberg

mu = (0)
sigma = (1.0)

def func1(x):
	from math import erf, exp, sqrt, pi
	return ((1.0)/((sqrt(2*pi))*sigma)) * exp(-(x - mu)**2 / ((2.0)*sigma**2))

def func2(x):
	from math import erf, exp, sqrt, pi
	return x*(1.0/(sqrt(2*pi)*sigma)) * exp(-(x - mu)**2 / (2.0*sigma**2))

def func3(x):
	from math import erf, exp, sqrt, pi
	return (x-mu)**2 * (1.0/(sqrt(2*pi)*sigma)) * exp(-(x - mu)**2 / (2.0*sigma**2))

X = []
Y1 = []
Y2 = []
Y3 = []
Z1 = []
Z2 = []
Z3 = []
Q2 = []
Q3 = []
N = 60
R = 20
for i in range(N):
  a = ((-R) * sigma + sigma * (i*2*R / float(N)))
  b = ((-R) * sigma + sigma * ((i+1)*2*R / float(N)))
  x = (b+a)/(2.0)
  y1, y2, y3, y4, y5 = func(a, b, mu, sigma)

  # z1 = romberg(func1, a, b)
  # z2 = romberg(func2, a, b)
  # z3 = romberg(func3, a, b)

  # print z1

  # y4 = z2 / z1
  # y5 = z3 / z1

  X.append(x)
  # Y1.append(y1)
  # Y2.append(y2)
  # Y3.append(y3)
  Q2.append(y4)
  Q3.append(y5)
# Z1.append(z1)
# Z2.append(z2)
# Z3.append(z3)
# if y1 > 0:
# 	Q2.append(y2 / y1)
# 	Q3.append(y3 / y1)
# else:
# 	print "NOOO"
# 	Q2.append(0)	
# 	Q3.append(0)	

from math import exp
print X[1] - X[0]
d = exp(-(X[1]-X[0])**2/(16*sigma**2))
print d

from matplotlib import pylab
#pylab.plot(X, Y1)
# pylab.plot(X, Y2)
# pylab.plot(X, Y3)
# pylab.plot(X, Z1)
# pylab.plot(X, Z2)
# pylab.plot(X, Z3)
# pylab.plot(X, Q2)
# pylab.plot(X, [xx for xx in X])
pylab.plot(X, Q3, label="real")
pylab.plot(X, [d*xx*xx for xx in X], label="approx")
pylab.legend()
pylab.show()
