
X = []
Y1 = []
Y2 = []
Y3 = []
for line in open("result.txt").readlines():
  t = map(float, line.split())
  X.append(t[0])
  Y1.append(t[1])
  Y2.append(t[2])
  Y3.append(t[3])

from matplotlib import pylab
pylab.plot(X, Y1, label="old")
pylab.plot(X, Y2, label="naive")
pylab.plot(X, Y3, label="new")
pylab.legend()
pylab.show()
