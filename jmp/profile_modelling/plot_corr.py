
X1, Y1 = zip(*map(lambda x: map(float, x.split()), open("KZ.txt")))
X2, Y2 = zip(*map(lambda x: map(float, x.split()), open("KlnZ.txt")))
X3, Y3 = zip(*map(lambda x: map(float, x.split()), open("sum.txt")))

Y1 = map(lambda x: -x, Y1)
Y2 = map(lambda x: -x, Y2)

Y4 = [y3 + y1 for y1, y3 in zip(Y1, Y3)]
Y5 = [y3 + y2 for y2, y3 in zip(Y2, Y3)]
Y6 = [y2 - y1 for y1, y2 in zip(Y2, Y1)]

from matplotlib import pylab
pylab.plot(X1, Y1, color='blue')
pylab.plot(X2, Y2, color='red')
pylab.plot(X3, Y3, color='black')
pylab.plot(X3, Y4, color='green')
pylab.plot(X3, Y5, color='orange')
pylab.plot(X3, Y6, color='yellow')
pylab.show()
