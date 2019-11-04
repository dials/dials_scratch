from __future__ import print_function

from math import sqrt, sin, cos, pi, exp


def normal_2d(x, mu, sigma):
    A = 1.0 / sqrt((2 * pi) ** 2 * sigma.determinant())
    B = ((x - mu).transpose() * sigma.inverse() * (x - mu))[0]
    return A * exp(-0.5 * B)


def find_maximum(mu1, mu2, s1, s2):

    from scitbx import matrix

    mu = matrix.col((mu1, mu2))
    sigma = matrix.sqr((s1 ** 2, 0, 0, s2 ** 2))

    t = []
    f = []

    a = -pi
    b = pi

    peak_theta = 0
    peak_value = 0

    for i in range(10000):
        theta = a + i * (b - a) / 10000
        x = matrix.col((cos(theta), sin(theta)))
        ft = normal_2d(x, mu, sigma)

        if ft > peak_value:
            peak_value = ft
            peak_theta = theta

        t.append(theta * 180 / pi)
        f.append(ft)

    x = cos(peak_theta)
    y = sin(peak_theta)

    print(x, y)


mu1 = 1.1 * 1 / sqrt(2)
mu2 = 1.1 * 1 / sqrt(2)

s1 = 0.1
s2 = 0.2


def func(x, mu1, mu2, s1, s2):
    A = -(s1 ** 4) * s2 ** 4
    B = 2 * (s1 ** 2 * s2 ** 2) * (s1 ** 2 + s2 ** 2)
    C = (
        (mu2 ** 2 * s1 ** 4 + mu1 ** 2 * s2 ** 4)
        - 4 * s1 ** 2 * s2 ** 2
        - (s1 ** 4 + s2 ** 4)
    )
    D = 2 * ((s1 ** 2 + s2 ** 2) - (s1 ** 2 * mu2 ** 2 + s2 ** 2 * mu1 ** 2))
    E = mu1 ** 2 + mu2 ** 2 - 1
    return A * x ** 4 + B * x ** 3 + C * x ** 2 + D * x + E


def df(x, mu1, mu2, s1, s2):
    A = -(s1 ** 4) * s2 ** 4
    B = 2 * (s1 ** 2 * s2 ** 2) * (s1 ** 2 + s2 ** 2)
    C = (
        (mu2 ** 2 * s1 ** 4 + mu1 ** 2 * s2 ** 4)
        - 4 * s1 ** 2 * s2 ** 2
        - (s1 ** 4 + s2 ** 4)
    )
    D = 2 * ((s1 ** 2 + s2 ** 2) - (s1 ** 2 * mu2 ** 2 + s2 ** 2 * mu1 ** 2))
    return 4 * A * x ** 3 + 3 * B * x ** 2 + 2 * C * x + D


l0 = 0

for t in range(10):
    l0 = l0 - func(l0, mu1, mu2, s1, s2) / df(l0, mu1, mu2, s1, s2)

    x1 = mu1 / (1 - l0 * s1 ** 2)
    x2 = mu2 / (1 - l0 * s2 ** 2)

    print(l0, x1, x2, x1 ** 2 + x2 ** 2)

find_maximum(mu1, mu2, s1, s2)

# X = []
# Y = []
# import numpy as np
# for x in np.arange(-100, 100, 0.01):

#   X.append(x)
#   Y.append(func(x, mu1, mu2, s1, s2))

# from matplotlib import pylab
# pylab.plot(X, Y)
# pylab.show()
