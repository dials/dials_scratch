from __future__ import division
from __future__ import print_function
from math import sqrt, exp, pi, cos, sin
from scitbx import matrix


def normal_pdf(x, mu, sigma):
    """
    Multi variate normal

    """
    N = len(mu)
    sigma_inv = sigma.inverse()
    A = 1.0 / (sqrt((2 * pi) ** N * sigma.determinant()))
    B = ((x - mu).transpose() * sigma_inv * (x - mu))[0]
    return A * exp(-0.5 * B)


def integrate_non_scaled(s0, r, sigma):

    R = s0.length()
    N = 1000

    s1 = s0 + r

    def func(theta):
        return R * normal_pdf(matrix.col((R * cos(theta), R * sin(theta))), s1, sigma)

    a = -pi
    b = pi

    I_f = func(a) / 2.0 + func(b) / 2.0
    for i in range(1, N):
        I_f += func(a + i * (b - a) / N)
    I_f *= (b - a) / N

    return I_f


def integrate_scaled(s0, r, sigma):

    l = 1.0 / s0.length()

    # sigma = l**2 * sigma
    sigma = l ** 2 * sigma

    r = r / s0.length()
    s0 = s0.normalize()
    s1 = s0 + r

    R = s0.length()
    N = 1000

    def func(theta):
        return R * normal_pdf(matrix.col((R * cos(theta), R * sin(theta))), s1, sigma)

    a = -pi
    b = pi

    I_f = func(a) / 2.0 + func(b) / 2.0
    for i in range(1, N):
        I_f += func(a + i * (b - a) / N)
    I_f *= (b - a) / N

    return I_f


def integrate_zero_non_scaled(s0, sigma):

    return integrate_non_scaled(s0, matrix.col((0, 0)), sigma)


def integrate_zero_scaled(s0, sigma):

    return integrate_scaled(s0, matrix.col((0, 0)), sigma)


def test_non_scaled(s0, r, sigma):
    return integrate_non_scaled(s0, r, sigma) / integrate_zero_non_scaled(s0, sigma)


def test_scaled(s0, r, sigma):
    return integrate_scaled(s0, r, sigma) / integrate_zero_scaled(s0, sigma)


s0 = matrix.col((10, 0))
s1 = matrix.col((1, 1)).normalize() * s0.length()
r = s1 * 1.1 - s0

sigma = matrix.sqr((0.5, 0, 0, 0.5))

print(test_non_scaled(s0, r, sigma))
print(test_scaled(s0, r, sigma))
