from __future__ import division
from __future__ import print_function


def simple_summation(data, background):
    """
  Perform a simple summation with background subtraction

  """
    signal = 0
    for C, B in zip(data, background):
        signal += C - B
    return signal


def estimate_signal(data, background, estimator="mean"):
    """
  :param data: The array of pixel data
  :param background: The array of known background values
  :param estimator: Mean or Maximum of Posterior distribution (mean, max)

  :returns: A bayesian estimate of the signal

  """
    if estimator == "mean":

        C = sum(data)
        B = sum(background)

        from scipy.special import gammaincc

        K = gammaincc(C + 1, B)
        return K * ((C + 1) * gammaincc(C + 2, B) - B * gammaincc(C + 1, B))

    elif estimator == "max":

        C = sum(data)
        B = sum(background)

        return max(0, C - B)

    else:
        raise RuntimeError("BANG!")


def gen_shapes(N):
    from math import exp

    mid = N / 2
    sigma = (N / 2) / 5
    bg = [1] * N
    fg = [0] * N
    mk = [False] * N
    for i in range(N):
        fg[i] = exp(-(i - mid) ** 2 / (2 * sigma ** 2))
        mk[i] = abs(i - mid) < 3 * sigma
    bg = [b / sum(bg) for b in bg]
    fg = [f / sum(fg) for f in fg]
    return bg, fg, mk


def simulate(B, S, bg, fg):
    from numpy.random import poisson

    data = [poisson(B * b + S * s) for b, s in zip(bg, fg)]
    return data


bg, fg, mk = gen_shapes(20)

B = 100.0
S = 1.0

Y1 = []
Y2 = []
Y3 = []

for t in range(2000):

    data = simulate(B, S, bg, fg)
    background = [B * bg[i] for i in range(len(data))]

    S_sum = simple_summation(data, background)
    S_mean = estimate_signal(data, background, estimator="mean")
    S_max = estimate_signal(data, background, estimator="max")

    print(S_sum, S_mean, S_max)

    Y1.append(S_sum)
    Y2.append(S_mean)
    Y3.append(S_max)

from matplotlib import pylab

pylab.hist(Y1, bins=20)
pylab.hist(Y2, bins=20)
# pylab.hist(Y3, bins=20)
pylab.show()
