from __future__ import division
from __future__ import print_function


def mpart(a):
    print("________________________________________________________________")
    print(a)
    print("id(a) =", id(a))
    a[:, :] = 0
    for pos in range(0, 4, 1):
        a[pos, pos] = pos
