import sys

import numpy as np

from matplotlib import pyplot as plt

from dials.array_family import flex


def cc(i1, i2):
    # def cc(r1, r2, val):
    #     m12 = r1.match(r2)

    #     r1 = r1.select(m12[0])
    #     r2 = r2.select(m12[1])

    #     i1 = r1[val]
    #     i2 = r2[val]

    lc = flex.linear_correlation(i1, i2)
    return lc.coefficient()


def prob(i1, i2, s1, s2):
    sq1 = flex.sqrt(s1)
    sq2 = flex.sqrt(s2)
    pp = (i1 - i2) / (sq1 + sq2)
    prod = 1
    for i in range(pp.size()):
        prod = prod * pp[i]
    return prod


def plot_cc(var):
    fig = plt.figure()
    plt.imshow(var)


#     plt.show()


def main(data):
    l = len(data)
    mat_I = np.full((l, l), 0.0)
    mat_s = np.full((l, l), 0.0)

    for a in range(l):
        for b in range(a + 1, l):
            m12 = data[a].match(data[b])
            r1 = data[a].select(m12[0])
            r2 = data[b].select(m12[1])
            I = cc(r1["intensity.scale.value"], r2["intensity.scale.value"])
            # I = cc(data[a], data[b], "intensity.scale.value")
            mat_I[(a, b)] = mat_I[(b, a)] = I
            s = cc(r1["intensity.scale.variance"], r2["intensity.scale.variance"])
            # s = cc(data[a], data[b], "intensity.scale.variance")
            mat_s[(a, b)] = mat_s[(b, a)] = s
            p = prob(
                r1["intensity.scale.value"],
                r2["intensity.scale.value"],
                r1["intensity.scale.variance"],
                r2["intensity.scale.variance"],
            )
            print(a, b, p)

    plot_cc(mat_I)
    plot_cc(mat_s)
    plt.show()


if __name__ == "__main__":
    data = [flex.reflection_table.from_file(arg) for arg in sys.argv[1:]]
    data = [d.select(d.get_flags(d.flags.scaled)) for d in data]

    main(data)
