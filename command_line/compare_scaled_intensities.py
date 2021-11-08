import sys

import numpy as np

from matplotlib import pyplot as plt

from dials.array_family import flex


def cc(i1, i2):
    lc = flex.linear_correlation(i1, i2)
    return lc.coefficient()


def prob(i1, i2, s1, s2):
    sq1 = flex.sqrt(s1)
    sq2 = flex.sqrt(s2)
    pp = flex.abs(i1 - i2) / (sq1 + sq2)
    # print(flex.min(pp), flex.max(pp))
    prod = 1
    for i in range(pp.size()):
        prod = prod * pp[i]
    return prod


def plot_cc(var, name):
    fig = plt.figure()
    plt.title(name)
    plt.imshow(var)
    plt.colorbar()


# TODO FIXME NB. Should I just delete the reflections that are <0. (like with a mask I guess) Just how many are they? Does that make sense?
def main(data):
    l = len(data)
    mat_I = np.full((l, l), 0.0)
    mat_s = np.full((l, l), 0.0)

    for a in range(l):
        mat_I[(a, a)] = 1.0
        mat_s[(a, a)] = 1.0
        for b in range(a + 1, l):
            m12 = data[a].match(data[b])
            r1 = data[a].select(m12[0])
            r2 = data[b].select(m12[1])
            I = cc(r1["intensity.scale.value"], r2["intensity.scale.value"])
            mat_I[(a, b)] = mat_I[(b, a)] = I
            s = cc(r1["intensity.scale.variance"], r2["intensity.scale.variance"])
            mat_s[(a, b)] = mat_s[(b, a)] = s
            p = prob(
                r1["intensity.scale.value"],
                r2["intensity.scale.value"],
                r1["intensity.scale.variance"],
                r2["intensity.scale.variance"],
            )
            print(a, b, p)

    plot_cc(mat_I, "Intensity correlation coefficient")
    plot_cc(mat_s, "Intensity variance correlation coefficient")
    plt.show()


if __name__ == "__main__":
    data = [flex.reflection_table.from_file(arg) for arg in sys.argv[1:]]
    data = [d.select(d.get_flags(d.flags.scaled)) for d in data]

    main(data)
