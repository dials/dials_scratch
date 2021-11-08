import sys
import math

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
    # pp = (i1 - i2) / (sq1 + sq2)
    prod = 1
    for i in range(pp.size()):
        prod = prod * pp[i]
    return prod


def t_test(i1, i2):
    # Degrees of freedom
    dof = i1.size() + i2.size() - 2

    # Calculate mean and standard error
    mv1 = flex.mean_and_variance(i1)
    m1 = mv1.mean()
    std1 = mv1.unweighted_sample_standard_deviation()
    se1 = std1 / math.sqrt(i1.size())

    mv2 = flex.mean_and_variance(i2)
    m2 = mv2.mean()
    std2 = mv2.unweighted_sample_standard_deviation()
    se2 = std2 / math.sqrt(i2.size())

    # Calculate standard error of the difference
    sed = math.sqrt(se1 ** 2 + se2 ** 2)

    # Find t
    t = (m1 - m2) / sed

    # TBC
    return t


def plot_cc(var, name):
    fig = plt.figure()
    plt.title(name)
    plt.imshow(var)
    plt.colorbar()


def main(data):
    l = len(data)
    mat_I = np.full((l, l), 0.0)
    mat_s = np.full((l, l), 0.0)

    for a in range(l):
        mat_I[(a, a)] = 1.0
        mat_s[(a, a)] = 1.0
        for b in range(a + 1, l):
            # Delete negative values
            mask1 = data[a]["intensity.scale.value"] <= 0.0
            mask2 = data[b]["intensity.scale.value"] <= 0.0
            data[a].del_selected(mask1)
            data[b].del_selected(mask2)
            m12 = data[a].match(data[b])
            r1 = data[a].select(m12[0])
            r2 = data[b].select(m12[1])
            # Compute correlation coefficient for I and sigma
            I = cc(r1["intensity.scale.value"], r2["intensity.scale.value"])
            mat_I[(a, b)] = mat_I[(b, a)] = I
            s = cc(r1["intensity.scale.variance"], r2["intensity.scale.variance"])
            mat_s[(a, b)] = mat_s[(b, a)] = s
            # p = prob(
            #     r1["intensity.scale.value"],
            #     r2["intensity.scale.value"],
            #     r1["intensity.scale.variance"],
            #     r2["intensity.scale.variance"],
            # )
            # print(a, b, p)
            t = t_test(r1["intensity.scale.value"], r2["intensity.scale.value"])
            print(a, b, t)

    plot_cc(mat_I, "Intensity correlation coefficient")
    plot_cc(mat_s, "Intensity variance correlation coefficient")
    plt.show()


if __name__ == "__main__":
    data = [flex.reflection_table.from_file(arg) for arg in sys.argv[1:]]
    data = [d.select(d.get_flags(d.flags.scaled)) for d in data]

    main(data)
