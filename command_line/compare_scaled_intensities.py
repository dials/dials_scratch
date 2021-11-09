import sys
import math

import numpy as np

from scipy.stats import t
from matplotlib import pyplot as plt

from dials.array_family import flex


def cc(i1, i2):
    lc = flex.linear_correlation(i1, i2)
    return lc.coefficient()


def compare(i1, i2, s1, s2):
    # |i1 - i2| <= beta*s
    beta = 2
    # beta = math.sqrt(2)

    # Calculate difference
    diff = i1 - i2
    diff = flex.abs(diff)

    # Calculate uncertainty s
    # s = sqrt(s1^2 + s2^2)
    s = flex.sqrt(s1 ** 2 + s2 ** 2)

    c = diff <= s * beta
    return c.as_numpy_array()


def chi_square():
    pass


def pseudo_chi(i1, i2, s1, s2):
    # This should be around 1
    n = i1.size()
    somm = 0
    for k in range(n):
        _d = i1[k] ** 2 + i2[k] ** 2
        somm += abs(_d) / (s1[k] ** 2 + s2[k] ** 2)
    return somm / n


def t_test_paired(i1, i2):
    # Applied to paired samples, 2-tailed
    assert i1.size() == i2.size()
    n = i1.size()

    # Degrees of freedom
    dof = n - 1

    # Calculate means
    m1 = flex.mean(i1)
    m2 = flex.mean(i2)

    # Average difference between observations
    d_avg = sum([(i1[i] - i2[i]) for i in range(n)]) / n

    # Calculate standard error of the average difference
    somm = 0
    for j in range(n):
        _d = i1[j] - i2[j]
        somm += (_d - d_avg) ** 2
    std = math.sqrt(somm / (n - 1))
    sed = std / math.sqrt(n)

    # Find t
    T = (m1 - m2) / sed

    # Considering 5% ...
    alpha = 0.05

    # Find p-value (for 2 tailed test)
    p_val = t.sf(abs(T), dof) * 2

    # If p > alpha => Null hp accepted => equal observations
    res = p_val > alpha
    return T, p_val, res


def plot_cc(var, name):
    fig = plt.figure()
    plt.title(name)
    plt.imshow(var)
    plt.colorbar()


def main(data):
    l = len(data)
    mat_I = np.full((l, l), 0.0)
    mat_s = np.full((l, l), 0.0)
    T_test = []

    for a in range(l):
        mat_I[(a, a)] = 1.0
        mat_s[(a, a)] = 1.0
        for b in range(a + 1, l):
            m12 = data[a].match(data[b])
            r1 = data[a].select(m12[0])
            r2 = data[b].select(m12[1])
            # Compute correlation coefficient for I and sigma
            I = cc(r1["intensity.scale.value"], r2["intensity.scale.value"])
            mat_I[(a, b)] = mat_I[(b, a)] = I
            s = cc(r1["intensity.scale.variance"], r2["intensity.scale.variance"])
            mat_s[(a, b)] = mat_s[(b, a)] = s
            # t-test for paired samples
            T, p_val, t_res = t_test_paired(
                r1["intensity.scale.value"], r2["intensity.scale.value"]
            )
            c_res = compare(
                r1["intensity.scale.value"],
                r2["intensity.scale.value"],
                r1["intensity.scale.variance"],
                r2["intensity.scale.variance"],
            )
            chi = pseudo_chi(
                r1["intensity.scale.value"],
                r2["intensity.scale.value"],
                r1["intensity.scale.variance"],
                r2["intensity.scale.variance"],
            )
            T_test.append(t_res)
            # print(a, b, t_res, c_res.count(False), chi)
            # what's an acceptable value for chi? 1.2 ok, but 1.5? still ok?

    t_fail = [i for i in T_test if i == False]
    print(
        f"Number of matches: {len(T_test)}, of which {len(t_fail)} fail the T-test."
    )  # It might be interesting to know which and what the lc is there...

    plot_cc(mat_I, "Intensity correlation coefficient")
    plot_cc(mat_s, "Intensity variance correlation coefficient")
    plt.show()


if __name__ == "__main__":
    data = [flex.reflection_table.from_file(arg) for arg in sys.argv[1:]]
    data = [d.select(d.get_flags(d.flags.scaled)) for d in data]

    main(data)
