import math
import argparse

import numpy as np
import pandas as pd

from pathlib import Path
from tabulate import tabulate
from scipy.stats import t, probplot

from matplotlib import pyplot as plt

from dials.array_family import flex

parser = argparse.ArgumentParser(description="Compare scaled intensities.")
parser.add_argument("refl_files", nargs="*")
parser.add_argument(
    "-o",
    "--output",
    help="Output directory to save results. If not passed, defaults to working directory.",
    type=str,
)

# Correlation coefficient
def lin_cc(x1, x2):
    lc = flex.linear_correlation(x1, x2)
    return lc.coefficient()


# Normal probability plot
def normal_probability_plot(x1, x2, v1, v2, filename):
    k = sum(x1 * x2) / sum(x2 ** 2)
    n = x1.size()
    dm_real = (x1 - k * x2) / flex.sqrt(v1 + k ** 2 * v2)
    dm_real_sorted = np.sort(dm_real.as_numpy_array())

    res = probplot(dm_real_sorted, plot=plt)
    # plt.show()
    if filename:
        plt.savefig(filename)
    plt.close()
    # print("slope: ", res[1][0])
    # print("intercept: ", res[1][1])
    return res[1]


# Paired T-test, 2-tailed
def paired_T_test(x1, x2):
    assert x1.size() == x2.size()
    n = x1.size()

    # Degrees of freedom
    dof = n - 1

    # Calculate means
    m1 = flex.mean(x1)
    m2 = flex.mean(x2)

    # Average difference between observations
    d_avg = sum([(x1[i] - x2[i]) for i in range(n)]) / n

    # Calculate standard error of the average difference
    somm = 0
    for j in range(n):
        _d = x1[j] - x2[j]
        somm += (_d - d_avg) ** 2
    std = math.sqrt(somm / (n - 1))
    sed = std / math.sqrt(n)

    # Find t
    T = (m1 - m2) / sed

    # Considering 5% confidence interval
    alpha = 0.05

    # Find p-value (for 2 tailed test)
    p_val = t.sf(abs(T), dof) * 2

    # If p > alpha => Null hp accepted => equal observations
    res = p_val > alpha
    return T, p_val, res


def compare(data, wdir):
    tab = []
    l = len(data)
    CC_I = np.full((l, l), 0.0)  # Correlation coefficient matrix for intensities
    # CC_s = np.full((l, l), 0.0)     # Correlation coefficient matrix for variances

    for a in range(l):
        CC_I[(a, a)] = 1.0
        tab.append(
            {
                "refl1": a,
                "refl2": a,
                "Correlation Coefficient": 1.0,
                "Normal Probability Plot [slope, intercept]": [0.0, 0.0],
            }
        )

        for b in range(a + 1, l):
            # Match reflections
            m12 = data[a].match(data[b])
            r1 = data[a].select(m12[0])
            r2 = data[b].select(m12[1])
            i1 = r1["intensity.scale.value"]
            i2 = r2["intensity.scale.value"]
            s1 = r1["intensity.scale.variance"]
            s2 = r2["intensity.scale.variance"]
            # Find correlation coefficient
            I = lin_cc(i1, i2)
            CC_I[(a, b)] = CC_I[(b, a)] = I
            # Paired T-test with 95% confidence interval
            T, p_val, T_res = paired_T_test(i1, i2)
            # Normal Probability Plot
            npp = normal_probability_plot(
                i1, i2, s1, s2, wdir / f"compare_refl{a}_and_{b}"
            )
            tab.append(
                {
                    "refl1": a,
                    "refl2": b,
                    "Correlation Coefficient": I,
                    "Normal Probability Plot [slope, intercept]": [npp[0], npp[1]],
                    "Paired T-test [T-statistic, p-value]": [T, p_val],
                    "T-test Null HP accepted": T_res,
                }
            )
    # Plot and save CC
    plt.title("Intensity correlation coefficient")
    plt.imshow(CC_I)
    plt.colorbar()
    plt.savefig(wdir / f"CC_intensities")
    plt.close()

    # Save results to a file
    df = pd.DataFrame(tab)
    with open(wdir / f"Comparison_results.txt", "w") as f:
        f.write("Summary of scaled intensities comparison.\n")
        f.write("Null HP for T-test: No significant change in.\n")
        f.write(tabulate(df, headers="keys", tablefmt="psql", showindex=False))


if __name__ == "__main__":
    args = parser.parse_args()
    data = [flex.reflection_table.from_file(arg) for arg in args.refl_files]
    data = [d.select(d.get_flags(d.flags.scaled)) for d in data]

    if args.output:
        wdir = Path(args.output).expanduser().resolve()
        wdir.mkdir(exist_ok=True)
    else:
        wdir = Path(".").expanduser().resolve()

    compare(data, wdir)
