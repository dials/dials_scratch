import sys
import time
import math
import argparse
import logging

import numpy as np
import pandas as pd

from pathlib import Path
from tabulate import tabulate
from scipy.stats import t, probplot

from matplotlib import pyplot as plt

from dials.array_family import flex

logger = logging.getLogger("CompareScaledIntensities")
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter("%(levelname)s    ||  %(message)s")
CH = logging.StreamHandler(sys.stdout)
CH.setLevel(logging.DEBUG)
CH.setFormatter(formatter)
logger.addHandler(CH)

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
    """Return the linear correlation coefficient."""
    lc = flex.linear_correlation(x1, x2)
    return lc.coefficient()


# Normal probability plot
def normal_probability_plot(x1, x2, v1, v2, filename):
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    k = sum(x1 * x2) / sum(x2 ** 2)
    # n = x1.size()
    dm_real = (x1 - k * x2) / flex.sqrt(v1 + k ** 2 * v2)
    dm_idx = np.argsort(dm_real.as_numpy_array())
    dm_real_sorted = dm_real.as_numpy_array()[dm_idx]
    # dm_real_sorted = np.sort(dm_real.as_numpy_array())

    res = probplot(dm_real_sorted, plot=ax)
    plt.text(0, 1, f"slope = {res[1][0]:.3f}", transform=ax.transAxes)

    if filename:
        plt.savefig(filename)
    plt.close()
    return res, dm_idx


# Paired T-test, 2-tailed, default confidence interval 5%
def paired_T_test(x1, x2, alpha=0.05):
    """Return the p-value for a two tailed paired T-test."""
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

    # Find p-value (for 2 tailed test)
    p_val = t.sf(abs(T), dof) * 2

    # If p > alpha => Null hp accepted => equal observations
    # res = p_val > alpha
    return p_val


def matcher(d0, d1):
    """Return the reflections which are common to d0, d1 in a common order"""

    d0 = d0.select(d0["partiality"] > 0.99)
    d1 = d1.select(d1["partiality"] > 0.99)

    id0 = d0["id"].as_numpy_array()

    hkl = d0["miller_index"].as_vec3_double().parts()
    h0 = hkl[0].iround().as_numpy_array()
    k0 = hkl[1].iround().as_numpy_array()
    l0 = hkl[2].iround().as_numpy_array()
    e0 = d0["entering"].as_numpy_array()
    n0 = np.array(range(len(e0)))

    x0 = {"h": h0, "k": k0, "l": l0, "e": e0, "id": id0, "n0": n0}

    p0 = pd.DataFrame(data=x0, columns=["h", "k", "l", "e", "id", "n0"])

    id1 = d1["id"].as_numpy_array()

    hkl = d1["miller_index"].as_vec3_double().parts()
    h1 = hkl[0].iround().as_numpy_array()
    k1 = hkl[1].iround().as_numpy_array()
    l1 = hkl[2].iround().as_numpy_array()
    e1 = d1["entering"].as_numpy_array()
    n1 = np.array(range(len(e1)))

    x1 = {"h": h1, "k": k1, "l": l1, "e": e1, "id": id1, "n1": n1}

    p1 = pd.DataFrame(data=x1, columns=["h", "k", "l", "e", "id", "n1"])

    merged = pd.merge(p0, p1, on=["h", "k", "l", "e", "id"], how="inner")

    n0 = merged["n0"].to_numpy()
    n1 = merged["n1"].to_numpy()

    # print(d0.size(), d1.size())

    d0 = d0.select(flex.size_t(n0))
    d1 = d1.select(flex.size_t(n1))

    return d0, d1


def extract_hkl_and_xyz(r1, r2, idx):
    hkl = []
    xyz1 = []
    xyz2 = []
    miller = r1[
        "miller_index"
    ].as_vec3_double()  # matching done by hkl so just one should be enough
    pos1 = r1["xyzobs.px.value"]
    pos2 = r2["xyzobs.px.value"]
    for i in idx:
        hkl.append(miller[i])
        xyz1.append(pos1[i])
        xyz2.append(pos2[i])
    return hkl, xyz1, xyz2


def compare(data, wdir):
    tab = []
    l = len(data)
    logger.info(f"Number of reflection files: {l}.")
    CC_I = np.full((l, l), 0.0)  # Correlation coefficient matrix for intensities
    # CC_s = np.full((l, l), 0.0)     # Correlation coefficient matrix for variances
    DF = pd.DataFrame()

    for a in range(l):
        CC_I[(a, a)] = 1.0
        tab.append(
            {
                "refl1": a,
                "refl2": a,
                "Correlation Coefficient": 1.0,
                "Normal Probability Plot [slope, intercept]": "-",
                "Paired T-test [p-value]": "-",
                "Number of matched reflections": "-",
            }
        )

        for b in range(a + 1, l):
            logger.info(f"Match reflections {a} with {b}.")
            # Match reflections
            r1, r2 = matcher(data[a], data[b])
            i1 = r1["intensity.scale.value"]
            i2 = r2["intensity.scale.value"]
            s1 = r1["intensity.scale.variance"]
            s2 = r2["intensity.scale.variance"]
            logger.info(f"Number of matched reflections: {i1.size()}")
            # Find correlation coefficient
            I = lin_cc(i1, i2)
            logger.info(f"Correlation coefficient: {I}")
            CC_I[(a, b)] = CC_I[(b, a)] = I
            # Paired T-test
            p_val = paired_T_test(i1, i2)
            logger.info(
                f"P-value for paired t-test with 5% confidence interval: {p_val}"
            )
            # Normal Probability Plot
            npp, idx = normal_probability_plot(
                i1, i2, s1, s2, wdir / f"compare_refl{a}_and_{b}"
            )
            logger.info(
                f"Normal Probability Plot \n slope: {npp[1][0]} \n intercept: {npp[1][1]} \n"
            )
            tab.append(
                {
                    "refl1": a,
                    "refl2": b,
                    "Correlation Coefficient": I,
                    "Normal Probability Plot [slope, intercept]": [
                        npp[1][0],
                        npp[1][1],
                    ],
                    "Paired T-test [p-value]": p_val,
                    "Number of matched reflections": i1.size(),
                }
            )

            # Get hkl and positions
            HKL, XYZ1, XYZ2 = extract_hkl_and_xyz(r1, r2, idx)
            # Create a dataframe with the data
            x = {
                "reflections": f"{a} - {b}",
                "I1": i1.as_numpy_array()[idx],
                "s1": s1.as_numpy_array()[idx],
                "I2": i2.as_numpy_array()[idx],
                "s2": s2.as_numpy_array()[idx],
                "index": idx,
                "hkl": HKL,
                "xyzobs1": XYZ1,
                "xyzobs2": XYZ2,
                "observations": npp[0][1],
                "quantiles": npp[0][0],
            }
            df = pd.DataFrame(data=x)
            DF = pd.concat([DF, df])

    # Plot and save CC
    plt.title("Intensity correlation coefficient")
    plt.imshow(CC_I)
    plt.colorbar()
    plt.savefig(wdir / f"CC_intensities")
    plt.close()

    # Save results to a file
    res = pd.DataFrame(tab)
    logger.info("Save results table to a txt file.")
    with open(wdir / f"Comparison_results.txt", "w") as f:
        f.write("Summary of scaled intensities comparison.\n")
        f.write(tabulate(res, headers="keys", tablefmt="psql", showindex=False))

    # Save dataframe to a csv file for further processing (can be opened with pd.read_csv(filename, index_col=0))
    logger.info("Save Reflections dataframe to a csv file.")
    DF.to_csv(wdir / "Reflections.csv")


if __name__ == "__main__":
    tic = time.process_time()
    logger.info("Look through reflection files and compare scaled intensities.")
    args = parser.parse_args()
    data = [flex.reflection_table.from_file(arg) for arg in args.refl_files]
    data = [d.select(d.get_flags(d.flags.scaled)) for d in data]

    if args.output:
        wdir = Path(args.output).expanduser().resolve()
        wdir.mkdir(exist_ok=True)
    else:
        wdir = Path(".").expanduser().resolve()
    logger.info(f"Results will be saved in {wdir}.")

    compare(data, wdir)
    toc = time.process_time()
    # print(f"Time taken: {toc-tic:.4f} s.")
    logger.info(f"Time taken: {toc-tic:.4f} s.")
