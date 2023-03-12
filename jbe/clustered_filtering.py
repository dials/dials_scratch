# load some scaled data

from dials.array_family import flex
import numpy as np
from dxtbx import flumpy
import logging
from dials.util import log
from math import pi
from scipy.optimize import minimize
from dials.algorithms.scaling.Ih_table import IhTable
from dials.util.version import dials_version

from dials.util.options import OptionParser, reflections_and_experiments_from_files
from random import uniform, seed

seed(0)

logger = logging.getLogger("dials.cluster_filter")


def gaussian_prob_vec(xs, vars, mu, sigma):
    comb_sigmasq = sigma**2 + vars
    u = (xs - mu) ** 2 / (2 * comb_sigmasq)
    return np.exp(-u) / np.sqrt(2 * pi * comb_sigmasq)


def calc_r_vectors(cluster_weights, mus, sigmas, intensities, variances):
    probs = np.empty(shape=(2, intensities.size))
    for idx, (w, mu, sig) in enumerate(zip(cluster_weights, mus, sigmas)):
        probs[idx] = w * gaussian_prob_vec(intensities, variances, mu, sig)
    sumprob = np.sum(probs, axis=0)
    probs /= sumprob
    return probs


def calc_ll(params, intensities, variances):

    cluster_weights = (params[0], 1.0 - params[0])
    mus = params[1:3]
    sigmas = params[3:5]

    probs = np.empty(shape=(2, intensities.size))
    for idx, (w, mu, sig) in enumerate(zip(cluster_weights, mus, sigmas)):
        probs[idx] = w * gaussian_prob_vec(intensities, variances, mu, sig)
    sumprob = np.sum(probs, axis=0)
    ll = -1.0 * np.sum(np.log(sumprob))
    return ll


def test_group(intensities, variances, miller_index):
    if len(intensities) < 10:
        return
    logger.info(f"Assessing {miller_index}")

    res = minimize(
        calc_ll,
        np.array([0.5, 0.0, 400.0, 20.0, 200.0]),
        args=(intensities, variances),
        method="L-BFGS-B",
        bounds=[
            (0.05, 0.95),
            (-40.0, 40.0),
            (np.min(intensities), np.max(intensities)),
            (1, 100),
            (1, 1000),
        ],
    )
    if not res.success:
        return

    cluster_weights = (res.x[0], 1.0 - res.x[0])
    mus = res.x[1:3]
    sigmas = res.x[3:5]

    logger.info(f"Cluster weights: {cluster_weights}")
    logger.info(f"Gaussian means: {mus}")
    logger.info(f"Gaussian sigmas: {sigmas}")
    # if abs(mus[0]) > sigmas[0]:
    #    logger.info("First gaussian not centred about zero")
    #    #return
    # else:
    #     logger.info("First gaussian centred about zero")
    if mus[1] > sigmas[1]:
        logger.info("Second gaussian seems like a significant signal")
    # else:
    #    logger.info("Second gaussian doesn't seems like a significant signal")
    #    return
    if cluster_weights[0] > 0.89:
        logger.info("Majority of signal in zero gaussian, unsuccessful separation")
        return
    # OK so the fit was okay, but was it better than a single gaussian fit?

    logger.info(f"Separation appears successful, taking statistical sample")
    # now take a statistical sample of the "real" gaussian
    final_rs = calc_r_vectors(cluster_weights, mus, sigmas, intensities, variances)

    in_real = flex.bool()
    real_xs = []
    for r, i in zip(final_rs.T, intensities):
        if r[0] < uniform(0, 1):
            in_real.append(True)
            real_xs.append(i)
        else:
            in_real.append(False)
    logger.info(f"Removed {len(in_real) - len(real_xs)}/{len(in_real)} reflections")
    return in_real, real_xs, res


def run():

    parser = OptionParser(
        read_experiments=True,
        read_reflections=True,
        check_format=False,
        epilog=__doc__,
    )
    params, _, args = parser.parse_args(show_diff_phil=False, return_unhandled=True)
    log.config(verbosity=1, logfile="dials.cluster_filter.log")
    logger.info(dials_version())

    diff_phil = parser.diff_phil.as_str()
    if diff_phil:
        logger.info("The following parameters have been modified:\n%s", diff_phil)

    reflections, expts = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )
    refls = reflections[0]

    refls["intensity"] = refls["intensity.scale.value"]
    refls["variance"] = refls["intensity.scale.variance"]
    refls["initial_index"] = flex.size_t_range(refls.size())
    good_refls = refls.select(refls.get_flags(refls.flags.scaled))

    Ih_table = IhTable(
        [good_refls],
        space_group=expts[0].crystal.get_space_group(),
        indices_lists=[good_refls["initial_index"]],
    )
    block = Ih_table.blocked_data_list[0]

    to_exclude = flex.size_t([])
    for group_idx in range(0, block.n_groups):

        sel = flex.bool(block.n_groups, False)
        sel[group_idx] = True
        sel_block = block.select_on_groups(sel)

        sel = sel_block.intensities / (sel_block.variances**0.5) > -1.0
        sel_block = sel_block.select(sel)
        if sel_block.size:
            I = sel_block.intensities / sel_block.inverse_scale_factors
            V = sel_block.variances / (sel_block.inverse_scale_factors**2)
            result = test_group(I, V, sel_block.asu_miller_index[0])
            if result:
                in_real = result[0]
                to_exclude.extend(
                    flumpy.from_numpy(
                        sel_block.Ih_table["loc_indices"].to_numpy()
                    ).select(~in_real)
                )

    logger.info(to_exclude.size())
    logger.info(refls.size())

    bad = flex.bool(refls.size(), False)
    bad.set_selected(to_exclude, True)
    refls = refls.select(~bad)
    refls.as_file("filtered.refl")
    logger.info("Done")


if __name__ == "__main__":
    run()
