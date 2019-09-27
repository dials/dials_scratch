from __future__ import absolute_import, division, print_function

import copy
import random
import sys

import scitbx.random
from dials.algorithms.indexing import assign_indices
from dials.algorithms.indexing.indexer import Indexer
from dials.array_family import flex
from dxtbx.model.detector_helpers import set_slow_fast_beam_centre_mm
from dxtbx.serialize import load
from libtbx.utils import time_log


def run(experiments, reflections, random_seed=42):
    scitbx.random.set_random_seed(random_seed)
    random.seed(random_seed)

    reflections["id"] = flex.int(len(reflections), 0)
    reflections = reflections.select(reflections.get_flags(reflections.flags.indexed))

    beam = experiments[0].beam
    detector = experiments[0].detector
    p_id, (x, y) = detector.get_ray_intersection(beam.get_s0())

    g = scitbx.random.variate(scitbx.random.normal_distribution(mean=0, sigma=2))

    n = 100
    shift_x = g(n)
    shift_y = g(n)

    expected_miller_indices = reflections["miller_index"]
    non_zero_sel = expected_miller_indices != (0, 0, 0)

    misindexed_global = flex.size_t()
    correct_global = flex.size_t()
    misindexed_local = flex.size_t()
    correct_local = flex.size_t()

    global_timer = time_log("global")
    local_timer = time_log("local")

    for d_x, d_y in zip(shift_x, shift_y):
        set_slow_fast_beam_centre_mm(detector, beam, (y + d_y, x + d_x), p_id)

        refl = Indexer.map_centroids_to_reciprocal_space(experiments, reflections)

        refl_global = copy.deepcopy(refl)
        refl_global["id"] = flex.int(len(refl), -1)
        global_timer.start()
        assign_indices.AssignIndicesGlobal()(refl_global, experiments)
        global_timer.stop()

        misindexed_global.append(
            (expected_miller_indices == refl_global["miller_index"])
            .select(non_zero_sel)
            .count(False)
        )
        correct_global.append(
            (expected_miller_indices == refl_global["miller_index"])
            .select(non_zero_sel)
            .count(True)
        )

        refl_local = copy.deepcopy(refl)
        refl_local["id"] = flex.int(len(refl), -1)
        local_timer.start()
        assign_indices.AssignIndicesLocal()(refl_local, experiments)
        local_timer.stop()

        misindexed_local.append(
            (expected_miller_indices == refl_local["miller_index"])
            .select(non_zero_sel)
            .count(False)
        )
        correct_local.append(
            (expected_miller_indices == refl_local["miller_index"])
            .select(non_zero_sel)
            .count(True)
        )

        print("Beam centre shift: (%.2f, %.2f)" % (d_x, d_y))
        print("Misindexed global: %i" % misindexed_global[-1])
        print("Correct global: %i" % correct_global[-1])
        print("Misindexed local: %i" % misindexed_local[-1])
        print("Correct local: %i" % correct_local[-1])
        print()

    print(global_timer.legend)
    print(global_timer.report())
    print(local_timer.report())

    vmax = max(flex.max(correct_global), flex.max(correct_local))

    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    fig, axes = plt.subplots(ncols=2, sharey=True, figsize=(15, 10))
    sc = axes[0].scatter(
        shift_x,
        shift_y,
        vmin=0,
        vmax=1,
        c=correct_global.as_double() / vmax,
        cmap="viridis",
    )
    sc = axes[1].scatter(
        shift_x,
        shift_y,
        vmin=0,
        vmax=1,
        c=correct_local.as_double() / vmax,
        cmap="viridis",
    )
    axes[0].set_title("global")
    axes[1].set_title("local")
    for ax in axes:
        ax.set_aspect("equal")
        ax.set_xlabel("beam centre shift (mm)")
    axes[0].set_ylabel("beam centre shift (mm)")

    cbar = plt.colorbar(sc, ax=axes, shrink=0.5)
    cbar.set_label("Fraction correctly indexed")
    plt.savefig("correctly_indexed.png")


if __name__ == "__main__":
    args = sys.argv[1:]
    assert len(args) == 2
    experiments_json, reflections_file = args
    experiments = load.experiment_list(experiments_json, check_format=False)
    reflections = flex.reflection_table.from_file(reflections_file)
    run(experiments, reflections)
