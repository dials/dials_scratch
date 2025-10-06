"""
This script identifies overlapping reflections in scaled data from two lattices
and rejects those that are closer than a user-defined minimum distance.
"""

from __future__ import annotations

import logging
import sys
import matplotlib.pyplot as plt


from scitbx.math import five_number_summary
from annlib_ext import AnnAdaptor
from cctbx import miller
import libtbx.phil
from dxtbx.model import ExperimentList
import dials.util
import dials.util.log
from dials.array_family import flex
from dials.util.options import ArgumentParser, reflections_and_experiments_from_files
from dials.util.version import dials_version
from dials.algorithms.scaling.scaling_library import (
    scaled_data_as_miller_array,
)
import dials.util

from dials.util.filter_reflections import filter_reflection_table

logger = logging.getLogger("dials.reject_overlaps")

# Define the master PHIL scope for this program.
phil_scope = libtbx.phil.parse(
    """
    output {
        reflections = filtered.refl
            .type = path
            .help = "Option to set filepath for output reflections."
        overlaps = None
            .type = path
            .help = "Option to set filepath for output overlap reflections."
                    "This can be useful to view the overlaps on the image."
        log = dials.reject_overlaps.log
            .type = path
        plot = None
            .type = path
            .help = "Option to set filepath for output plot of I/Ih vs overlap distance."
    }
    minimum_distance = 2.0
        .type = float
        .help = "Minimum distance in pixels below which reflections are rejected."
    plot_bin_width = 0.2
        .type = float
        .help = "Bin width in pixels for plot of I/Ih vs overlap distance."
    """
)


def reject_overlaps(
    experiments: ExperimentList,
    reflections: flex.reflection_table,
    params: libtbx.phil.scope_extract,
) -> (ExperimentList, flex.reflection_table):
    """
    Identify potentially overlapping reflections in multi-lattice data and reject
    those that are closer than a user-defined minimum distance.

    Args:
        experiments:  An experiment list.
        reflections:  A reflection table.
        params:       Program parameters, in the form of a scope_extract object,
                      which is the usable form of a parsed PHIL scope.
    Returns:
        A tuple containing the modified reflection table and the overlaps table.
    """

    # Keep track of some columns that will be deleted by the filtering step
    reflections["original_index"] = flex.size_t(range(len(reflections)))
    intensity_sum_value = reflections["intensity.sum.value"].deep_copy()
    intensity_sum_variance = reflections["intensity.sum.variance"].deep_copy()
    intensity_prf_value = reflections["intensity.prf.value"].deep_copy()
    intensity_prf_variance = reflections["intensity.prf.variance"].deep_copy()

    # Filter bad reflections using dials.util.filter methods
    reflections = filter_reflection_table(
        reflections,
        intensity_choice=["scale"],
        d_min=None,
        d_max=None,
        combine_partials=True,
        partiality_threshold=0.4,
    )

    # Put back the intensity columns that were removed by filtering
    reflections["intensity.sum.value"] = intensity_sum_value.select(
        reflections["original_index"]
    )
    reflections["intensity.sum.variance"] = intensity_sum_variance.select(
        reflections["original_index"]
    )
    reflections["intensity.prf.value"] = intensity_prf_value.select(
        reflections["original_index"]
    )
    reflections["intensity.prf.variance"] = intensity_prf_variance.select(
        reflections["original_index"]
    )
    del reflections["original_index"]

    # Scale factor has been applied, so now set to 1.0 so it does not get
    # set again. Create a new column to keep the recalculated scale factor
    # (which will be different as it is simple ratio of intensity to merged
    # intensity)
    reflections["inverse_scale_factor"] = flex.double(reflections.size(), 1.0)
    reflections["recalc_scale_factor"] = flex.double(reflections.size(), 1.0)

    # Now extract the miller array
    scaled_array = scaled_data_as_miller_array(
        [reflections], experiments, anomalous_flag=False
    )

    # Sort data into merging groups
    asu_set = miller.set.map_to_asu(scaled_array)
    perm = asu_set.sort_permutation(by_value="packed_indices")
    reflections = reflections.select(perm)
    hkl = asu_set.indices().select(perm)
    iobs = scaled_array.data().select(perm)
    sigmas = scaled_array.sigmas().select(perm)

    # Get the merged data
    merged = scaled_array.merge_equivalents(use_internal_variance=False)
    merged_array = merged.array()

    # Outer loop over the merged reflections.
    j = 0
    for i, (hmerge, imerge) in enumerate(
        zip(merged_array.indices(), merged_array.data())
    ):

        # Collect the observations of the merging group, and their indices
        obs = []
        sel = []
        while j < len(hkl) and hkl[j] == hmerge:
            obs.append(iobs[j])
            sel.append(j)
            j += 1

        # Loop over the reflections for this merging group, calculating the
        # scale factor between this reflection and the merged intensity
        for k in sel:
            reflections["recalc_scale_factor"][k] = (
                reflections["intensity.scale.value"][k] / imerge
            )

    # Alternate algorithm - gives essentially the same results once the filtering
    # has been applied
    # reflections["intensity"] = reflections["intensity.scale.value"]
    # reflections["variance"] = reflections["intensity.scale.variance"]
    # from dials.algorithms.scaling.Ih_table import IhTable
    # space_group = experiments[0].crystal.get_space_group()
    # ih = IhTable([reflections], space_group, nblocks=1)
    # ih = ih.blocked_data_list[0].as_reflection_table()
    # ih.sort("loc_indices")
    # reflections["recalc_scale_factor2"] = reflections["intensity.scale.value"] / ih["Ih_values"]

    # Separate reflections for each lattice
    r1 = reflections.select(reflections["id"] == 0)
    r2 = reflections.select(reflections["id"] == 1)

    # Create the KD Tree and for every reflection in r1 find the distance to the
    # nearest neighbour in r2. The indices of the nearest neighbours are then in
    # ann.nn, but we won't use them here.
    ann = AnnAdaptor(r2["xyzcal.px"].as_double(), dim=3, k=1)
    ann.query(r1["xyzcal.px"].as_double())
    r1["distances"] = flex.sqrt(ann.distances)

    # Now repeat for the second lattice
    ann = AnnAdaptor(r1["xyzcal.px"].as_double(), dim=3, k=1)
    ann.query(r2["xyzcal.px"].as_double())
    r2["distances"] = flex.sqrt(ann.distances)

    # Rejoin the reflection tables
    reflections = flex.reflection_table.concat((r1, r2))

    overlaps = reflections.select(reflections["distances"] < params.minimum_distance)
    logger.info(
        f"Found {len(overlaps)} overlaps within {params.minimum_distance} pixels"
    )

    # Create plot of I/Ih vs overlap distance
    if params.output.plot is not None:
        bins = reflections["distances"] / params.plot_bin_width
        bins = flex.floor(bins).iround()
        medians = []
        distances = []
        for i in range(max(bins) + 1):
            sel = bins == i
            # Only consider bins with at least 100 reflections
            if sel.count(True) < 100:
                continue
            distances.append((i + 1) * params.plot_bin_width)
            medians.append(
                five_number_summary(reflections["recalc_scale_factor"].select(sel))[2]
            )

        fig, ax = plt.subplots()
        ax.scatter(distances, medians)
        ax.axhline(y=1.0, color="r", linestyle="-")
        ax.set_ylabel("I/Ih")
        ax.set_xlabel("Overlap distance (pixels)")
        plt.savefig(params.output.plot)

    # Reject the overlaps
    sel = reflections["distances"] >= params.minimum_distance
    reflections = reflections.select(sel)

    return reflections, overlaps


@dials.util.show_mail_on_error()
def run(args: list[str] = None, phil: libtbx.phil.scope = phil_scope) -> None:
    """
    Check command-line input and call other functions to do the legwork.

    Run the script, parsing arguments found in 'args' and using the PHIL scope
    defined in 'phil'.

    Args:
        args: The arguments supplied by the user (default: sys.argv[1:])
        phil: The PHIL scope definition (default: phil_scope, the master PHIL scope
        for this program).
    """
    usage = "dials.reject_overlaps [options] scaled.expt scaled.refl"

    parser = ArgumentParser(
        usage=usage,
        phil=phil,
        read_reflections=True,
        read_experiments=True,
        check_format=False,
        epilog=__doc__,
    )

    params, options = parser.parse_args(args=args, show_diff_phil=False)

    # Configure the logging.
    dials.util.log.config(options.verbose, logfile=params.output.log)

    # Log the dials version
    logger.info(dials_version())

    # Log the difference between the PHIL scope definition and the active PHIL scope,
    # which will include the parsed user inputs.
    diff_phil = parser.diff_phil.as_str()
    if diff_phil:
        logger.info("The following parameters have been modified:\n%s", diff_phil)

    reflections, experiments = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )

    # Try to load the models and data
    nexp = len(experiments)
    if nexp == 0 or len(reflections) == 0:
        parser.print_help()
        return
    if len(reflections) > 1:
        sys.exit("Only one reflections list can be imported at present")
    reflections = reflections[0]

    # Run the algorithm
    reflections, overlaps = reject_overlaps(experiments, reflections, params)

    # File output
    logger.info("Writing the reflection table to %s", params.output.reflections)
    reflections.as_file(params.output.reflections)

    if params.output.overlaps is not None:
        logger.info(
            "Writing the overlaps reflection table to %s", params.output.overlaps
        )
        overlaps.as_file(params.output.overlaps)


if __name__ == "__main__":
    run()
