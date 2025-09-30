"""
This script implements "smart" outlier detection and handling for scaled data
using Robert Von Dreele's method from GSAS-II. This involves calculating a
Z-score for each intensity within its merging group, and for those with high
values of |Z|, the intensity is replaced with the merged intensity for that
merging group. Optionally, the outliers can be rejected instead.

This is not a statistically rigorous method, but has proven practical benefit
in the case of 3D ED data, where dynamical scattering can lead to
sporadically outlying intensity values on some images.
"""

from __future__ import annotations

import logging
import sys
import numpy as np

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

logger = logging.getLogger("dials.smart_outlier")

# Define the master PHIL scope for this program.
phil_scope = libtbx.phil.parse(
    """
    output {
        reflections = smart.refl
            .type = path
            .help = "Option to set filepath for output reflections."
        log = dials.smart_outlier.log
            .type = path
    }
    reject = False
      .type = bool
        .help = "Whether to actually exclude outliers from the output or just"
                "flatten their intensities using the Von Dreele method."
    """
)


def smart_outlier(
    experiments: ExperimentList,
    reflections: flex.reflection_table,
    params: libtbx.phil.scope_extract,
) -> (ExperimentList, flex.reflection_table):
    """
    Write the behaviour of the program as functions and classes outside run().

    Don't include file output here, remember that this function may be re-used
    elsewhere by someone who doesn't need the output written immediately to file.

    It can be especially helpful to document any expected exceptions that might be
    raised, in order to keep track of what needs to be handled in any code that
    re-uses this function.

    Args:
        experiments:  An experiment list.
        reflections:  A reflection table.
        params:       Program parameters, in the form of a scope_extract object,
                      which is the usable form of a parsed PHIL scope.
    Returns:
        The modified reflection table.
    """

    # Filter bad reflections using dials.util.filter methods
    reflections = filter_reflection_table(
        reflections,
        intensity_choice=["scale"],
        d_min=None,
        d_max=None,
        combine_partials=True,
        partiality_threshold=0.4,
    )

    # Scale factor has been applied, so now set to 1.0 so it does not get
    # set again
    reflections["inverse_scale_factor"] = flex.double(reflections.size(), 1.0)

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
    outliers_indices = []
    j = 0
    for i, (hmerge, imerge, sigmerge) in enumerate(
        zip(merged_array.indices(), merged_array.data(), merged_array.sigmas())
    ):

        # Collect the observations of the merging group, and their indices
        obs = []
        sel = []
        while j < len(hkl) and hkl[j] == hmerge:
            obs.append(iobs[j])
            sel.append(j)
            j += 1

        if len(obs) < 2:
            continue

        # Calculate the observed standard deviation of the merging group
        # (not necessarily the same as the merged sigma)
        std = np.std(obs)

        # Loop over the reflections for this merging group, setting outlier
        # intensities to the merged intensity, where an outlier is one with
        # an absolute z-score > 2.0
        for k in sel:
            abs_z = abs((reflections[k]["intensity.scale.value"] - imerge)) / std
            if abs_z > 2.0:
                outliers_indices.append(k)
                reflections["intensity.scale.value"][k] = imerge

    outliers = reflections.select(flex.size_t(outliers_indices))
    logger.info("Total number of outliers identified: %d", len(outliers))

    image_numbers = outliers["xyzcal.px"].parts()[2].iround()
    for iexp, exp in enumerate(experiments):
        images_this_exp = image_numbers.select(outliers["id"] == iexp)
        hist = np.histogram(
            images_this_exp,
            bins=list(range(min(images_this_exp), max(images_this_exp) + 1)),
        )
        header = ["Image", "Number of outliers"]
        rows = zip([str(e) for e in hist[1][:-1]], [str(e) for e in hist[0]])
        table = dials.util.tabulate(rows, header)
        logger.info(f"Experiment {iexp}, outliers per image:\n{table}")

    if params.reject:
        mask = flex.bool(len(reflections), True)
        mask.set_selected(flex.size_t(outliers_indices), False)
        reflections = reflections.select(mask)
        logger.info("Rejected %d reflections", len(outliers))

    return reflections


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
    usage = "dials.smart_outlier [options] scaled.expt scaled.refl"

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
    reflections = smart_outlier(experiments, reflections, params)

    # File output
    logger.info("Writing the reflection table to %s", params.output.reflections)
    reflections.as_file(params.output.reflections)


if __name__ == "__main__":
    run()
