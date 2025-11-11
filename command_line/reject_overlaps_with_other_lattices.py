"""
This script identifies overlapping reflections in an integrated dataset with
any integrated reflection that belongs to a different lattice.
"""

from __future__ import annotations

import logging
import sys

from annlib_ext import AnnAdaptor
import libtbx.phil
from dxtbx.model import ExperimentList
import dials.util
import dials.util.log
from dials.array_family import flex
from dials.util.options import ArgumentParser, reflections_and_experiments_from_files
from dials.util.version import dials_version
import dials.util


logger = logging.getLogger("dials.reject_overlaps_with_other_lattices")

# Define the master PHIL scope for this program.
phil_scope = libtbx.phil.parse(
    """
    output {
        experiments = filtered.expt
            .type = path
            .help = "Option to set filepath for output experiments."
        reflections = filtered.refl
            .type = path
            .help = "Option to set filepath for output reflections."
        overlaps = None
            .type = path
            .help = "Option to set filepath for output overlap reflections."
                    "This can be useful to view the overlaps on the image."
        log = dials.reject_overlaps.log
            .type = path
    }
    reference_id = 0
        .type = int(value_min=0)
        .help = "Experiment id of the chosen lattice from which reflections will be rejected"
    minimum_distance = 2.0
        .type = float
        .help = "Minimum distance in pixels below which reflections are rejected."
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

    selected_experiments = ExperimentList()
    selected_experiments.append(experiments[params.reference_id])

    r1 = reflections.select(reflections["id"] == params.reference_id)
    r2 = reflections.select(reflections["id"] != params.reference_id)

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

    overlaps = r1.select(r1["distances"] < params.minimum_distance)
    logger.info(
        f"Found {len(overlaps)} overlaps within {params.minimum_distance} pixels"
    )

    # Reject the overlaps
    sel = r1["distances"] >= params.minimum_distance
    r1 = r1.select(sel)

    return r1, overlaps, selected_experiments


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
    reflections, overlaps, experiments = reject_overlaps(
        experiments, reflections, params
    )

    # File output
    logger.info("Writing the reflection table to %s", params.output.reflections)
    reflections.as_file(params.output.reflections)

    logger.info("Writing the experiments to %s", params.output.experiments)
    experiments.as_file(params.output.experiments)

    if params.output.overlaps is not None:
        logger.info(
            "Writing the overlaps reflection table to %s", params.output.overlaps
        )
        overlaps.as_file(params.output.overlaps)


if __name__ == "__main__":
    run()
