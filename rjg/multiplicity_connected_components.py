"""
A docstring

This can double as a helpful message which explains how the program is run.
"""

from __future__ import absolute_import, division, print_function

import logging
import sys

from cctbx import uctbx
import libtbx.phil

from dials.array_family import flex

import dials.util

import dials.util.log

from dials.util.options import OptionParser, flatten_experiments, flatten_reflections

from dxtbx.model import ExperimentList

from typing import List


# Define a logger.
logger = logging.getLogger("dials.boilerplate")

# Define the master PHIL scope for this program.
phil_scope = libtbx.phil.parse(
    """
    anomalous = True
        .type = bool
    expand_to_p1 = True
        .type = bool
    """
)


def do_connected_components(
    experiments,  # type: ExperimentList
    reflection_tables,  # type: flex.reflection_table
    expand_to_p1=True,  # type: bool
    anomalous=True, # type: bool
):  # type: (...) -> [{}]
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
        params:       Some parameters, in the form of a scope_extract object,
                      which is the usable form of a parsed PHIL scope.
    """

    from dials.report.analysis import scaled_data_as_miller_array

    miller_array = scaled_data_as_miller_array(reflection_tables, experiments).primitive_setting()
    unique = miller_array.unique_under_symmetry().map_to_asu()
    missing_set = unique.complete_set().lone_set(unique)
    if anomalous:
        missing_set = missing_set.as_anomalous_set()
    else:
        missing_set = missing_set.as_non_anomalous_set()
    if expand_to_p1:
        missing_set = missing_set.expand_to_p1()

    from annlib_ext import AnnAdaptor
    mi = missing_set.indices().as_vec3_double().as_double()
    k = 6
    ann = AnnAdaptor(data=mi, dim=3, k=k)
    ann.query(mi)

    import networkx as nx
    G = nx.Graph()
    distance_cutoff = 2**0.5
    for i in range(missing_set.size()):
        ik = i * k
        for i_ann in range(k):
            if ann.distances[ik + i_ann] <= distance_cutoff:
                j = ann.nn[ik + i_ann]
                G.add_edge(i, j)

    conn = sorted(nx.connected_components(G), key=len, reverse=True)
    conn = [c for c in conn if len(c) > 25]
    for c in conn:
        c_ms = missing_set.select(flex.size_t(list(c)))
        d_max, d_min = (uctbx.d_star_sq_as_d(ds2) for ds2 in c_ms.min_max_d_star_sq())
        logger.info("%s reflections: %.2f-%.2f Å" % (len(c), d_max, d_min))
    return conn

def run(args=None, phil=phil_scope):  # type: (List[str], libtbx.phil.scope) -> None
    """
    Check command-line input and call other functions to do the legwork.

    Run the script, parsing arguments found in 'args' and using the PHIL scope
    defined in 'phil'.

    Args:
        args: The arguments supplied by the user (default: sys.argv[1:])
        phil: The PHIL scope definition (default: phil_scope, the master PHIL scope
        for this program).
    """
    usage = "dials.command_name [options] imported.expt strong.refl"

    parser = OptionParser(
        usage=usage,
        phil=phil,
        read_reflections=True,
        read_experiments=True,
        check_format=False,
        epilog=__doc__,
    )

    params, options = parser.parse_args(args=args, show_diff_phil=False)

    # Configure the logging.
    dials.util.log.config(options.verbose)

    # Log the difference between the PHIL scope definition and the active PHIL scope,
    # which will include the parsed user inputs.
    diff_phil = parser.diff_phil.as_str()
    if diff_phil:
        logger.info("The following parameters have been modified:\n%s", diff_phil)

    experiments = flatten_experiments(params.input.experiments)
    reflections = flatten_reflections(params.input.reflections)

    if not experiments or not reflections:
        parser.print_help()
        return
    if len(reflections) != 1 and len(experiments) != len(reflections):
        sys.exit("Number of experiments must equal the number of reflection tables")

    from dials.util.multi_dataset_handling import (
        assign_unique_identifiers,
        parse_multiple_datasets,
    )
    reflections = parse_multiple_datasets(reflections)
    experiments, reflections = assign_unique_identifiers(experiments, reflections)

    # Do whatever this program is supposed to do.
    do_connected_components(experiments, reflections,
                            anomalous=params.anomalous,
                            expand_to_p1=params.expand_to_p1)


if __name__ == "__main__":
    with dials.util.show_mail_on_error():
        run()
