"""
A docstring

This can double as a helpful message which explains how the program is run.
"""

from __future__ import absolute_import, division, print_function

import logging
import sys

import networkx as nx

from annlib_ext import AnnAdaptor
from cctbx import uctbx
import libtbx.phil

from dxtbx.model import ExperimentList
import dials.util.log
from dials.array_family import flex
from dials.report.analysis import scaled_data_as_miller_array
from dials.util.options import OptionParser, flatten_experiments, flatten_reflections

from typing import List


# Define a logger.
logger = logging.getLogger("dials.boilerplate")

# Define the master PHIL scope for this program.
phil_scope = libtbx.phil.parse(
    """
    min_component_size = 50
      .type = int(value_min=0)
    """
)


def do_connected_components(
    experiments,  # type: ExperimentList
    reflection_tables,  # type: flex.reflection_table
    min_component_size=50, # type: int
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
    """

    miller_array = scaled_data_as_miller_array(reflection_tables, experiments, anomalous_flag=False).primitive_setting()
    unique = miller_array.unique_under_symmetry().map_to_asu()
    unique = unique.generate_bijvoet_mates()
    complete_set = unique.complete_set()
    missing_set = complete_set.lone_set(unique)
    missing_set = missing_set.expand_to_p1().customized_copy(crystal_symmetry=missing_set.crystal_symmetry())

    mi = missing_set.indices().as_vec3_double().as_double()
    k = 6
    ann = AnnAdaptor(data=mi, dim=3, k=k)
    ann.query(mi)

    G = nx.Graph()
    distance_cutoff = 2**0.5
    for i in range(missing_set.size()):
        ik = i * k
        for i_ann in range(k):
            if ann.distances[ik + i_ann] <= distance_cutoff:
                j = ann.nn[ik + i_ann]
                G.add_edge(i, j)

    conn = sorted(nx.connected_components(G), key=len, reverse=True)
    unique_mi = []
    unique_ms = []
    for i, c in enumerate(conn):
        ms = missing_set.select(flex.size_t(list(c))).customized_copy(
            crystal_symmetry=miller_array
        ).as_non_anomalous_set().map_to_asu()
        ms = ms.unique_under_symmetry()
        mi = set(ms.indices())
        if mi not in unique_mi:
            unique_ms.append(ms)
            unique_mi.append(mi)

    n_expected = unique.as_non_anomalous_set().complete_set().size()
    unique_ms = sorted(unique_ms, key=lambda ms: ms.size(), reverse=True)
    unique_ms = [ms for ms in unique_ms if ms.size() > min_component_size]
    for ms in unique_ms:
        d_max, d_min = (uctbx.d_star_sq_as_d(ds2) for ds2 in ms.min_max_d_star_sq())
        logger.info("%i reflections (%.1f%%): %.2f-%.2f Å" % (ms.size(), 100 * ms.size()/n_expected, d_max, d_min))
    return unique_ms


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
                            min_component_size=params.min_component_size,
                            )


if __name__ == "__main__":
    with dials.util.show_mail_on_error():
        run()
