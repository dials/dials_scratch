# LIBTBX_SET_DISPATCHER_NAME dev.dials.disjoint_split_reflections
"""Split a reflection table into N disjoint sets where observations for each
unique Miller index are split evenly. When N=2 this is similar to the procedure
used in the CC1/2 calculation. By cycling through the set assignments, this
version produces sets that differ in length by only a single reflection"""

from __future__ import absolute_import, division
from __future__ import print_function
import sys
import logging
from itertools import cycle

import cctbx.crystal
import cctbx.miller

import libtbx.phil

import dials.util
import dials.util.log
from dials.util.options import OptionParser, flatten_experiments, flatten_reflections
from dials.array_family import flex

# Define a logger. __name__ cannot be used as this script is called directly.
# Omit the dev prefix to use the DIALS logger
logger = logging.getLogger("dials.disjoint_split_reflections")

# Define the master PHIL scope for this program
phil_scope = libtbx.phil.parse(
    """
output {
    log = dev.dials.disjoint_split_reflections.log
        .type = path

    reflections_prefix = split
        .type = str
        .help = "Filename prefix for the split reflections"
}
N = 2
    .type = int(value_min=2)
"""
)


def split_reflections(experiments, reflections, params):
    """Split reflections into equal sized groups, ensuring that each unique
    Miller index is split proportionately"""

    # Check the space group is the same for all experiments
    sg = experiments[0].crystal.get_space_group()
    if not all((e.crystal.get_space_group() == sg for e in experiments)):
        raise RuntimeError("Not all space groups are equal")

    # create new column containing the reduced Miller index
    xl = experiments[0].crystal
    symm = cctbx.crystal.symmetry(xl.get_unit_cell(), space_group=sg)
    hkl_set = cctbx.miller.set(symm, reflections["miller_index"])
    asu_set = hkl_set.map_to_asu()
    reflections["asu_miller_index"] = asu_set.indices()

    # Set up a cycle through sets and a column to store values
    cyc = cycle(range(params.N))
    reflections["disjoint_set"] = flex.size_t(len(reflections))

    # Loop over unique Miller indices and split matching reflections into sets
    for uniq in set(asu_set.indices()):
        sel = (reflections["miller_index"] == uniq).iselection()
        vals = flex.size_t(next(cyc) for i in sel)
        vals = vals.select(flex.random_permutation(len(vals)))
        reflections["disjoint_set"].set_selected(sel, vals)

    # Now split the table
    splits = []
    for i in range(params.N):
        subset = reflections.select(reflections["disjoint_set"] == i)
        del subset["disjoint_set"]
        splits.append(subset)

    # Write out the splits
    for i, refls in enumerate(splits):
        fname = (
            params.output.reflections_prefix
            + "_{:d}_of_{:d}".format(i + 1, len(splits))
            + ".refl"
        )
        logger.info("Writing {:d} reflections to {}".format(len(refls), fname))
        refls.as_file(fname)


def run(args=None, phil=phil_scope):

    usage = (
        "dev.dials.disjoint_split_reflections [options] integrated.expt integrated.refl"
    )
    parser = OptionParser(
        usage=usage,
        phil=phil,
        read_reflections=True,
        read_experiments=True,
        check_format=False,
        epilog=__doc__,
    )

    params, options = parser.parse_args(args=args, show_diff_phil=False)

    # Configure the logging
    dials.util.log.config(verbosity=options.verbose, logfile=params.output.log)

    # Log the PHIL diff. Done _after_ log config so written to log
    diff_phil = parser.diff_phil.as_str()
    if diff_phil:
        logger.info("The following parameters have been modified:\n%s", diff_phil)

    experiments = flatten_experiments(params.input.experiments)
    reflections = flatten_reflections(params.input.reflections)

    # Do any simple loading validation here
    if len(reflections) != 1:
        sys.exit("Error: Exactly one reflection file needed")

    reflections = reflections[0]
    logger.info("Loaded {0} reflections".format(len(reflections)))
    split_reflections(experiments, reflections, params)


if __name__ == "__main__":
    with dials.util.show_mail_on_error():
        run()
