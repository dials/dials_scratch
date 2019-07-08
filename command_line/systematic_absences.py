#!/usr/bin/env python
# coding: utf-8
"""
Command line script to assess the space group symmetry (for MX datasets)

Currently only tests axial absences to determine which of the 65 MX space
groups is most likely.
"""
from __future__ import absolute_import, division, print_function
import logging
import sys
from dials.util import log, show_mail_on_error
from dials.util.options import OptionParser, flatten_reflections, flatten_experiments
from dials.util.version import dials_version
from dials_scratch.jbe.sys_abs.laue_groups_info import laue_groups, score_screw_axes, score_space_groups
from dials_scratch.jbe.sys_abs.screw_axes import ScrewAxisObserver
from dials.util.filter_reflections import ScaleIntensityReducer

from libtbx import phil

help_message = """Program to assess space group symmetry (for MX datasets)."""

logger = logging.getLogger("dials")
info_handle = log.info_handle(logger)
phil_scope = phil.parse("""
""")


def run_sys_abs_checks(params, experiments, reflections):

    reflections = ScaleIntensityReducer.reduce_on_intensities(reflections[0])
    print("Number of reflections in dataset: %s" % reflections.size())
    reflections['intensity'] = reflections['intensity.scale.value']
    reflections['variance'] = reflections['intensity.scale.variance']

    # Get the laue class from the space group
    laue_group = str(experiments[0].crystal.get_space_group().build_derived_patterson_group().info())
    print("Laue group: %s" % laue_group)
    if not laue_group in laue_groups:
        print("No absences to check for this laue group")
        return
    screw_axis_scores = score_screw_axes(laue_groups[laue_group], reflections)
    best_sg = score_space_groups(screw_axis_scores, laue_groups[laue_group])

    # now set space group

    observer = ScrewAxisObserver()
    #print(observer.data)


    # also check for lattice centering?

    # option to choose setting convention (.g C2 or I2  maybe standard, reference?)

    # also check for obverse/reverse twinning in rhombohedral?

    # give some warnings about potential NCS? etc
    pass

def plots(params, experiments, reflections):
    pass


def run(args=None):
    """Run the script from the command-line."""
    usage = """Usage: dials.systematic_absences scaled.refl scaled.expt [options]"""

    parser = OptionParser(
        usage=usage,
        read_experiments=True,
        read_reflections=True,
        phil=phil_scope,
        check_format=False,
        epilog=help_message,
    )
    params, _ = parser.parse_args(args=args, show_diff_phil=False)

    if not params.input.experiments or not params.input.reflections:
        parser.print_help()
        sys.exit()

    reflections = flatten_reflections(params.input.reflections)
    experiments = flatten_experiments(params.input.experiments)

    #log.config(verbosity=1, info=params.output.log, debug=params.output.debug.log)
    logger.info(dials_version())

    diff_phil = parser.diff_phil.as_str()
    if diff_phil != "":
        logger.info("The following parameters have been modified:\n")
        logger.info(diff_phil)

    ### Assert that all data have been scaled with dials - should only be
    # able to input one reflection table and experimentlist that are
    # matching and scaled together.

    run_sys_abs_checks(params, experiments, reflections)


if __name__ == "__main__":
    with show_mail_on_error():
        run()
