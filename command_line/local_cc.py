#!/usr/bin/env python
#
# dials.local_cc.py
#
#  Copyright (C) 2019 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

# LIBTBX_SET_DISPATCHER_NAME dev.dials.local_cc

from __future__ import absolute_import, division
import libtbx.load_env
from libtbx.utils import Sorry
from cctbx import miller
from libtbx.phil import parse
from iotbx.reflection_file_reader import any_reflection_file
from dials.util import log
from dials.util.options import OptionParser
from dials.util.options import flatten_reflections, flatten_experiments
from dials.array_family import flex
from dials_scratch.jmp.local_statistics.local_cc import plot_local_cc_half
from dials_scratch.jmp.local_statistics.local_cc import plot_local_cc_vs_ref
from scitbx import matrix
import logging


logger = logging.getLogger("dials")

help_message = """

This script computes local correlations in intensities

"""

# Create the phil scope
phil_scope = parse(
    """

    input {

        hklin = None
            .type = str
            .help = "The input MTZ file"
        
        hklref = None
            .type = str
            .help = "The reference MTZ file"

    }

    output {

        record = False
            .type = bool
            .help = "Record an animation"

        directory = "output"
            .type = str
            .help = "The output directory"
    }

    intensity = *sum prf
        .type = choice
        .help = "In the case of integrated data which intensity to use"

    kernel_size = 2
        .type = int
        .help = "The kernel size to compute the local CC 1/2 (2*kernel_size+1)"

    point_size=0.0025
        .type = float
        .help = "The size of the plotted points"

    """,
    process_includes=True,
)


class Script(object):
    def __init__(self):

        # The script usage
        usage = "usage: %s [options] experiment.json" % libtbx.env.dispatcher_name

        # Create the parser
        self.parser = OptionParser(
            usage=usage,
            phil=phil_scope,
            epilog=help_message,
            read_experiments=True,
            read_reflections=True,
            check_format=False,
        )

    def run(self):
        from time import time

        # Check the number of arguments is correct
        start_time = time()

        # Parse the command line
        params, options = self.parser.parse_args(show_diff_phil=False)
        reflections = flatten_reflections(params.input.reflections)
        experiments = flatten_experiments(params.input.experiments)
        if params.input.hklin is None:
            if len(reflections) == 0 or len(experiments) == 0:
                self.parser.print_help()
                return
            elif len(reflections) != 1:
                raise Sorry("more than 1 reflection file was given")
            elif len(experiments) == 0:
                raise Sorry("no experiment list was specified")
            reflections = reflections[0]
            reflections, A = self.process_reflections_and_experiments(
                reflections, experiments, params.intensity
            )
            if params.input.hklref is not None:
                raise Sorry("hklin can only be set along side hklin")
        else:
            if params.input.hklref is None:
                reflections, A = self.process_mtzfile(params.input.hklin)
            else:
                reflections, A = self.process_mtzfile(
                    params.input.hklin, merge_equivalents=True
                )

        # Read the reference data
        if params.input.hklref is None:
            reference = None
        else:
            reference, _ = self.process_mtzfile(
                params.input.hklref, merge_equivalents=True
            )

        # Configure logging
        log.config(info="dials.potato.log", debug="dials.potato.debug.log")

        # Log the diff phil
        diff_phil = self.parser.diff_phil.as_str()
        if diff_phil is not "":
            logger.info("The following parameters have been modified:\n")
            logger.info(diff_phil)

        # Plot the local cc
        if reference is None:
            plot_local_cc_half(
                reflections,
                A,
                kernel_size=params.kernel_size,
                record=params.output.record,
                directory=params.output.directory,
                point_size=params.point_size,
            )
        else:
            plot_local_cc_vs_ref(
                reflections,
                reference,
                A,
                kernel_size=params.kernel_size,
                record=params.output.record,
                directory=params.output.directory,
                point_size=params.point_size,
            )

    def process_reflections_and_experiments(
        self, reflections, experiments, which_intensity
    ):
        """
        Get the ASU miller indices and set an intensity and variance to use

        """
        # Get only those reflections that have been indexed
        if "intensity.scaled.value" in reflections:
            selection = reflections.get_flags(reflections.flags.scaled)
            reflections["intensity"] = (
                reflections["intensity.scale.value"]
                / reflections["inverse_scale_factor"]
            )
            reflections["variance"] = (
                reflections["intensity.scale.variance"]
                / reflections["inverse_scale_factor"]
            )
        else:
            if which_intensity == "sum":
                selection = reflections.get_flags(reflections.flags.integrated_sum)
                reflections["intensity"] = reflections["intensity.sum.value"]
                reflections["variance"] = reflections["intensity.sum.variance"]
            elif which_intensity == "prf":
                selection = reflections.get_flags(reflections.flags.integrated_prf)
                reflections["intensity"] = reflections["intensity.prf.value"]
                reflections["variance"] = reflections["intensity.prf.variance"]
            else:
                raise RuntimeError("Unknown intensity: %s" % which_intensity)
        reflections = reflections.select(selection)

        # Get the ASU miller index
        A = experiments[0].crystal.get_B()
        cs = experiments[0].crystal.get_crystal_symmetry()
        ms = miller.set(cs, reflections["miller_index"], anomalous_flag=False)
        ms_asu = ms.map_to_asu()
        reflections["asu_miller_index"] = ms_asu.indices()

        # Return the reflections
        return reflections, A

    def process_mtzfile(self, filename, merge_equivalents=False):
        """
        Process the mtz file

        """
        reader = any_reflection_file(filename)

        # Get the columns as miller arrays
        miller_arrays = reader.as_miller_arrays(merge_equivalents=merge_equivalents)

        # Select the desired columns
        intensities = None
        miller_set = None
        for array in miller_arrays:
            if array.info().labels == ["I", "SIGI"]:
                intensities = array
                miller_set = array.set(anomalous_flag=False)
        if intensities is None and merge_equivalents == True:
            for array in miller_arrays:
                if array.info().labels == ["IMEAN", "SIGIMEAN"]:
                    intensities = array
                    miller_set = array.set(anomalous_flag=False)
        if intensities is None:
            raise RuntimeError("No I/SIGI column in mtz file")

        # Create a reflection table
        reflections = flex.reflection_table()
        if not merge_equivalents:
            reflections["miller_indices"] = intensities.indices()
            reflections["intensity"] = intensities.data()
            reflections["variance"] = intensities.sigmas() ** 2
            miller_set_asu = miller_set.map_to_asu()
            reflections["asu_miller_index"] = miller_set_asu.indices()
        else:
            reflections["asu_miller_index"] = intensities.indices()
            reflections["intensity"] = intensities.data()
            reflections["variance"] = intensities.sigmas() ** 2

        A = matrix.col(
            miller_set.crystal_symmetry().unit_cell().fractionalization_matrix()
        ).transpose()

        # Return reflection table
        return reflections, A


if __name__ == "__main__":
    from dials.util import halraiser

    try:
        script = Script()
        script.run()
    except Exception as e:
        halraiser(e)
