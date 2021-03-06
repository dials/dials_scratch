#!/usr/bin/env python
#
# dials.merge_stills.py
#
#  Copyright (C) 2018 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

# LIBTBX_SET_DISPATCHER_NAME dials.merge_stills

from __future__ import absolute_import, division
import libtbx.load_env
from libtbx.utils import Sorry
from libtbx.phil import parse
from iotbx.reflection_file_reader import any_reflection_file
from dials.util import log
from dials.util.options import OptionParser
from dials.util.options import flatten_reflections, flatten_experiments
from dials_scratch.jmp.merge.merge import scale_and_merge
from dials.array_family import flex
import logging


logger = logging.getLogger(libtbx.env.dispatcher_name)

help_message = """

This script does merging for stills

"""

# Create the phil scope
phil_scope = parse(
    """

  input {
    reference = None
      .type = str
      .help = "The reference reflections"
  }

  include scope dials_scratch.jmp.merge.merge.phil_scope

""",
    process_includes=True,
)


class Script(object):
    """
    The integration program.

    """

    def __init__(self):
        """
        Initialise the script.

        """

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
        """
        Perform the integration.

        """
        from time import time

        # Check the number of arguments is correct
        start_time = time()

        # Parse the command line
        params, options = self.parser.parse_args(show_diff_phil=False)
        reflections = flatten_reflections(params.input.reflections)
        experiments = flatten_experiments(params.input.experiments)
        if len(reflections) == 0 or len(experiments) == 0:
            self.parser.print_help()
            return
        elif len(reflections) != 1:
            raise Sorry("more than 1 reflection file was given")
        elif len(experiments) == 0:
            raise Sorry("no experiment list was specified")
        reflections = reflections[0]

        # Configure logging
        log.config(info="dials.merge_stills.log", debug="dials.merge_stills.debug.log")

        # Log the diff phil
        diff_phil = self.parser.diff_phil.as_str()
        if diff_phil is not "":
            logger.info("The following parameters have been modified:\n")
            logger.info(diff_phil)

        # Read the reference
        if params.input.reference:
            reference = self._read_reference(params.input.reference)
        else:
            reference = None

        # Do the merging
        reflections = scale_and_merge(experiments, reflections, params, reference)

    def _read_reference(self, filename):

        # Read the MTZ file
        reader = any_reflection_file(filename)

        # Get the columns as miller arrays
        miller_arrays = reader.as_miller_arrays(merge_equivalents=False)

        # Select the desired columns
        intensities = None
        for array in miller_arrays:
            if array.info().labels == ["IMEAN", "SIGIMEAN"]:
                intensities = array
        assert intensities is not None
        indices = intensities.indices()
        I = intensities.data()
        V = intensities.sigmas() ** 2
        reference = flex.reflection_table()
        reference["miller_index"] = indices
        reference["intensity.ref.value"] = I
        reference["intensity.ref.variance"] = V
        return reference


if __name__ == "__main__":
    from dials.util import halraiser

    try:
        script = Script()
        script.run()
    except Exception as e:
        halraiser(e)
