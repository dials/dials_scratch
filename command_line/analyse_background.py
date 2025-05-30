#!/usr/bin/env python
#
# dials.analyse_background.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import absolute_import, division, print_function

from matplotlib import pylab

from dials.util import show_mail_on_error
from dials.util.options import flatten_experiments, OptionParser
import libtbx.load_env
from libtbx.phil import parse

# DIALS_ENABLE_COMMAND_LINE_COMPLETION
# LIBTBX_SET_DISPATCHER_NAME dev.dials.analyse_background

help_message = ""

# Create the phil scope
phil_scope = parse("", process_includes=True)


class Script(object):
    """The integration program."""

    def __init__(self):
        """Initialise the script."""
        # The script usage
        usage = "usage: %s [options] experiment.expt" % libtbx.env.dispatcher_name

        # Create the parser
        self.parser = OptionParser(
            usage=usage, phil=phil_scope, epilog=help_message, read_experiments=True
        )

    def run(self):
        """Analyse the background"""
        # Parse the command line
        params, options = self.parser.parse_args(show_diff_phil=False)
        experiments = flatten_experiments(params.input.experiments)
        if len(experiments) == 0:
            self.parser.print_help()
            return

        assert len(experiments) == 1

        # Get the imageset
        imageset = experiments[0].imageset

        total_image = None
        for i in range(len(imageset)):
            print(i)
            image = imageset.get_raw_data(i)
            if total_image is None:
                total_image = image[0]
            else:
                total_image += image[0]
        total_image /= len(imageset)
        print(min(total_image))
        print(max(total_image))
        print(sum(total_image) / len(total_image))

        pylab.imshow(total_image.as_numpy_array(), vmin=0, vmax=2)
        pylab.show()


if __name__ == "__main__":
    with show_mail_on_error():
        script = Script()
        script.run()
