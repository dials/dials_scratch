#!/usr/bin/env dials.python

#
#  Copyright (C) (2017) Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

# LIBTBX_SET_DISPATCHER_NAME dev.dials.show_mask

"""
Just display the mask for the specified image

"""

from __future__ import absolute_import, division, print_function

from dials.util.options import flatten_experiments
from libtbx.phil import parse
from dials.util import show_mail_on_error

# The phil scope
phil_scope = parse(
    """

  image = 0
    .type = int
    .help = "Which image to show"

""",
    process_includes=True,
)

help_message = """

Utility to just display the mask for the desired image

Examples::

  dev.dials.show_mask models.expt image=1

"""


class Script(object):
    def __init__(self):
        """Initialise the script."""
        from dials.util.options import OptionParser
        import libtbx.load_env

        # The script usage
        usage = (
            "usage: %s [options] [param.phil] "
            "experiments1.expt experiments2.expt reflections1.refl "
            "reflections2.refl..." % libtbx.env.dispatcher_name
        )

        # Create the parser
        self.parser = OptionParser(
            usage=usage, phil=phil_scope, read_experiments=True, epilog=help_message
        )

    def run(self):
        """Execute the script."""

        # Parse the command line
        params, options = self.parser.parse_args(show_diff_phil=True)

        experiments = flatten_experiments(params.input.experiments)
        assert len(experiments) == 1
        imageset = experiments[0].imageset

        mask = imageset.get_mask(params.image)

        assert len(mask) == 1

        print("Num True: %d" % mask[0].count(True))
        print("Num False: %d" % mask[0].count(False))

        from matplotlib import pylab

        pylab.imshow(mask[0].as_numpy_array(), interpolation="none")
        pylab.show()


if __name__ == "__main__":
    with show_mail_on_error():
        script = Script()
        script.run()
