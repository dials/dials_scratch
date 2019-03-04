#!/usr/bin/env python
#
# dials.integrate_still.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import absolute_import, division
from __future__ import print_function
import libtbx.load_env
import logging

logger = logging.getLogger(libtbx.env.dispatcher_name)

help_message = """

"""

# Create the phil scope
from libtbx.phil import parse

phil_scope = parse(
    """

  output {

    experiments = 'integrated_experiments.json'
      .type = str
      .help = "The experiments output filename"

    reflections = 'integrated.pickle'
      .type = str
      .help = "The integrated output filename"

    model = 'model.pickle'
      .type = str
      .help = "The image model filename"

    model_image = 'model.png'
      .type = str
      .help = "The image model filename"

  }

  include scope dials_scratch.jmp.stills.custard.phil_scope

""",
    process_includes=True,
)


class Script(object):
    """ The integration program. """

    def __init__(self):
        """Initialise the script."""
        from dials.util.options import OptionParser
        import libtbx.load_env

        # The script usage
        usage = "usage: %s [options] experiment.json" % libtbx.env.dispatcher_name

        # Create the parser
        self.parser = OptionParser(
            usage=usage, phil=phil_scope, epilog=help_message, read_experiments=True
        )

    def run(self):
        """ Perform the integration. """

        from dials.util.options import flatten_reflections, flatten_experiments
        from dxtbx.model.experiment_list import ExperimentListDumper
        from dials_scratch.jmp.stills.custard import Integrator
        from dials.array_family import flex
        from scitbx import matrix
        import cPickle as pickle
        import sys

        # Parse the command line
        params, options = self.parser.parse_args(show_diff_phil=False)
        experiments = flatten_experiments(params.input.experiments)
        assert len(experiments) == 1

        # Do the integration
        integrator = Integrator(experiments[0], params)
        integrator.process()

        # Get the results
        reflections = integrator.reflections
        experiments[0] = integrator.experiment

        # Get the predicted image model
        image_pred = integrator.image_pred

        # Compute min/max partiality
        selection = reflections.get_flags(reflections.flags.integrated_sum)
        partiality = reflections["partiality"].select(selection)
        min_partiality = flex.min(partiality)
        max_partiality = flex.max(partiality)

        # Get the mosaicity
        mosaicity = integrator.mosaicity

        print("")
        print("Mosaicity: %f" % mosaicity)
        print("")

        # The mosaicity matrix in q-space
        Mq = matrix.sqr((mosaicity, 0, 0, 0, mosaicity, 0, 0, 0, mosaicity))
        print("Mosacity matrix in q-space")
        print(Mq.as_numpy_array())
        print("")

        # The mosaicity matrix in hkl
        A = matrix.sqr(experiments[0].crystal.get_A())
        Mh = A.inverse() * Mq
        print("Mosacity matrix in h-space")
        print(Mh.as_numpy_array())
        print("")

        # Print partiality
        print(
            "Min partiality: %f, Max partiality: %f" % (min_partiality, max_partiality)
        )
        print("")

        # Save the reflections
        reflections.as_pickle(params.output.reflections)

        # Save the experiments
        ExperimentListDumper(experiments).as_json(params.output.experiments)

        # Save the image model
        pickle.dump(image_pred, open(params.output.model, "w"))

        # Save an image of the image model
        from matplotlib import pylab

        fig, ax = pylab.subplots(nrows=1, ncols=1)
        ax.imshow(image_pred.as_numpy_array(), interpolation="none")
        fig.savefig(params.output.model_image, dpi=600)
        pylab.close(fig)


if __name__ == "__main__":
    from dials.util import halraiser

    try:
        script = Script()
        script.run()
    except Exception as e:
        halraiser(e)
