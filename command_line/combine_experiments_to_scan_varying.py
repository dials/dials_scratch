#!/usr/bin/env python

# LIBTBX_SET_DISPATCHER_NAME dev.dials.combine_experiments_to_scan_varying

"""
Combine multiple experiments with narrow scans into a single experiment with
a wider scan. The static model for the resulting experiment will be taken from
the first experiment of the input. A scan-varying crystal model for the result
will be constructed by concatenating the static crystal models of the input
experiments. The scans of input experiments must be in order of increasing
rotation angle (and image number), with no gaps or overlaps between them.

This script supports pseudo-scan-varying refinement of a crystal model with
a protocol like this:

# Chunk a scan into equal-width blocks
dials.slice_sequence static.json static.pickle block_size=5.0

dials.combine_experiments static_sliced.json static_sliced.pickle \
  reference_from_experiment.beam=0 \
  reference_from_experiment.detector=0 \
  reference_from_experiment.goniometer=0

# Refine crystal model in static blocks 0f 5 deg with a single beam and detector
dials.refine combined_experiments.json combined_reflections.pickle

# Recombine refined result into a single discontinuously-varying model
dev.dials.combine_experiments_to_scan_varying \
  refined_experiments.json refined.pickle
"""

from __future__ import division, print_function, absolute_import
import sys
from libtbx.utils import Sorry
from dials.util.options import flatten_reflections, flatten_experiments, OptionParser
from libtbx.table_utils import simple_table
from scitbx import matrix
from dxtbx.model.experiment_list import Experiment
from dxtbx.model.experiment_list import ExperimentList
from libtbx.phil import parse
from copy import deepcopy

# The phil scope
phil_scope = parse(
    """
output {
  experiments = combined_to_sv_experiments.json
    .type = str
    .help = "The filename for the combined experiment"

  reflections = combined_to_sv_reflections.pickle
    .type = str
    .help = "The filename for combined reflections"
}
""",
    process_includes=True,
)


class Script(object):
    """A class for running the script."""

    def __init__(self):
        """Initialise the script."""
        import libtbx.load_env

        # The script usage
        import __main__

        usage = (
            "usage: dials.python {0} refined_experiments.json " "refined.pickle"
        ).format(__main__.__file__)

        # Create the parser
        self.parser = OptionParser(
            usage=usage,
            phil=phil_scope,
            read_experiments=True,
            read_reflections=True,
            check_format=False,
            epilog=__doc__,
        )

        return

    @staticmethod
    def extended_scan(experiments):
        """Return a single scan that encompasses the full scan width of all of
    the inputs"""
        scan = deepcopy(experiments[0].scan)

        for exp in experiments[1:]:
            scan.append(exp.scan)

        return scan

    @staticmethod
    def combine_crystals(experiments, scan):
        """Create a single crystal model using the static models of each of the
    input models to define a scan-varying A matrix"""

        crystal = deepcopy(experiments[0].crystal)
        crystal.reset_scan_points()

        # Extract A matrices from the experiments
        A = []
        for exp in experiments:
            A.extend([exp.crystal.get_A()] * exp.scan.get_num_images())
        assert len(A) == scan.get_num_images()

        # There should be one more scan point than the number of images in the
        # scan so that the A matrix is defined at the beginning and end of the
        # scan and at every image boundary. By adding to the end of the A list
        # the values at the discontinuities favour the later scan
        A += [A[-1]]

        # Set the batched static model as a 'scan-varying' model of the crystal
        crystal.set_A_at_scan_points(A)

        return crystal

    def run(self):
        """Execute the script."""

        # Parse the command line
        self.params, options = self.parser.parse_args(show_diff_phil=True)

        if not self.params.input.experiments:
            self.parser.print_help()
            sys.exit()

        # Try to load the models
        experiments = flatten_experiments(self.params.input.experiments)
        nexp = len(experiments)
        if nexp == 0:
            print("No Experiments found in the input")
            self.parser.print_help()
            return

        ref_beam = experiments[0].beam
        ref_goniometer = experiments[0].goniometer
        ref_detector = experiments[0].detector

        scan = self.extended_scan(experiments)

        crystal = self.combine_crystals(experiments, scan)

        experiment = Experiment(
            beam=ref_beam,
            detector=ref_detector,
            scan=scan,
            goniometer=ref_goniometer,
            crystal=crystal,
        )

        experiments = ExperimentList([experiment])

        # Reset experiment IDs in the reflections
        reflections = flatten_reflections(self.params.input.reflections)
        assert len(reflections) == 1
        reflections = reflections[0]
        reflections["id"] *= 0

        # Save the experiments to file
        print(
            "Saving the combined experiment to {0}".format(
                self.params.output.experiments
            )
        )
        from dxtbx.model.experiment_list import ExperimentListDumper

        dump = ExperimentListDumper(experiments)
        dump.as_json(self.params.output.experiments)

        # Save the reflections to file
        print(
            "Saving the combined reflections to {0}".format(
                self.params.output.reflections
            )
        )
        reflections.as_pickle(self.params.output.reflections)


if __name__ == "__main__":
    from dials.util import halraiser

    try:
        script = Script()
        script.run()
    except Exception as e:
        halraiser(e)
