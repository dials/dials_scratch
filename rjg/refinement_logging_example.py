from __future__ import absolute_import, division, print_function

import logging
import sys

from libtbx import phil

from dxtbx.serialize import load
from dials.array_family import flex
from dials.algorithms.refinement import RefinerFactory
from dials.algorithms.refinement.refiner import phil_scope


def run(args):
    assert len(args) == 2

    expts = load.experiment_list(args[0], check_format=False)
    refl = flex.reflection_table.from_file(args[1])
    params = phil_scope.fetch(source=phil.parse("")).extract()
    params.refinement.verbosity = 12

    # no output by default
    print("Starting refinement #1")
    refiner = RefinerFactory.from_parameters_data_experiments(params, refl, expts)
    refiner.run()
    print("Finished refinement #1")

    # configure logging
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )

    # we should now get logging output from refinement
    print("Starting refinement #2")
    refiner = RefinerFactory.from_parameters_data_experiments(params, refl, expts)
    refiner.run()
    print("Finished refinement #2")

    # switch off logging for dials.algorithms.refinement
    logging.getLogger("dials.algorithms.refinement").setLevel(logging.ERROR)
    print("Starting refinement #3")
    refiner = RefinerFactory.from_parameters_data_experiments(params, refl, expts)
    refiner.run()
    print("Finished refinement #3")


if __name__ == "__main__":
    run(sys.argv[1:])
