"""
Simple script to run the stills indexer on the spotfinding results from a
still sequence i.e. SSX data.

Saves the indexed result into a single experiment list/reflection table
with a joint detector and beam model.
"""

import logging
import sys, os
from io import StringIO

from libtbx import phil

from dials.util import log, show_mail_handle_errors
from dials.util.options import OptionParser, reflections_and_experiments_from_files
from dials.util.version import dials_version

from dxtbx.model import ExperimentList
from dials.array_family import flex
from dials.algorithms.indexing.indexer import Indexer
try:
    from typing import List
except ImportError:
    pass

logger = logging.getLogger("dials")

from dials.command_line.index import phil_scope
from dials.algorithms.indexing import DialsIndexError
from dials.util.log import LoggingContext


program_defaults_phil_str = """
indexing {
  method = fft1d
}
refinement {
  parameterisation {
    auto_reduction {
      min_nref_per_parameter = 1
      action = fix
    }
    beam.fix = all
    detector.fix = all
  }
  reflections {
    weighting_strategy.override = stills
    outlier.algorithm = null
  }
}
"""

phil_scope = phil.parse("""
include scope dials.command_line.index.phil_scope
""", process_includes=True).fetch(
    phil.parse(program_defaults_phil_str)
)

def index(experiments, observed, params):
    params.refinement.parameterisation.scan_varying = False
    params.indexing.stills.indexer = "stills"

    reflections = observed.split_by_experiment_id()
    indexed_experiments = ExperimentList()
    indexed_reflections = flex.reflection_table()

    def run_with_disabled_logs(fn, fnargs):
        sys.stdout = open(os.devnull, 'w') # block printing from rstbx
        log1 = logging.getLogger("dials.algorithms.refinement.reflection_manager")
        log2 = logging.getLogger("dials.algorithms.refinement.refiner")
        log3 = logging.getLogger("dials.algorithms.indexing.stills_indexer")
        log4 = logging.getLogger("dials.algorithms.indexing.nave_parameters")
        with LoggingContext(log1, level=logging.WARNING):
            with LoggingContext(log2, level=logging.WARNING):
                with LoggingContext(log3, level=logging.WARNING):
                    with LoggingContext(log4, level=logging.WARNING):
                        fn(*fnargs)
        sys.stdout = sys.__stdout__ # restore printing

    def index_all(experiments, reflections, params):
        n_found = 0
        for i, (expt, refl) in enumerate(zip(experiments, reflections)):
            elist = ExperimentList([expt])
            refl["imageset_id"] = flex.int(refl.size(), 0) #needed for centroid_px_to_mm
            refl.centroid_px_to_mm(elist)
            refl.map_centroids_to_reciprocal_space(elist)
            idxr = Indexer.from_parameters(
                refl,
                elist,
                params=params,
            )
            try:
                idxr.index()
            except DialsIndexError as e:
                logger.info(f"Image {i+1}: Failed to index, error: {e}")
            else:
                #renumber numerical ids in the table
                ids_map = dict(idxr.refined_reflections.experiment_identifiers())
                for k in idxr.refined_reflections.experiment_identifiers().keys():
                    del idxr.refined_reflections.experiment_identifiers()[k]
                idxr.refined_reflections["id"] += n_found
                for k,v in ids_map.items():
                    idxr.refined_reflections.experiment_identifiers()[k+n_found] = v
                n_found += len(ids_map.keys())
                logger.info(f"Image {i+1}: Indexed {idxr.refined_reflections.size()} spots")
                indexed_reflections.extend(idxr.refined_reflections)
                indexed_experiments.extend(idxr.refined_experiments)


    run_with_disabled_logs(index_all, (experiments, reflections, params))
    logger.info(f"{indexed_reflections.size()} spots indexed across {len(indexed_experiments)} images")

    from dials.command_line.combine_experiments import CombineWithReference
    # combine beam and detector models if not already
    if (len(indexed_experiments.detectors())) > 1 or (len(indexed_experiments.beams())) > 1:
        combine = CombineWithReference(
            detector=indexed_experiments[0].detector,
            beam=indexed_experiments[0].beam,
        )
        elist = ExperimentList()
        for expt in indexed_experiments:
            elist.append(combine(expt))
        indexed_experiments = elist
    expt_filename = "indexed.expt"
    refl_filename = "indexed.refl"
    logger.info(f"Saving indexed experiments to {expt_filename}")
    indexed_experiments.as_file(expt_filename)
    logger.info(f"Saving indexed reflections to {refl_filename}")
    indexed_reflections.as_file(refl_filename)

@show_mail_handle_errors()
def run(args: List[str] = None, phil: phil.scope = phil_scope) -> None:
    parser = OptionParser(
        usage="",
        read_experiments=True,
        read_reflections=True,
        phil=phil_scope,
        check_format=False,
        epilog=__doc__,
    )
    params, options = parser.parse_args(args=args, show_diff_phil=False)

    if not params.input.experiments or not params.input.reflections:
        parser.print_help()
        sys.exit()

    reflections, experiments = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )
    log.config(verbosity=1, logfile="dials.ssx_index.log")
    logger.info(dials_version())

    index(experiments, reflections[0], params)

if __name__ == "__main__":
    run()

