"""
Simple script to run the stills indexer on the spotfinding results from a
still sequence i.e. SSX data.

Saves the indexed result into a single experiment list/reflection table
with a joint detector and beam model.
"""

import logging
import sys, os
import math
from io import StringIO

from libtbx import phil

from dials.util import log, show_mail_handle_errors, tabulate
from dials.util.options import OptionParser, reflections_and_experiments_from_files
from dials.util.version import dials_version

from dxtbx.model import ExperimentList, Experiment
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

RAD2DEG = 180 / math.pi

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

phil_scope = phil.parse(
    """
include scope dials.command_line.index.phil_scope
""",
    process_includes=True,
).fetch(phil.parse(program_defaults_phil_str))


def index(experiments, observed, params):
    params.refinement.parameterisation.scan_varying = False
    params.indexing.stills.indexer = "stills"

    reflections = observed.split_by_experiment_id()
    indexed_experiments = ExperimentList()
    indexed_reflections = flex.reflection_table()

    def run_with_disabled_logs(fn, fnargs):
        sys.stdout = open(os.devnull, "w")  # block printing from rstbx
        log1 = logging.getLogger("dials.algorithms.refinement.reflection_manager")
        log2 = logging.getLogger("dials.algorithms.refinement.refiner")
        log3 = logging.getLogger("dials.algorithms.indexing.stills_indexer")
        log4 = logging.getLogger("dials.algorithms.indexing.nave_parameters")
        log5 = logging.getLogger("dials.algorithms.indexing.basis_vector_search.real_space_grid_search")
        log6 = logging.getLogger("dials.algorithms.indexing.basis_vector_search.combinations")
        with LoggingContext(log1, level=logging.ERROR):
            with LoggingContext(log2, level=logging.ERROR):
                with LoggingContext(log3, level=logging.ERROR):
                    with LoggingContext(log4, level=logging.ERROR):
                        with LoggingContext(log5, level=logging.ERROR):
                            with LoggingContext(log6, level=logging.ERROR):
                                fn(*fnargs)
        sys.stdout = sys.__stdout__  # restore printing

    if params.indexing.known_symmetry.unit_cell:
        method_list = ["fft1d", "real_space_grid_search"]
    else:
        method_list = [params.indexing.method]
    methods = ", ".join(method_list)
    pl = "s" if (len(method_list) > 1) else ""
    def plural(method):
        return f" with {method} method" if (len(method_list) > 1) else ""

    logger.info(f"Attempting indexing with {methods} method{pl}")

    def index_all(experiments, reflections, params):
        n_found = 0
        overall_summary_header = ["Image", "lattice", "expt_id", "n_indexed", "dX", "dY", "dPsi"]
        overall_summary_rows = []
        # toggle on n_lattice and method depending on if multiple or not?
        for i, (expt, refl) in enumerate(zip(experiments, reflections)):
            ### Calculate necessary quantities
            elist = ExperimentList([expt])
            refl["imageset_id"] = flex.int(
                refl.size(), 0
            )  # needed for centroid_px_to_mm
            refl.centroid_px_to_mm(elist)
            refl.map_centroids_to_reciprocal_space(elist)
            path = expt.imageset.get_path(i)

            ### Loop over methods and attempt indexing
            for method in method_list:
                params.indexing.method = method
                idxr = Indexer.from_parameters(refl, elist, params=params)
                try:
                    idxr.index()
                except DialsIndexError as e:
                    logger.info(f"Image {i+1}: Failed to index{plural(method)}, error: {e}")
                else:
                    ### Extract rmsd values for output
                    ids_map = dict(idxr.refined_reflections.experiment_identifiers())
                    logger.info(f"Image {i+1}: Indexed with {method} method.")
                    header = ["Crystal #", "n indexed", "RMSD_X (px)", "RMSD_Y (px)", "RMSD_DeltaPsi (deg)"]
                    rows = []
                    n_strong = refl.size()
                    for n_cryst, id_ in enumerate(ids_map.keys()):
                        selr = idxr.refined_reflections.select(idxr.refined_reflections["id"] == id_)
                        calx, caly, calz = selr["xyzcal.px"].parts()
                        obsx, obsy, obsz = selr["xyzobs.px.value"].parts()
                        delpsi = selr["delpsical.rad"]
                        rmsd_x = flex.mean((calx - obsx) ** 2)**0.5
                        rmsd_y = flex.mean((caly - obsy) ** 2)**0.5
                        rmsd_z = flex.mean(((delpsi) * RAD2DEG) ** 2)**0.5
                        n_id_ = calx.size()
                        n_indexed = f"{n_id_}/{n_strong} ({100*n_id_/n_strong:2.1f}%)"
                        rows.append([
                            f"{n_cryst + 1}",
                            n_indexed,
                            f"{rmsd_x:.3f}",
                            f"{rmsd_y:.3f}",
                            f" {rmsd_z:.4f}",
                        ])
                        overall_summary_rows.append([
                            path.split("/")[-1],
                            n_cryst + 1,
                            n_found + n_cryst,
                            n_indexed,
                            f"{rmsd_x:.3f}",
                            f"{rmsd_y:.3f}",
                            f" {rmsd_z:.4f}",
                        ])
                    logger.info(tabulate(rows, header))

                    ### renumber numerical ids in the table
                    for k in idxr.refined_reflections.experiment_identifiers().keys():
                        del idxr.refined_reflections.experiment_identifiers()[k]
                    idxr.refined_reflections["id"] += n_found
                    for k, v in ids_map.items():
                        idxr.refined_reflections.experiment_identifiers()[k + n_found] = v
                    n_found += len(ids_map.keys())

                    ### Add data to overall list/table
                    indexed_reflections.extend(idxr.refined_reflections)
                    for jexpt in idxr.refined_experiments:
                        indexed_experiments.append(
                            Experiment(
                                identifier=jexpt.identifier,
                                beam=jexpt.beam,
                                detector=jexpt.detector,
                                scan=jexpt.scan,
                                goniometer=jexpt.goniometer,
                                crystal=jexpt.crystal,
                                imageset=jexpt.imageset[i:i+1],
                            )
                        )
                    break

        # now give overall summary
        logger.info(tabulate(overall_summary_rows, overall_summary_header))


    run_with_disabled_logs(index_all, (experiments, reflections, params))
    n_images = len(set(e.imageset.get_path(0) for e in indexed_experiments))
    logger.info(
        f"{indexed_reflections.size()} spots indexed on {n_images} images"
    )

    from dials.command_line.combine_experiments import CombineWithReference

    # combine beam and detector models if not already
    if (len(indexed_experiments.detectors())) > 1 or (
        len(indexed_experiments.beams())
    ) > 1:
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

    diff_phil = parser.diff_phil.as_str()
    if diff_phil:
        logger.info("The following parameters have been modified:\n%s", diff_phil)

    index(experiments, reflections[0], params)


if __name__ == "__main__":
    run()
