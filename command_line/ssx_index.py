"""
Simple script to run the stills indexer on the spotfinding results from a
still sequence i.e. SSX data.

Saves the indexed result into a single experiment list/reflection table
with a joint detector and beam model.
"""

import logging
import time
import sys, os
import math
import numpy as np
from io import StringIO
import concurrent.futures
from collections import defaultdict

from libtbx import phil, Auto

from dials.util import log, show_mail_handle_errors, tabulate
from dials.util.options import OptionParser, reflections_and_experiments_from_files
from dials.util.version import dials_version

from dxtbx.model import ExperimentList, Experiment
from dials.array_family import flex
from dials.algorithms.indexing.indexer import Indexer
from dials.algorithms.indexing.max_cell import find_max_cell

from xfel.clustering.cluster import Cluster
from xfel.clustering.cluster_groups import unit_cell_info
from cctbx import crystal
from dials.command_line.combine_experiments import CombineWithReference
from dials.util.mp import multi_node_parallel_map

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
output.log = dials.ssx_index.log
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
method = *fft1d *real_space_grid_search
    .type = choice(multi=True)
include scope dials.command_line.index.phil_scope
""",
    process_includes=True,
).fetch(phil.parse(program_defaults_phil_str))


def _index_one(experiment, refl, params, method_list, expt_no):
    elist = ExperimentList([experiment])
    for method in method_list:
        params.indexing.method = method
        idxr = Indexer.from_parameters(refl, elist, params=params)
        try:
            idxr.index()
        except (DialsIndexError, AssertionError) as e:
            logger.info(
                f"Image {expt_no+1}: Failed to index with {method} method, error: {e}"
            )
            if method == method_list[-1]:
                return None, None, expt_no
        else:
            logger.info(
                f"Image {expt_no+1}: Indexed {idxr.refined_reflections.size()}/{refl.size()} spots with {method} method."
            )
            return idxr.refined_experiments, idxr.refined_reflections, expt_no


def run_with_disabled_logs(fn, fnargs):
    sys.stdout = open(os.devnull, "w")  # block printing from rstbx
    log1 = logging.getLogger("dials.algorithms.refinement.reflection_processor")
    log2 = logging.getLogger("dials.algorithms.refinement.refiner")
    log8 = logging.getLogger("dials.algorithms.refinement.reflection_manager")
    log3 = logging.getLogger("dials.algorithms.indexing.stills_indexer")
    log4 = logging.getLogger("dials.algorithms.indexing.nave_parameters")
    log5 = logging.getLogger(
        "dials.algorithms.indexing.basis_vector_search.real_space_grid_search"
    )
    log6 = logging.getLogger(
        "dials.algorithms.indexing.basis_vector_search.combinations"
    )
    log7 = logging.getLogger("dials.algorithms.indexing.indexer")
    with LoggingContext(log1, level=logging.ERROR):
        with LoggingContext(log2, level=logging.ERROR):
            with LoggingContext(log3, level=logging.ERROR):
                with LoggingContext(log4, level=logging.ERROR):
                    with LoggingContext(log5, level=logging.ERROR):
                        with LoggingContext(log6, level=logging.ERROR):
                            with LoggingContext(log7, level=logging.ERROR):
                                with LoggingContext(log8, level=logging.ERROR):
                                    result = fn(*fnargs)
                                    sys.stdout = sys.__stdout__  # restore printing
                                    return result


class Processor(object):
    # Wrap some functions into a class to allow multiprocessing with
    # multi_node_parallel_map

    def __init__(self, experiments, reflections, params):
        self.experiments = experiments
        self.reflections = reflections
        self.params = params
        self._results_order = np.array([], dtype=np.int32)
        self._results = defaultdict(list)
        self._n_strong = np.array([table.size() for table in self.reflections])
        self._n_found = 0
        self._tables_list = []
        self._expts_list = []
        self.indexed_experiments = ExperimentList()
        self.indexed_reflections = flex.reflection_table()
        self.summary_table = ""

    def process_output(self, result):
        idx_expts, idx_refl, index = result[0], result[1], result[2]
        if idx_expts:
            self._results_order = np.append(self._results_order, [index])
            ids_map = dict(idx_refl.experiment_identifiers())
            path = idx_expts[0].imageset.get_path(index)
            for n_cryst, id_ in enumerate(ids_map.keys()):
                selr = idx_refl.select(idx_refl["id"] == id_)
                calx, caly, calz = selr["xyzcal.px"].parts()
                obsx, obsy, obsz = selr["xyzobs.px.value"].parts()
                delpsi = selr["delpsical.rad"]
                rmsd_x = flex.mean((calx - obsx) ** 2) ** 0.5
                rmsd_y = flex.mean((caly - obsy) ** 2) ** 0.5
                rmsd_z = flex.mean(((delpsi) * RAD2DEG) ** 2) ** 0.5
                n_id_ = calx.size()
                n_indexed = f"{n_id_}/{self._n_strong[index]} ({100*n_id_/self._n_strong[index]:2.1f}%)"
                self._results[index].append(
                    [
                        path.split("/")[-1],
                        n_indexed,
                        f"{rmsd_x:.3f}",
                        f"{rmsd_y:.3f}",
                        f" {rmsd_z:.4f}",
                    ]
                )
            self._n_found += len(ids_map.keys())
            self._tables_list.append(idx_refl)
            elist = ExperimentList()
            for jexpt in idx_expts:
                elist.append(
                    Experiment(
                        identifier=jexpt.identifier,
                        beam=jexpt.beam,
                        detector=jexpt.detector,
                        scan=jexpt.scan,
                        goniometer=jexpt.goniometer,
                        crystal=jexpt.crystal,
                        imageset=jexpt.imageset[index : index + 1],
                    )
                )
            self._expts_list.append(elist)

    def process(self, method_list):
        inputs = []
        for i, (expt, refl) in enumerate(zip(self.experiments, self.reflections)):
            inputs.append((expt, refl, self.params, method_list, i))
        mp_nproc = self.params.indexing.nproc
        mp_njobs = 1
        mp_method = "multiprocessing"

        def execute_task(input_):
            return _index_one(*input_)

        if mp_nproc > 1:
            multi_node_parallel_map(
                func=execute_task,
                iterable=inputs,
                njobs=mp_njobs,
                nproc=mp_nproc,
                callback=self.process_output,
                cluster_method=mp_method,
                preserve_order=True,
            )
        else:
            for input_ in inputs:
                self.process_output(execute_task(input_))

    def finalize(self):
        # Results came in a non-deterministic order. So sort them out
        # to match the input order, adjusting the ids as appropriate.
        n_tot = 0
        for i in np.argsort(self._results_order):
            table = self._tables_list[i]
            expts = self._expts_list[i]
            self.indexed_experiments.extend(expts)
            ids_map = dict(table.experiment_identifiers())
            for k in table.experiment_identifiers().keys():
                del table.experiment_identifiers()[k]
            table["id"] += n_tot
            for k, v in ids_map.items():
                table.experiment_identifiers()[k + n_tot] = v
            n_tot += len(ids_map.keys())
            self.indexed_reflections.extend(table)
        self.indexed_reflections.assert_experiment_identifiers_are_consistent(
            self.indexed_experiments
        )

        # Add a few extra useful items to the summary table.
        overall_summary_header = [
            "Image",
            "expt_id",
            "n_indexed",
            "RMSD_X",
            "RMSD_Y",
            "RMSD_dPsi",
        ]

        rows = []
        total = 0
        if self.params.indexing.multiple_lattice_search.max_lattices > 1:
            show_lattices = True
            overall_summary_header.insert(1, "lattice")
        else:
            show_lattices = False
        for i, k in enumerate(sorted(self._results.keys())):
            for j, cryst in enumerate(self._results[k]):
                cryst.insert(1, total)
                if show_lattices:
                    cryst.insert(1, j + 1)
                rows.append(cryst)
                total += 1

        self.summary_table = tabulate(rows, overall_summary_header)


def index(experiments, observed, params):
    params.refinement.parameterisation.scan_varying = False
    params.indexing.stills.indexer = "stills"

    reflections = observed.split_by_experiment_id()
    # Calculate necessary quantities
    for refl, experiment in zip(reflections, experiments):
        elist = ExperimentList([experiment])
        refl["imageset_id"] = flex.int(refl.size(), 0)  # needed for centroid_px_to_mm
        refl.centroid_px_to_mm(elist)
        refl.map_centroids_to_reciprocal_space(elist)

    if (params.indexing.max_cell is Auto) and (
        not params.indexing.known_symmetry.unit_cell
    ):
        max_cells = []
        for refl in reflections:
            try:
                max_cells.append(find_max_cell(refl).max_cell)
            except (DialsIndexError, AssertionError):
                pass
        logger.info(f"Setting max cell to {max(max_cells):.1f} " + "\u212B")
        params.indexing.max_cell = max(max_cells)

    method_list = params.method
    if "real_space_grid_search" in method_list:
        if not params.indexing.known_symmetry.unit_cell:
            logger.info("No unit cell given, real_space_grid_search will not be used")
            method_list.remove("real_space_grid_search")
    methods = ", ".join(method_list)
    pl = "s" if (len(method_list) > 1) else ""
    logger.info(f"Attempting indexing with {methods} method{pl}")

    def index_all(experiments, reflections, params):
        processor = Processor(experiments, reflections, params)
        processor.process(method_list)
        processor.finalize()
        return (
            processor.indexed_experiments,
            processor.indexed_reflections,
            processor.summary_table,
        )

    indexed_experiments, indexed_reflections, summary = run_with_disabled_logs(
        index_all, (experiments, reflections, params)
    )

    # print some clustering information
    crystal_symmetries = [
        crystal.symmetry(
            unit_cell=expt.crystal.get_unit_cell(),
            space_group=expt.crystal.get_space_group(),
        )
        for expt in indexed_experiments
    ]
    ucs = Cluster.from_crystal_symmetries(crystal_symmetries)
    clusters, cluster_axes = ucs.ab_cluster(
        5000, log=None, write_file_lists=False, doplot=False
    )
    logger.info("\nUnit cell clustering analysis\n" + unit_cell_info(clusters))
    logger.info("\nSummary of images sucessfully indexed\n" + summary)

    n_images = len(set(e.imageset.get_path(0) for e in indexed_experiments))
    logger.info(f"{indexed_reflections.size()} spots indexed on {n_images} images")

    # combine beam and detector models if not already
    if (len(indexed_experiments.detectors())) > 1 or (
        len(indexed_experiments.beams())
    ) > 1:
        combine = CombineWithReference(
            detector=indexed_experiments[0].detector, beam=indexed_experiments[0].beam
        )
        elist = ExperimentList()
        for expt in indexed_experiments:
            elist.append(combine(expt))
        indexed_experiments = elist

    logger.info(f"Saving indexed experiments to {params.output.experiments}")
    indexed_experiments.as_file(params.output.experiments)
    logger.info(f"Saving indexed reflections to {params.output.reflections}")
    indexed_reflections.as_file(params.output.reflections)


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
    log.config(verbosity=1, logfile=params.output.log)
    logger.info(dials_version())

    diff_phil = parser.diff_phil.as_str()
    if diff_phil:
        logger.info("The following parameters have been modified:\n%s", diff_phil)

    st = time.time()
    index(experiments, reflections[0], params)
    logger.info(f"Total time: {time.time() - st:.2f}s")


if __name__ == "__main__":
    run()
