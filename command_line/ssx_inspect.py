"""Inspect integrated ssx data"""
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

phil_scope = phil.parse("""
d_min = None
    .type = float
threshold_i_over_sigma = None
    .type = float
""")

def inspect(experiments, reflections, params):

    Ioversig_per_crystal = flex.double()
    if params.d_min:
        reflections = reflections.select(reflections["d"] > params.d_min)
        if not reflections:
            return
    I = reflections["intensity.sum.value"]
    sig = reflections["intensity.sum.variance"] ** 0.5
    Ioversig = I/sig
    images = []

    for (id_, identifier) in dict(reflections.experiment_identifiers()).items():
        sel = reflections["id"] == id_
        Ioversig_per_crystal.append(flex.mean(Ioversig.select(sel)))
        expt = (experiments.identifiers() == identifier).iselection()[0]
        images.append(experiments[expt].imageset.get_path(0).split("/")[-1])

    header = ["image", "expt_id", ("I/sigma" + f" (d_min={params.d_min})" if params.d_min else "I/sigma")]
    rows = []

    for i, (image, Iovers) in enumerate(zip(images, Ioversig_per_crystal)):
        rows.append([f"{image}", f"{i}", f"{Iovers:.2f}"])

    logger.info(tabulate(rows, header))

    if params.threshold_i_over_sigma:
        logger.info(f"Integrated images with I/sigma > {params.threshold_i_over_sigma}")
        header = ["image", "expt_id", ("I/sigma" + f" (d_min={params.d_min})" if params.d_min else "I/sigma")]
        rows = []

        for i, (image, Iovers) in enumerate(zip(images, Ioversig_per_crystal)):
            if Iovers > params.threshold_i_over_sigma:
                rows.append([f"{image}", f"{i}", f"{Iovers:.2f}"])

        logger.info(tabulate(rows, header))


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
    log.config(verbosity=1, logfile="dials.ssx_inspect.log")
    logger.info(dials_version())

    diff_phil = parser.diff_phil.as_str()
    if diff_phil:
        logger.info("The following parameters have been modified:\n%s", diff_phil)

    inspect(experiments, reflections[0], params)


if __name__ == "__main__":
    run()
