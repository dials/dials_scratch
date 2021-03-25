# LIBTBX_SET_DISPATCHER_NAME dev.dials.bivoet_group
"""
Group reflections by symmetry equivalents and print in csv format

Usage: dials.python bijvoet_group scaled.expt scaled.refl
"""

import os
import sys
import logging
import csv
from dials.util.options import OptionParser, flatten_experiments, flatten_reflections
import libtbx.phil
import dials.util
import dials.util.log
from cctbx import miller

# Define a logger. __name__ cannot be used as this script is called directly.
# Omit the dev prefix to use the DIALS logger
logger = logging.getLogger("dials.bijvoet_group")

# No PHIL options for this script, yet
phil_scope = libtbx.phil.parse(
    """
output {
    log = dev.dials.bijvoet_group.log
        .type = path

    csv_file = bijvoet_group.csv
        .type = path
        .help = "Filename for the output table in CSV format. If there are"
                "multiple experiments present each will be written to an"
                "individual file with a serial number appended"
}
"""
)


def bijvoet_table(experiments, reflections, params):

    csv_file = params.output.csv_file
    for experiment in experiments:
        rt = reflections.select_on_experiment_identifiers(
            [
                experiment.identifier,
            ]
        )
        exp_id = list(set(rt["id"]))
        assert len(exp_id) == 1
        exp_id = exp_id[0]
        symm = experiment.crystal.get_crystal_symmetry()

        ms = miller.set(symm, anomalous_flag=False, indices=rt["miller_index"])
        rt["miller_index_asu"] = ms.map_to_asu().indices()
        rt.sort("miller_index_asu")

        if len(experiments) > 1:
            (root, ext) = os.path.splitext(params.output.csv_file)
            csv_file = root + f"_{exp_id}" + ext

        logger.info(
            f"Writing table for dataset: {experiment.imageset.get_template()} to {csv_file}"
        )
        with open(csv_file, "w", newline="") as f:
            table_writer = csv.writer(f)
            table_writer.writerow(
                ["ASU_HKL", "HKL", "X", "Y", "Z", "partiality", "intensity"]
            )
            for ref in rt.rows():
                line = [
                    " ".join(str(e) for e in ref["miller_index_asu"]),
                    " ".join(str(e) for e in ref["miller_index"]),
                    f"{ref['xyzcal.px'][0]:.3f}",
                    f"{ref['xyzcal.px'][1]:.3f}",
                    f"{ref['xyzcal.px'][2]:.3f}",
                    f"{ref['partiality']:.3f}",
                    f"{ref['intensity.scale.value']:.6f}",
                ]
                table_writer.writerow(line)


@dials.util.show_mail_on_error()
def run(args=None):
    usage = "dials.python bijvoet_group [options] scaled.expt scaled.refl"

    parser = OptionParser(
        usage=usage,
        phil=phil_scope,
        read_reflections=True,
        read_experiments=True,
        check_format=False,
        epilog=__doc__,
    )

    params, options = parser.parse_args(args=args, show_diff_phil=False)
    dials.util.log.config(verbosity=options.verbose, logfile=params.output.log)

    experiments = flatten_experiments(params.input.experiments)
    reflections = flatten_reflections(params.input.reflections)
    if len(reflections) != 1:
        sys.exit("Exactly one reflection file needed.")

    bijvoet_table(experiments, reflections[0], params)


if __name__ == "__main__":
    run()
