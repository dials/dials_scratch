"""
Group reflections by symmetry equivalents and print in csv format

Usage: dials.python bijvoet_group scaled.expt scaled.refl
"""

import sys
from math import ceil
from dials.util.options import OptionParser, flatten_experiments, flatten_reflections
import libtbx.phil
import dials.util
from cctbx import miller

# No PHIL options for this script, yet
phil_scope = libtbx.phil.parse("")


def bijvoet_table(experiments, reflections, params):
    for experiment in experiments:
        rt = reflections.select_on_experiment_identifiers(
            [
                experiment.identifier,
            ]
        )
        exp_id = list(set(rt["id"]))
        assert len(exp_id) == 1
        symm = experiment.crystal.get_crystal_symmetry()

        ms = miller.set(symm, anomalous_flag=False, indices=rt["miller_index"])
        rt["miller_index_asu"] = ms.map_to_asu().indices()
        rt.sort("miller_index_asu")

        print(f"\nDataset: {experiment.imageset.get_template()}")
        print("ASU_HKL,HKL,image,intensity")
        for ref in rt.rows():
            line = [
                " ".join(str(e) for e in ref["miller_index_asu"]),
                " ".join(str(e) for e in ref["miller_index"]),
                str(ceil(ref["xyzcal.px"][2])),
                str(ref["intensity.scale.value"]),
            ]
            print(",".join(line))


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

    experiments = flatten_experiments(params.input.experiments)
    reflections = flatten_reflections(params.input.reflections)

    if len(reflections) != 1:
        sys.exit("Exactly one reflection file needed.")

    bijvoet_table(experiments, reflections[0], params)


if __name__ == "__main__":
    run()
