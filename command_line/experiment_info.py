from __future__ import division, print_function


def experiment_info(expt):
    from scitbx import matrix

    detector = expt.detector
    crystal = expt.crystal
    assert len(detector) == 1
    panel = detector[0]
    fast = matrix.col(panel.get_fast_axis()).normalize()
    slow = matrix.col(panel.get_slow_axis()).normalize()
    origin = matrix.col(panel.get_origin())
    normal = fast.cross(slow)
    distance = abs(origin.dot(normal))
    unit_cell = crystal.get_unit_cell().parameters()
    print(
        distance,
        unit_cell[0],
        unit_cell[1],
        unit_cell[2],
        unit_cell[3],
        unit_cell[4],
        unit_cell[5],
    )


if __name__ == "__main__":

    from dials.util.options import OptionParser

    parser = OptionParser(check_format=False, read_experiments=True)
    from dials.util.options import flatten_experiments

    params, options = parser.parse_args(show_diff_phil=True)
    experiments = flatten_experiments(params.input.experiments)
    for expt in experiments:
        experiment_info(expt)
