# make dx, dy distortion maps from two detector models, a moving and a reference
# from hierarchy_level=0 refinement and =1 respectively - at the moment it
# assumes P6M detector with exactly 60 panels
from __future__ import absolute_import, division, print_function

from scitbx.array_family import flex
from scitbx import matrix
from libtbx import phil

phil_scope = phil.parse(
    """
map_prefix = ''
  .type = str
  .help = 'Prefix for output dx, dy pickle files'
"""
)


def offset_dx_dy_p6m(detector0, detector1):
    """Compute pixel-remap function between detector1 and detector0, assuming
  panels are arranged [0][1][2][3][4] in the fast direction first etc. This
  also assumes that the rotation of pixels is locally small, but on the scale
  of a panel (80 mm ish) the rotation matters. N.B. Also assumes that the
  panel gap in the map is 7 pixels in fast direction and 17 in slow direction.

  Offsets are with respect to the fast and slow directions on the original
  coordinate system.
  """

    assert len(detector0) == 60
    assert len(detector1) == 60

    dx = flex.double(flex.grid(2527, 2463), 0.0)
    dy = flex.double(flex.grid(2527, 2463), 0.0)

    for k, (p0, p1) in enumerate(zip(detector0, detector1)):

        delta_x = 487 + 7
        delta_y = 195 + 17

        panel_x = k % 5
        panel_y = k // 5

        nx0 = panel_x * delta_x
        ny0 = panel_y * delta_y

        o0 = matrix.col(p0.get_origin())
        f0n = matrix.col(p0.get_fast_axis())
        s0n = matrix.col(p0.get_slow_axis())
        f0 = matrix.col(p0.get_fast_axis()) * 0.172
        s0 = matrix.col(p0.get_slow_axis()) * 0.172

        o1 = matrix.col(p1.get_origin())
        f1 = matrix.col(p1.get_fast_axis()) * 0.172
        s1 = matrix.col(p1.get_slow_axis()) * 0.172

        for j in range(195):
            for i in range(487):
                d = (o1 + f1 * i + s1 * j) - (o0 + f0 * i + s0 * j)
                dx[(ny0 + j, nx0 + i)] = -d.dot(f0n) / 0.172
                dy[(ny0 + j, nx0 + i)] = -d.dot(s0n) / 0.172

    import cPickle as pickle

    with open("dx.pickle", "w") as f:
        pickle.dump(dx, f)

    with open("dy.pickle", "w") as f:
        pickle.dump(dy, f)


def load_experiment(experiment_file):
    from dxtbx.model.experiment_list import ExperimentListFactory

    return ExperimentListFactory.from_json_file(experiment_file, check_format=False)


if __name__ == "__main__":
    import sys

    experiment0 = load_experiment(sys.argv[1])[0]
    experiment1 = load_experiment(sys.argv[2])[0]

    offset_dx_dy_p6m(experiment0.detector, experiment1.detector)
