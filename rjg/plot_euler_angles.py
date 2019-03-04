from __future__ import absolute_import, division
from cctbx.array_family import flex

import iotbx.phil

help_message = """
"""


phil_scope = iotbx.phil.parse(
    """
"""
)


def unwrap_angles(angles, deg=False):
    import math

    pi = math.pi
    two_pi = math.pi * 2
    if deg:
        pi = 180
        two_pi = 360
    prev_angle = angles[0]
    unwrapped_angles = flex.double([angles[0]])
    for i, angle in enumerate(angles):
        if i == 0:
            continue
        while (prev_angle - angle) > pi:
            angle += two_pi
        while (angle - prev_angle) > pi:
            angle -= two_pi
        unwrapped_angles.append(angle)
        prev_angle = angle
    return unwrapped_angles


def run(args):

    from dials.util.options import OptionParser
    from dials.util.options import flatten_experiments
    import libtbx.load_env

    usage = "%s [options] experiments.json" % libtbx.env.dispatcher_name

    parser = OptionParser(
        usage=usage,
        phil=phil_scope,
        read_experiments=True,
        check_format=False,
        epilog=help_message,
    )

    params, options = parser.parse_args(show_diff_phil=True)
    experiments = flatten_experiments(params.input.experiments)
    if len(experiments) <= 1:
        parser.print_help()
        return

    import scitbx.matrix

    euler_angles = flex.vec3_double()
    expt_ids = flex.int()
    for i, experiment in enumerate(experiments):
        crystal = experiment.crystal
        U = scitbx.matrix.sqr(crystal.get_U())
        euler = U.r3_rotation_matrix_as_x_y_z_angles(deg=True)
        euler_angles.append(euler)
        expt_ids.append(i)

    x, y, z = [unwrap_angles(angles, deg=True) for angles in euler_angles.parts()]

    colour_map = "jet"
    import matplotlib
    from matplotlib import pyplot
    from mpl_toolkits.mplot3d import Axes3D

    cmap = matplotlib.cm.get_cmap(colour_map)
    colours = expt_ids.as_double()
    colours -= flex.min(colours)
    colours /= flex.max(colours)
    colours = [cmap(c) for c in colours]
    fig = pyplot.figure()

    ax = pyplot.scatter(x, y, c=colours)
    pyplot.xlabel("X")
    pyplot.ylabel("Y")
    pyplot.axes().set_aspect("equal")
    pyplot.savefig("xy.png")
    pyplot.clf()

    ax = pyplot.scatter(y, z, c=colours)
    pyplot.xlabel("Y")
    pyplot.ylabel("Z")
    pyplot.axes().set_aspect("equal")
    pyplot.savefig("yz.png")
    pyplot.clf()

    ax = pyplot.scatter(z, x, c=colours)
    pyplot.xlabel("Z")
    pyplot.ylabel("X")
    pyplot.axes().set_aspect("equal")
    pyplot.savefig("zx.png")
    pyplot.clf()

    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(x, y, z, c=colours)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_aspect("equal")
    pyplot.show()


if __name__ == "__main__":
    import sys

    run(sys.argv[1:])
