from scitbx import matrix
import sys


def read_expt(filename):
    from dials.util.options import flatten_datablocks
    from dials.util.phil import DataBlockConverters

    converter = DataBlockConverters(check_format=False)
    return flatten_datablocks([converter.from_string(filename)])


def write_expt(experiments, filename):
    from dxtbx.datablock import DataBlockDumper

    dump = DataBlockDumper(experiments)
    dump.as_file(filename)


expts = read_expt(sys.argv[1])
expt = expts[0]

sweeps = expt.extract_sweeps()

expt = sweeps[0]

gonio = expt.get_goniometer()

e1 = matrix.col((1, 0, 0))
e2 = matrix.col((0.6691306063588582, 0.7431448254773942, 0))
e3 = matrix.col((1, 0, 0))
z = matrix.col((0, 0, 1))

from dials.algorithms.refinement import rotation_decomposition

import math

a = math.pi / 3

Rplus = z.axis_and_angle_as_r3_rotation_matrix(a)
Rminus = z.axis_and_angle_as_r3_rotation_matrix(-a)

Splus = rotation_decomposition.solve_r3_rotation_for_angles_given_axes(
    Rplus, e1, e2, e3, return_both_solutions=True, deg=True
)
Sminus = rotation_decomposition.solve_r3_rotation_for_angles_given_axes(
    Rminus, e1, e2, e3, return_both_solutions=True, deg=True
)

s = Splus[0]
F = e2.axis_and_angle_as_r3_rotation_matrix(
    s[1]
) * e3.axis_and_angle_as_r3_rotation_matrix(s[2])

gonio.set_fixed_rotation(F.elems)
gonio.set_angles(s)

write_expt(expts, sys.argv[2])
