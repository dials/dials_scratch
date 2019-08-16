from __future__ import print_function
from scitbx import matrix


def read_expt(filename):
    from dials.phil import ExperimentListConverters
    from dials.util.options import flatten_experiments

    converter = ExperimentListConverters(check_format=False)
    return flatten_experiments([converter.from_string(filename)])


def write_expt(experiments, filename):
    experiments.as_file(filename)


import sys

expts = read_expt(sys.argv[1])
expt = expts[0]

scan = expt.scan
gonio = expt.goniometer

assert len(expt.detector) == 1

# crud trying to work out correct min two theta

p = expt.detector[0]

dims = p.get_image_size()
pixel = p.get_pixel_size()

size = tuple([_d * _p for _d, _p in zip(dims, pixel)])
fast = matrix.col(p.get_fast_axis())
slow = matrix.col(p.get_slow_axis())
origin = matrix.col(p.get_origin())

s0n = matrix.col(expt.beam.get_s0()).normalize()
panel, xy = expt.detector.get_ray_intersection(s0n)
zero = origin + xy[0] * fast + xy[1] * slow

resolution = expt.detector.get_max_inscribed_resolution(expt.beam.get_s0())
import math

# this is wrong https://github.com/dials/dials/issues/348
# however right enough for this...

theta = math.asin(expt.beam.get_wavelength() / (2 * resolution))
print("Using two-theta: %.3f" % (2 * theta * 180.0 / math.pi))

epochs = scan.get_epochs()
exposure_times = scan.get_exposure_times()
image_range = scan.get_image_range()
oscillation = scan.get_oscillation()

current = 1 + image_range[1] - image_range[0]
turn = int(round(360.0 / oscillation[1]))
extra = turn - current

for j in range(extra):
    epochs.append(0.0)
    exposure_times.append(0.0)

image_range = image_range[0], image_range[1] + extra

scan.set_image_range(image_range)
scan.set_epochs(epochs)
scan.set_exposure_times(exposure_times)

write_expt(expts, sys.argv[2])

# now for amusement try decomposing rotation of 90 degrees about beam to
# measure blind region - computer says no if mini kappa :(

e1 = matrix.col((1, 0, 0))
e2 = matrix.col((0.914, 0.279, -0.297))
e3 = matrix.col((1, 0, 0))

R = matrix.sqr((0, 1, 0, -1, 0, 0, 0, 0, 1))

from dials.algorithms.refinement import rotation_decomposition

solutions = rotation_decomposition.solve_r3_rotation_for_angles_given_axes(
    R, e1, e2, e3, return_both_solutions=True, deg=True
)

assert solutions is None

# now try getting a rotation of two-theta about the beam - this should (i)
# be possible? and (ii) move the blind region into somewhere we can actually
# record...

R_tt = s0n.axis_and_angle_as_r3_rotation_matrix(2 * theta)

s = rotation_decomposition.solve_r3_rotation_for_angles_given_axes(
    R_tt, e1, e2, e3, return_both_solutions=False, deg=False
)

# use solution

print("Using angles: %.3f %.3f" % (180 * s[1] / math.pi, 180 * s[2] / math.pi))

F = e2.axis_and_angle_as_r3_rotation_matrix(
    s[1]
) * e3.axis_and_angle_as_r3_rotation_matrix(s[2])

gonio.set_fixed_rotation(F.elems)
write_expt(expts, sys.argv[3])
