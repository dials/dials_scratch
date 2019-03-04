'''Comparison of execution times for two functions that do the same thing.

Result: dR_from_axis_and_angle is approximately 3 times faster than
matrix.rec.axis_and_angle_as_r3_derivative_wrt_angle.
'''
from __future__ import print_function

from scitbx import matrix
from scitbx.array_family import flex #required dependency
import random
from math import pi
from libtbx.test_utils import approx_equal
from libtbx.development.timers import Timer
from dials.algorithms.refinement.refinement_helpers import \
  dR_from_axis_and_angle as dR_from_axis_and_angle_cpp
from dials.algorithms.refinement.refinement_helpers import \
  dR_from_axis_and_angle_py

trials = []

for i in range(10000):
  # generate random axis and angle
  trials.append((matrix.col.random(3, -1, 1).normalize(),
                 random.uniform(0, 2*pi)))

t = Timer('dR_from_axis_and_angle in C++')
dR_0 = []
for trial in trials:
  dR_0.append(dR_from_axis_and_angle_cpp(trial[0], trial[1]))
del(t)

print()

t = Timer('dR_from_axis_and_angle in Python')
dR_1 = []
for trial in trials:
  dR_1.append(dR_from_axis_and_angle_py(trial[0], trial[1]))
del(t)

print()

t = Timer('axis_and_angle_as_r3_derivative_wrt_angle')
dR_2 = []
for trial in trials:
  dR_2.append(trial[0].axis_and_angle_as_r3_derivative_wrt_angle(trial[1]))
del (t)

print("check results are the same")
for (a, b) in zip(dR_1, dR_2):
  assert approx_equal(a, b)
print("OK")
