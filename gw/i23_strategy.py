from scitbx import matrix

e1 = matrix.col((1, 0, 0))
e2 = matrix.col((0.6691306063588582, 0.7431448254773942, 0))
e3 = matrix.col((1, 0, 0))
z = matrix.col((0, 0, 1))

from dials.algorithms.refinement import rotation_decomposition

Rplus = z.axis_and_angle_as_r3_rotation_matrix(1.1344640137963142)
Rminus = z.axis_and_angle_as_r3_rotation_matrix(-1.1344640137963142)

Splus = rotation_decomposition.solve_r3_rotation_for_angles_given_axes(
  Rplus, e1, e2, e3, return_both_solutions=True, deg=True)
Sminus = rotation_decomposition.solve_r3_rotation_for_angles_given_axes(
  Rminus, e1, e2, e3, return_both_solutions=True, deg=True)

for s in Splus:
  print '%7.2f %7.2f %7.2f' % s
for s in Sminus:
  print '%7.2f %7.2f %7.2f' % s
