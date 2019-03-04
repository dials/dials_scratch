from __future__ import print_function
from dials_scratch.lbl_feb_2017.apple import Apple
import sys

apple = Apple(sys.argv[1], sys.argv[2])
for image in sys.argv[3:]:
    apple.load(image)
    spots = apple.find_spots()
    indexed = apple.index(spots)
    apple.refine()
    U = apple.crystal.get_U()
    from scitbx import matrix

    rx, ry, rz = U.r3_rotation_matrix_as_x_y_z_angles()
    print("%.5f %.5f %.5f" % (rx, ry, rz), indexed.size())
