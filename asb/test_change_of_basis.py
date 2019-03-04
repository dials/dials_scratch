from __future__ import division
from __future__ import print_function
from cctbx.crystal_orientation import crystal_orientation, basis_type
from scitbx.matrix import sqr, col
from cctbx.sgtbx import change_of_basis_op, space_group_info
from cctbx.uctbx import unit_cell
from dxtbx.model import Crystal
from libtbx.test_utils import approx_equal

"""
Jiffy script to show that crystal_orientation.change_basis and
Crystal.change_basis work properly! Just need to add a transpose
call when using change_of_basis_op with the crystal_orientation
object.
"""


def print_matrix(m):
    print("%7.4f %7.4f %7.4f" % m[0:3])
    print("%7.4f %7.4f %7.4f" % m[3:6])
    print("%7.4f %7.4f %7.4f" % m[6:9])
    print()


# Cubic
# uc = unit_cell((100,100,100,90,90,90))
# sg = space_group_info("P 21 3").group()

# Thermolysin
# uc = unit_cell((93,93,131,90,90,120.0))
# sg = space_group_info("P 61 2 2").group()

# Lysozyme
uc = unit_cell((77, 77, 37, 90, 90, 90))
sg = space_group_info("P 43 21 2").group()

print(uc)
print(sg.info())

direct_matrix = sqr(uc.orthogonalization_matrix()).transpose()
crystal = Crystal(
    direct_matrix * col((1, 0, 0)),
    direct_matrix * col((0, 1, 0)),
    direct_matrix * col((0, 0, 1)),
    sg,
)

co = crystal_orientation(crystal.get_A(), basis_type.reciprocal)

print(crystal)
print(co)

ok_ops = []
bad_ops = []


def test_op(op):
    print("=" * 80)
    print("COB:", op)
    print("Crystal A")
    print_matrix(crystal.get_A())
    print("cctbx   A")
    print_matrix(co.reciprocal_matrix())
    dxtbx_a = sqr(crystal.change_basis(op).get_A())
    cctbx_a = sqr(
        co.change_basis(
            sqr(op.c().as_double_array()[0:9]).transpose()
        ).reciprocal_matrix()
    )
    print("Crystal A COB")
    print_matrix(dxtbx_a)
    print("cctbx   A COB")
    print_matrix(cctbx_a)

    good_op = approx_equal(dxtbx_a.elems, cctbx_a.elems, out=None)
    print("A matrices approx equal:", good_op)
    if good_op:
        ok_ops.append(op)
    else:
        bad_ops.append(op)


print("All possible ops")
for rot in crystal.get_space_group():
    op = change_of_basis_op(rot)
    test_op(op)

print("Ops that passed")
for op in ok_ops:
    print(op)
print("Ops that failed")
for op in bad_ops:
    print(op)
