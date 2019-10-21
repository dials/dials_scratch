from __future__ import absolute_import, division, print_function

from cctbx import crystal, sgtbx, uctbx

symmetries = []

for unit_cell, space_group in [
    ((39.7413, 183.767, 140.649, 90, 90, 90), "C 2 2 21"),
    ((40.2365, 180.881, 143.335, 90, 90, 90), "C 2 2 21"),
    ((40.3094, 181.286, 143.971, 90, 90, 90), "C 2 2 21"),
    ((40.3785, 181.635, 144.537, 90, 90, 90), "C 2 2 21"),
    ((40.1506, 180.407, 142.743, 90, 90, 90), "C 2 2 21"),
    ((40.1721, 180.663, 142.75, 90, 90, 90), "C 2 2 21"),
    ((40.16, 142.899, 92.4167, 90, 102.48, 90), "P 1 2 1"),
    ((180.613, 40.1558, 142.737, 90, 90.0174, 90), "C 1 2 1"),
    ((40.1561, 180.55, 142.684, 90, 90, 90), "C 2 2 21"),
    ((180.138, 40.1064, 142.443, 90, 90.1121, 90), "C 1 2 1"),
    ((40.2263, 180.891, 143.374, 90, 89.9327, 90), "C 1 2 1"),
    ((40.166, 180.615, 142.794, 90, 90, 90), "C 2 2 21"),
]:
    symmetries.append(
        crystal.symmetry(
            unit_cell=uctbx.unit_cell(unit_cell),
            space_group_info=sgtbx.space_group_info(space_group),
        )
    )

for symmetry in symmetries:
    print(
        " ".join(["%s"] * 6) % symmetry.unit_cell().parameters()
        + " "
        + str(symmetry.space_group_info()).replace(" ", "")
    )
print()

for cs in symmetries:
    cs = cs.primitive_setting()
    print(cs.unit_cell(), cs.space_group_info())

print()

for cs in symmetries:
    cs = cs.minimum_cell()
    print(cs.unit_cell(), cs.space_group_info())

print()

for cs in symmetries:
    cs = cs.niggli_cell()
    print(cs.unit_cell(), cs.space_group_info())
