from __future__ import print_function

from cctbx import crystal
from cctbx.sgtbx.lattice_symmetry import metric_subgroups

ucs_1 = [
    (
        4.805202948916906,
        12.808064769657364,
        16.544899201125446,
        106.45808502003258,
        90.0065567098825,
        100.77735674275475,
    ),
    (
        4.808011343212577,
        12.821894835790472,
        16.557339561965573,
        106.48431244651402,
        90.0252848479048,
        100.77252933676507,
    ),
    (
        4.8096632137789985,
        12.815648858527567,
        16.55931712239122,
        106.48990701341536,
        90.01703141314147,
        100.80397887485773,
    ),
    (
        4.807294085194974,
        12.822386757910516,
        16.560411742466663,
        106.43185845358086,
        90.02067929544215,
        100.79522302759383,
    ),
]
sgs_1 = ["P1"] * 4

print("Test case 1:")
for uc, sg in zip(ucs_1, sgs_1):
    cs = crystal.symmetry(unit_cell=uc, space_group_symbol=sg)
    groups = metric_subgroups(
        cs, max_delta=5, enforce_max_delta_for_generated_two_folds=True
    )
    group = groups.result_groups[0]
    print(cs)
    print("Minimum cell:", cs.minimum_cell().unit_cell())
    print("Best cell:", group["best_subsym"].unit_cell())
    print("Minimum cell (via best):", group["best_subsym"].minimum_cell().unit_cell())
    print()


ucs_2 = [
    (39.7413, 183.767, 140.649, 90, 90, 90),
    (40.16, 142.899, 92.4167, 90, 102.48, 90),
    (180.613, 40.1558, 142.737, 90, 90.0174, 90),
]
sgs_2 = ["C 2 2 21", "P 1 2 1", "C 1 2 1"]

print("Test case 2:")
for uc, sg in zip(ucs_2, sgs_2):
    cs = crystal.symmetry(unit_cell=uc, space_group_symbol=sg)
    groups = metric_subgroups(
        cs, max_delta=5, enforce_max_delta_for_generated_two_folds=True
    )
    group = groups.result_groups[0]
    print(cs)
    print("Minimum cell:", cs.minimum_cell().unit_cell())
    print("Best cell:", group["best_subsym"].unit_cell())
    print("Minimum cell (via best):", group["best_subsym"].minimum_cell().unit_cell())
    print()
