from cctbx import crystal, sgtbx, uctbx
import scitbx.matrix
from dxtbx.model import Crystal


def plot_cells(
    crystals, ax=None, facecolor="cyan", linewidth=1, edgecolor="r", alpha=0.1, **kwargs
):
    import numpy as np
    from mpl_toolkits.mplot3d import Axes3D
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
    import matplotlib.pyplot as plt

    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
    else:
        fig = None

    for cryst in crystals:
        A = scitbx.matrix.sqr(cryst.get_A()).inverse()
        a = scitbx.matrix.row(A[:3])
        b = scitbx.matrix.row(A[3:6])
        c = scitbx.matrix.row(A[6:])

        Z = np.array(
            [
                # front face clockwise starting from origin
                (0, 0, 0),
                a,
                a + b,
                b,
                # back face clockwise
                c,
                a + c,
                a + b + c,
                b + c,
            ]
        )

        # plot vertices
        ax.scatter3D(Z[:, 0], Z[:, 1], Z[:, 2])

        # list of sides' polygons of figure
        verts = [
            [Z[0], Z[1], Z[2], Z[3]],
            [Z[4], Z[5], Z[6], Z[7]],
            [Z[0], Z[1], Z[5], Z[4]],
            [Z[2], Z[3], Z[7], Z[6]],
            [Z[1], Z[2], Z[6], Z[5]],
            [Z[4], Z[7], Z[3], Z[0]],
        ]

        # plot sides
        pc = Poly3DCollection(verts, linewidth=linewidth, edgecolor=edgecolor, **kwargs)
        # https://github.com/matplotlib/matplotlib/issues/10237/
        pc.set_alpha(alpha)
        pc.set_facecolor(facecolor)
        ax.add_collection3d(pc)

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")

    if fig is not None:
        plt.show()
        plt.close()


def crystal_symmetry_as_dxtbx_crystal(crystal_symmetry):
    B = scitbx.matrix.sqr(
        crystal_symmetry.unit_cell().fractionalization_matrix()
    ).transpose()
    return Crystal(B, crystal_symmetry.space_group())


if __name__ == "__main__":
    uc_min = uctbx.unit_cell(
        (
            44.66208170999999,
            53.12629402999999,
            62.53397660267584,
            115.13670293012744,
            101.7265610491002,
            90.0,
        )
    )
    cs_min = crystal.symmetry(unit_cell=uc_min, space_group=sgtbx.space_group())
    from cctbx.sgtbx.lattice_symmetry import metric_subgroups

    subgroups = metric_subgroups(cs_min, max_delta=5)
    best_subgroup = subgroups.result_groups[0]
    print(f"best_sybsym {best_subgroup['best_subsym']}\n")
    print(f"ref_sybsym {best_subgroup['ref_subsym']}\n")

    cs_ref = crystal.symmetry(
        # unit_cell=(113.2236274, 53.12629403, 44.66208171, 90, 102.9736126, 90),
        unit_cell=(112.90, 53.14, 44.39, 90.00, 103.04, 90.00),
        space_group_symbol="C 1 2/m 1",
    )

    cs_best = cs_ref.best_cell()
    cb_ref_best = cs_ref.change_of_basis_op_to_best_cell()
    cb_ref_primitive = cs_ref.change_of_basis_op_to_primitive_setting()
    cb_best_primitive = cs_best.change_of_basis_op_to_primitive_setting()
    cb_best_ref = cs_best.change_of_basis_op_to_reference_setting()

    cryst = crystal_symmetry_as_dxtbx_crystal(cs_ref)
    cryst_best = cryst.change_basis(cb_ref_best)
    cryst_ref_primitive = cryst.change_basis(cb_ref_primitive)
    cryst_best_primitive = cryst_best.change_basis(cb_best_primitive)
    cryst_best_ref = cryst_best.change_basis(cb_best_ref)

    print(f"Reference symmetry:\n{cryst}\n")
    print(f"Best cell symmetry:\n{cryst_best}\n")
    print(f"Reference -> primitive symmetry:\n{cryst_ref_primitive}\n")
    print(f"Best -> primitive symmetry:\n{cryst_best_primitive}\n")
    print(f"Best -> reference symmetry:\n{cryst_best_ref}\n")

    from matplotlib import pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d", proj_type="ortho")
    plot_cells([cryst], ax=ax, facecolor="green")
    plot_cells([cryst_best], ax=ax, facecolor="blue")
    plot_cells([cryst_ref_primitive], ax=ax, facecolor="cyan")
    plot_cells([cryst_best_primitive], ax=ax, facecolor="red")
    plot_cells([cryst_best_ref], ax=ax, facecolor="purple")

    plt.show()
    plt.close()
