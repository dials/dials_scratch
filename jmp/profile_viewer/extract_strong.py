from __future__ import print_function


def extract_strong(directory):
    from glob import glob
    import os.path
    import shutil

    from dials.array_family import flex

    reflections = flex.reflection_table()

    filenames = glob(os.path.join(directory, "shoeboxes_*_*.pickle"))
    assert len(filenames) > 0

    for filename in sorted(filenames):
        r = flex.reflection_table.from_pickle(filename)
        print("Read %d reflections from %s" % (len(r), filename))
        selection = r.get_flags(r.flags.used_in_refinement)
        r = r.select(selection)
        print("Selected %d strong reflections" % len(r))
        reflections.extend(r)

    print("Writing reflections to reflections.pickle")
    reflections.as_pickle("reflections.pickle")

    print("Copying experiments")
    shutil.copyfile(
        os.path.join(directory, "integrated_experiments.json"), "experiments.json"
    )


def generate_profile(sample=None, grid_size=15, grid_range=0.25):
    from dxtbx.model.experiment.experiment_list import ExperimentListFactory
    from dials.array_family import flex

    experiments = ExperimentListFactory.from_json_file("experiments.json")
    reflections = flex.reflection_table.from_pickle("reflections.pickle")

    import random

    if sample is not None:
        selection = flex.size_t(
            random.sample(flex.size_t(range(len(reflections))), sample)
        )
        reflections = reflections.select(selection)
        print("Selected %d random reflections" % len(reflections))

    beam = experiments[0].beam
    detector = experiments[0].detector
    goniometer = experiments[0].goniometer
    scan = experiments[0].scan
    crystal = experiments[0].crystal

    from scitbx import matrix

    s0 = matrix.col(beam.get_s0())
    m2 = matrix.col(goniometer.get_rotation_axis())
    Rot = matrix.sqr(goniometer.get_fixed_rotation())

    image_range = scan.get_array_range()
    if crystal.num_scan_points > 0:
        print("Using scan varying")
        assert crystal.num_scan_points == scan.get_num_images() + 1
        A = dict(
            (z, crystal.get_A_at_scan_point(i))
            for i, z in enumerate(range(*image_range))
        )
    else:
        print("Using scan static")
        A = dict((z, crystal.get_A()) for z in range(*image_range))

    def get_coord(x, y, z):
        A_matrix = A[int(z)]
        s1 = (
            matrix.col(detector[0].get_pixel_lab_coord((x, y))).normalize()
            * s0.length()
        )
        r0 = s1 - s0
        phi = scan.get_angle_from_array_index(z, deg=False)
        R = m2.axis_and_angle_as_r3_rotation_matrix(phi, deg=False)
        h = (R * Rot * A_matrix).inverse() * r0
        return h

    from dials.array_family import flex

    assert grid_range <= 0.5 and grid_range > 0
    assert grid_size > 0
    num = 2 * grid_size + 1
    grid = flex.double(flex.grid(num, num, num))
    min_x = -grid_range
    max_x = grid_range
    scale1 = num / (max_x - min_x)
    scale0 = -min_x * scale1

    from dials.algorithms.shoebox import MaskCode

    for idx, r in enumerate(reflections):
        hkl = r["miller_index"]
        sbox = r["shoebox"]
        data = sbox.data
        bgrd = sbox.background
        mask = sbox.mask
        H = get_coord(*r["xyzcal.px"])
        H = H[0] - hkl[0], H[1] - hkl[1], H[2] - hkl[2]
        print(idx, len(reflections), H)

        x, y, z = r["xyzcal.px"]
        coords = []
        x0, x1, y0, y1, z0, z1 = sbox.bbox

        data = data - bgrd
        min_data = flex.max(data) * 0.02

        for k in range(z1 - z0):
            for j in range(y1 - y0):
                for i in range(x1 - x0):
                    n = data[k, j, i]
                    if mask[k, j, i] == 5 and n > min_data:
                        for kk in range(5):
                            for jj in range(5):
                                for ii in range(5):
                                    kkk = z0 + k + kk / 5.0
                                    jjj = y0 + j + jj / 5.0
                                    iii = x0 + i + ii / 5.0
                                    h = get_coord(iii, jjj, kkk)
                                    d = n / (5.0 * 5.0 * 5.0)
                                    hf = h[0] - hkl[0]
                                    kf = h[1] - hkl[1]
                                    lf = h[2] - hkl[2]
                                    i0 = int(scale0 + scale1 * hf)
                                    j0 = int(scale0 + scale1 * kf)
                                    k0 = int(scale0 + scale1 * lf)
                                    if (
                                        i0 >= 0
                                        and i0 < num
                                        and j0 >= 0
                                        and j0 < num
                                        and k0 >= 0
                                        and k0 < num
                                    ):
                                        grid[k0, j0, i0] += d

    import cPickle as pickle

    print("Saving grid to grid.pickle")
    pickle.dump(grid, open("grid.pickle", "w"))


def show_profile(grid_range):

    from dials.array_family import flex

    import cPickle as pickle

    grid = pickle.load(open("grid.pickle"))
    min_x = -grid_range
    max_x = grid_range

    print(flex.max(grid), flex.min(grid))

    import numpy

    vmax = max(grid)
    grid = grid.as_numpy_array()

    from dials_scratch.jmp.viewer import show_image_stack_multi_view

    show_image_stack_multi_view(
        grid, vmax=vmax, extent=[min_x, max_x, min_x, max_x], axis_names=["L", "K", "H"]
    )


def crystal_parameters():
    from dxtbx.model.experiment.experiment_list import ExperimentListFactory
    from dials.array_family import flex
    from scitbx import matrix

    experiments = ExperimentListFactory.from_json_file("experiments.json")
    reflections = flex.reflection_table.from_pickle("reflections.pickle")

    beam = experiments[0].beam
    detector = experiments[0].detector
    goniometer = experiments[0].goniometer
    scan = experiments[0].scan
    crystal = experiments[0].crystal

    from dials.algorithms.refinement.parameterisation.crystal_parameters import (
        CrystalUnitCellParameterisation,
    )

    parameters = CrystalUnitCellParameterisation(crystal)

    print("")
    print("Parameters")
    for param in parameters.get_params():
        print(param.name, param.value)

    print("")
    print("db/dp")
    for db_dp in parameters.get_ds_dp():
        print(db_dp.round(5))

    print("")
    print("B matrix")
    B = crystal.get_B().round(7)
    print(B)

    print(B.transpose() * B)

    # from rstbx.symmetry.constraints.parameter_reduction import symmetrize_reduce_enlarge

    # temp = symmetrize_reduce_enlarge(crystal.get_space_group())
    # temp.set_orientation(crystal.get_B())
    # independent = temp.forward_independent_parameters()
    # print independent

    # new_mm = temp.enlarge(independent)
    # temp.Bconverter.validate_and_setG(new_mm)
    # B = temp.Bconverter.back_as_orientation()
    # print matrix.sqr(B.reciprocal_matrix()).round(5)

    # from cctbx import sgtbx
    # from rstbx.symmetry.constraints import AGconvert

    # constraints = sgtbx.tensor_rank_2_constraints(space_group=crystal.get_space_group(), reciprocal_space=True)

    # params = constraints.all_params(independent_params=tuple(p.value for p in
    #                                                       parameters.get_params()))

    # from cctbx.crystal_orientation import crystal_orientation

    # converter = AGconvert()
    # converter.forward(crystal_orientation(crystal.get_B(), False))
    # converter.validate_and_setG(params)
    # orientation = converter.back_as_orientation()

    # print orientation
    # print ""
    print("SUM(p * dp/dp)")
    MAT = [
        p.value * db_dp
        for p, db_dp in zip(parameters.get_params(), parameters.get_ds_dp())
    ]
    COV = sum(MAT[1:], MAT[0])
    print(2 * COV.round(7))
    print((2 * COV - crystal.get_B()).round(7))
