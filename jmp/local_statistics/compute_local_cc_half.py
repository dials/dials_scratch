from __future__ import print_function


from dials.array_family import flex
from dxtbx.model.experiment_list import ExperimentListFactory
from sklearn.neighbors import NearestNeighbors
from random import shuffle
from cctbx import miller
from matplotlib import pylab
from math import sqrt


def compute_cchalf(mean, var):
    """
    Compute the CC 1/2 using the formular from Assmann, Brehm and Diederichs 2016

    :param mean: The list of mean intensities
    :param var: The list of variances on the half set of mean intensities
    :returns: The CC 1/2

    """
    assert len(mean) == len(var)
    n = len(mean)
    mean_of_means = sum(mean) / n
    sigma_e = sum(var) / n
    sigma_y = sum([(m - mean_of_means) ** 2 for m in mean]) / (n - 1)
    cchalf = (sigma_y - sigma_e) / (sigma_y + sigma_e)
    return cchalf


def compute_local_cchalf(reflections, kernel_size):

    # indices = list(set(reflections['miller_index']))

    H, K, L = zip(*list(reflections["miller_index"]))
    I = reflections["intensity.sum.value"]

    min_H, max_H = min(H), max(H)
    min_K, max_K = min(K), max(K)
    min_L, max_L = min(L), max(L)

    n_H = max_H - min_H + 1
    n_K = max_K - min_K + 1
    n_L = max_L - min_L + 1

    num_array = flex.int(flex.grid(n_L, n_K, n_H))
    sum_array = flex.double(num_array.accessor())
    sum_sq_array = flex.double(num_array.accessor())
    for h, k, l, intensity in zip(H, K, L, I):
        x = h - min_H
        y = k - min_K
        z = l - min_L
        assert x >= 0 and y >= 0 and z >= 0
        num_array[z, y, x] += 1
        sum_array[z, y, x] += intensity
        sum_sq_array[z, y, x] += intensity ** 2

    indices = []
    for k in range(num_array.all()[0]):
        for j in range(num_array.all()[1]):
            for i in range(num_array.all()[2]):
                if num_array[k, j, i] > 0:
                    indices.append((k, j, i))
                # else:
                #   n = 0
                #   if k > 0 and num_array[k-1,j,i] > 0:
                #     n += 1
                #   if j > 0 and num_array[k,j-1,i] > 0:
                #     n += 1
                #   if i > 0 and num_array[k,j,i-1] > 0:
                #     n += 1
                #   if k < num_array.all()[0]-1 and num_array[k+1,j,i] > 0:
                #     n += 1
                #   if j < num_array.all()[1]-1 and num_array[k,j+1,i] > 0:
                #     n += 1
                #   if i < num_array.all()[2]-1 and num_array[k,j,i+1] > 0:
                #     n += 1
                #   if n > 0:
                #     indices.append((k,j,i))

    mask = flex.bool(num_array.accessor())
    mean_array = flex.double(mask.accessor())
    var_array = flex.double(mask.accessor())
    for k in range(mask.all()[0]):
        for j in range(mask.all()[1]):
            for i in range(mask.all()[2]):

                n = num_array[k, j, i]
                sum_x = sum_array[k, j, i]
                sum_x2 = sum_sq_array[k, j, i]

                if n > 1:
                    mean = sum_x / n
                    var = (sum_x2 - (sum_x) ** 2 / n) / (n - 1)
                    var = var / n

                    mean_array[k, j, i] = mean
                    var_array[k, j, i] = var
                    mask[k, j, i] = True

    cchalf_array = flex.double(mask.accessor())
    for kk, jj, ii in indices:
        # print ii+min_H, jj+min_K, kk+min_L
        # for h, k, l in indices:
        #       ii = h - min_H
        #       jj = k - min_K
        #       kk = l - min_L
        # for k in range(mask.all()[0]):
        #   print k
        #   for j in range(mask.all()[1]):
        #     for i in range(mask.all()[2]):

        k0 = kk - kernel_size
        k1 = kk + kernel_size + 1
        j0 = jj - kernel_size
        j1 = jj + kernel_size + 1
        i0 = ii - kernel_size
        i1 = ii + kernel_size + 1
        if k0 < 0:
            k0 = 0
        if j0 < 0:
            j0 = 0
        if i0 < 0:
            i0 = 0
        if k1 > mask.all()[0]:
            k1 = mask.all()[0]
        if j1 > mask.all()[1]:
            j1 = mask.all()[1]
        if i1 > mask.all()[2]:
            i1 = mask.all()[2]

        mean_list = []
        var_list = []
        for kk in range(k0, k1):
            for jj in range(j0, j1):
                for ii in range(i0, i1):
                    if mask[kk, jj, ii]:
                        mean_list.append(mean_array[kk, jj, ii])
                        var_list.append(var_array[kk, jj, ii])

        if len(mean_list) > 1:
            # cchalf = compute_cchalf(mean_list, var_list)
            Imean = sum(mean_list) / len(mean_list)
            Ivar = sum((I - Imean) ** 2 for I in mean_list) / len(mean_list)
            cchalf = Imean / sqrt(Ivar)
            cchalf_array[kk, jj, ii] = cchalf

    return mask, cchalf_array, (min_H, max_H, min_K, max_K, min_L, max_L)


def display_local_cchalf(mask, cchalf, hrange):

    h = -hrange[0]  # cchalf.all()[2] // 2
    k = -hrange[2]  # cchalf.all()[1] // 2
    l = -hrange[4]  # cchalf.all()[0] // 2

    print(hrange)

    fig, (ax1, ax2, ax3) = pylab.subplots(ncols=3)
    ax1.imshow(
        cchalf.as_numpy_array()[l, :, :],
        vmin=-1,
        vmax=1,
        extent=(hrange[0], hrange[1], hrange[2], hrange[3]),
    )
    ax1.set_xlabel("H")
    ax1.set_ylabel("K")

    ax2.imshow(
        cchalf.as_numpy_array()[:, k, :],
        vmin=-1,
        vmax=1,
        extent=(hrange[0], hrange[1], hrange[4], hrange[5]),
    )
    ax2.set_xlabel("H")
    ax2.set_ylabel("L")

    ax3.imshow(
        cchalf.as_numpy_array()[:, :, h],
        vmin=-1,
        vmax=1,
        extent=(hrange[2], hrange[3], hrange[4], hrange[5]),
    )
    ax3.set_xlabel("K")
    ax3.set_ylabel("L")

    # fig.colorbar(im)
    pylab.show()


if __name__ == "__main__":

    import sys

    # Get the filename
    reflections_filename = sys.argv[1]
    experiments_filename = sys.argv[2]

    # Set the number of neighbours
    if len(sys.argv) > 3:
        kernel_size = int(sys.argv[3])
    else:
        kernel_size = 1

    # Read the reflections
    reflections = flex.reflection_table.from_pickle(reflections_filename)
    experiments = ExperimentListFactory.from_json_file(
        experiments_filename, check_format=False
    )

    selection = reflections.get_flags(reflections.flags.integrated)
    reflections = reflections.select(selection)

    cs = experiments[0].crystal.get_crystal_symmetry()
    ms = miller.set(cs, reflections["miller_index"])
    ms_asu = ms.map_to_asu()

    reflections["miller_index"] = ms_asu.indices()

    # Compute the local cchalf
    mask, cchalf, hrange = compute_local_cchalf(reflections, kernel_size)

    print(min(cchalf))
    print(max(cchalf))

    from dials_scratch.jmp.viewer import show_image_stack_multi_view

    show_image_stack_multi_view(
        cchalf.as_numpy_array(),
        interpolation="none",
        # vmin = min(cchalf),
        # vmax = max(cchalf),
        vmin=-1,
        vmax=1,
        axis_names=["H", "K", "L"],
    )

    # display_local_cchalf(mask, cchalf, hrange)
