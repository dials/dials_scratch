from __future__ import print_function
from dials.array_family import flex
from dxtbx.model.experiment_list import ExperimentListFactory
from sklearn.neighbors import NearestNeighbors
from random import shuffle
from cctbx import miller
from matplotlib import pylab
from math import sqrt
from scipy.ndimage import convolve
from matplotlib import cm
import pptk


def compute_local_mean_over_variance(reflections, kernel_size):

    # Get the miller indices, intensities and variances
    H, K, L = zip(*list(reflections["miller_index"]))
    I = reflections["intensity"]
    V = reflections["variance"]

    # Get the range of miller indices
    min_H, max_H = min(H), max(H)
    min_K, max_K = min(K), max(K)
    min_L, max_L = min(L), max(L)

    n_H = max_H - min_H + 1
    n_K = max_K - min_K + 1
    n_L = max_L - min_L + 1

    # Construct 3D arrays of items
    N_array = flex.int(flex.grid(n_L, n_K, n_H))
    I_array = flex.double(N_array.accessor())
    V_array = flex.double(N_array.accessor())
    H_array = flex.double(N_array.accessor())
    K_array = flex.double(N_array.accessor())
    L_array = flex.double(N_array.accessor())
    for h, k, l, intensity, variance in zip(H, K, L, I, V):
        x = h - min_H
        y = k - min_K
        z = l - min_L
        assert x >= 0 and y >= 0 and z >= 0
        N_array[z, y, x] += 1
        I_array[z, y, x] += intensity
        V_array[z, y, x] += variance
        H_array[z, y, x] = h
        K_array[z, y, x] = k
        L_array[z, y, x] = l

    # Compute the Intensity / Variance
    # mask = N_array > 0
    mask = V_array > 0
    selection = mask.as_1d()
    I_sub = I_array.as_1d().select(selection)
    V_sub = V_array.as_1d().select(selection)
    R_sub = I_sub / V_sub
    data = flex.double(len(mask))
    data.set_selected(selection, R_sub)
    data.reshape(mask.accessor())

    # Create an integer count array
    mask = mask.as_1d().as_int().as_double()
    mask.reshape(data.accessor())

    # Convolve the count array and intensity/variance array
    kernel = flex.double(flex.grid(kernel_size, kernel_size, kernel_size), 1.0)
    convolved_data = convolve(data.as_numpy_array(), kernel.as_numpy_array())
    convolved_mask = convolve(mask.as_numpy_array(), kernel.as_numpy_array())

    # Select data points
    convolved_data = flex.double(convolved_data)
    convolved_mask = flex.double(convolved_mask)
    selection = (convolved_mask > 0.5 * kernel_size ** 3).as_1d()
    R_sub = convolved_data.as_1d().select(selection)
    N_sub = convolved_mask.as_1d().select(selection)
    R_sub /= N_sub
    H_sub = H_array.as_1d().select(selection)
    K_sub = K_array.as_1d().select(selection)
    L_sub = L_array.as_1d().select(selection)
    return H_sub, K_sub, L_sub, R_sub


def compute_local_mean_over_variance2(reflections, kernel_size):

    # Get the miller indices, intensities and variances
    H, K, L = zip(*list(reflections["miller_index"]))
    I = reflections["intensity"]

    # Get the range of miller indices
    min_H, max_H = min(H), max(H)
    min_K, max_K = min(K), max(K)
    min_L, max_L = min(L), max(L)

    n_H = max_H - min_H + 1
    n_K = max_K - min_K + 1
    n_L = max_L - min_L + 1

    # Construct 3D arrays of items
    N_array = flex.int(flex.grid(n_L, n_K, n_H))
    I_array = flex.double(N_array.accessor())
    S_array = flex.double(N_array.accessor())
    H_array = flex.double(N_array.accessor())
    K_array = flex.double(N_array.accessor())
    L_array = flex.double(N_array.accessor())
    for h, k, l, intensity in zip(H, K, L, I):
        x = h - min_H
        y = k - min_K
        z = l - min_L
        assert x >= 0 and y >= 0 and z >= 0
        N_array[z, y, x] += 1
        I_array[z, y, x] += intensity
        S_array[z, y, x] += intensity ** 2
        H_array[z, y, x] = h
        K_array[z, y, x] = k
        L_array[z, y, x] = l

    # Convolve the count array and intensity/variance array
    kernel = flex.double(flex.grid(kernel_size, kernel_size, kernel_size), 1.0)
    convolved_I = convolve(I_array.as_numpy_array(), kernel.as_numpy_array())
    convolved_S = convolve(S_array.as_numpy_array(), kernel.as_numpy_array())
    convolved_N = convolve(
        N_array.as_double().as_numpy_array(), kernel.as_numpy_array()
    )

    # Select data points
    convolved_I = flex.double(convolved_I)
    convolved_S = flex.double(convolved_S)
    convolved_N = flex.double(convolved_N)
    print(max(convolved_N))
    selection = (convolved_N >= 0.5 * kernel_size ** 3).as_1d()
    I_sub = convolved_I.as_1d().select(selection)
    S_sub = convolved_S.as_1d().select(selection)
    N_sub = convolved_N.as_1d().select(selection)
    V_sub = S_sub / N_sub - (I_sub / N_sub) ** 2
    R_sub = (I_sub / N_sub) / V_sub
    H_sub = H_array.as_1d().select(selection)
    K_sub = K_array.as_1d().select(selection)
    L_sub = L_array.as_1d().select(selection)
    return H_sub, K_sub, L_sub, R_sub


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

    # Get the miller indices, intensities and variances
    H, K, L = zip(*list(reflections["asu_miller_index"]))
    I = reflections["intensity"]

    # Get the range of miller indices
    min_H, max_H = min(H), max(H)
    min_K, max_K = min(K), max(K)
    min_L, max_L = min(L), max(L)

    n_H = max_H - min_H + 1
    n_K = max_K - min_K + 1
    n_L = max_L - min_L + 1

    num_array = flex.int(flex.grid(n_L, n_K, n_H))
    sum_array = flex.double(num_array.accessor())
    sum_sq_array = flex.double(num_array.accessor())
    H_array = flex.double(num_array.accessor())
    K_array = flex.double(num_array.accessor())
    L_array = flex.double(num_array.accessor())
    for h, k, l, intensity in zip(H, K, L, I):
        x = h - min_H
        y = k - min_K
        z = l - min_L
        assert x >= 0 and y >= 0 and z >= 0
        num_array[z, y, x] += 1
        sum_array[z, y, x] += intensity
        sum_sq_array[z, y, x] += intensity ** 2
        H_array[z, y, x] = h
        K_array[z, y, x] = k
        L_array[z, y, x] = l

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
    mask_out = flex.bool(mask.accessor())
    for k, j, i in indices:
        # print ii+min_H, jj+min_K, kk+min_L
        # for h, k, l in indices:
        #       ii = h - min_H
        #       jj = k - min_K
        #       kk = l - min_L
        # for k in range(mask.all()[0]):
        #   print k
        #   for j in range(mask.all()[1]):
        #     for i in range(mask.all()[2]):

        k0 = k - kernel_size
        k1 = k + kernel_size + 1
        j0 = j - kernel_size
        j1 = j + kernel_size + 1
        i0 = i - kernel_size
        i1 = i + kernel_size + 1
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

        if len(mean_list) > 0.25 * (2 * kernel_size + 1) ** 3:
            cchalf = compute_cchalf(mean_list, var_list)
            # Imean = sum(mean_list) / len(mean_list)
            # Ivar = sum((I-Imean)**2 for I in mean_list) / len(mean_list)
            # cchalf = Imean / sqrt(Ivar)
            cchalf_array[k, j, i] = cchalf
            mask_out[k, j, i] = True
        else:
            mask_out[k, j, i] = False

    selection = (mask_out == True).as_1d()
    H_sub = H_array.as_1d().select(selection)
    K_sub = K_array.as_1d().select(selection)
    L_sub = L_array.as_1d().select(selection)
    R_sub = cchalf_array.as_1d().select(selection)

    return H_sub, K_sub, L_sub, R_sub


if __name__ == "__main__":

    import sys

    # Get the filename
    reflections_filename = sys.argv[1]
    experiments_filename = sys.argv[2]

    # Set the number of neighbours
    if len(sys.argv) > 3:
        kernel_size = int(sys.argv[3])
    else:
        kernel_size = 5

    # Read the reflections
    reflections = flex.reflection_table.from_pickle(reflections_filename)
    experiments = ExperimentListFactory.from_json_file(
        experiments_filename, check_format=False
    )

    # Get only those reflections that have been indexed
    if "intensity.scaled.value" in reflections:
        selection = reflections.get_flags(reflections.flags.scaled)
        reflections["intensity"] = (
            reflections["intensity.scale.value"] / reflections["inverse_scale_factor"]
        )
        reflections["variance"] = (
            reflections["intensity.scale.variance"]
            / reflections["inverse_scale_factor"]
        )
    else:
        selection = reflections.get_flags(reflections.flags.integrated)
        reflections["intensity"] = reflections["intensity.sum.value"]
        reflections["variance"] = reflections["intensity.sum.variance"]
    reflections = reflections.select(selection)

    # Get the ASU miller index
    A = experiments[0].crystal.get_A()
    cs = experiments[0].crystal.get_crystal_symmetry()
    ms = miller.set(cs, reflections["miller_index"], anomalous_flag=False)
    ms_asu = ms.map_to_asu()
    reflections["asu_miller_index"] = ms_asu.indices()

    # Compute the local mean / variance
    # H, K, L, result = compute_local_mean_over_variance(reflections, kernel_size)
    H, K, L, result = compute_local_cchalf(reflections, kernel_size)

    # Get the coordinats of points to show
    X, Y, Z = (flex.mat3_double(len(H), A) * flex.vec3_double(H, K, L)).parts()

    # Get the colours
    colors = cm.plasma.colors
    print(len(result), len(H), len(K), len(L))
    print(min(result))
    print(max(result))
    min_r = 0
    max_r = 1.0  # max(result)
    # max_r = max(result)
    B = 255.0 / (max_r - min_r)
    A = -B * min_r
    rgba = []
    for r in result:
        i = A + B * r
        if i < 0:
            i = 0
        if i > 255:
            i = 255
        # alpha = i/255.0
        alpha = 0.1 + 0.9 * (1.0 - i / 255.0)
        rgba.append(colors[int(i)] + [alpha])

    print(min(X), max(X))
    print(min(Y), max(Y))
    print(min(Z), max(Z))

    # Display the points
    xyz = list(zip(X, Y, Z))
    v = pptk.viewer(xyz, rgba)
    v.set(show_grid=False)
    v.set(bg_color=[0, 0, 0, 0])
    v.set(point_size=0.0025)
    v.set(lookat=[0, 0, 0])
