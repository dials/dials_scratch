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
import numpy as np
import logging
import os
import os.path

logger = logging.getLogger(__name__)

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
    '''
    Compute the local CC 1/2

    '''
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


def plot_local_cc_half(reflections, A, kernel_size=5, record=False, directory="output", point_size=0.0025):
    '''
    Plot the local CC 1/2
    '''

    # Compute the local mean / variance
    logger.info("Compute local CC 1/2")
    H, K, L, result = compute_local_cchalf(reflections, kernel_size)

    # Get the coordinats of points to show
    logger.info("Generating coordinates")
    X, Y, Z = (flex.mat3_double(len(H), A) * flex.vec3_double(H, K, L)).parts()

    # Get the colours
    logger.info("Map CC to colours")
    logger.info("Max CC 1/2: %f" % max(result))
    logger.info("Min CC 1/2: %f" % min(result))
    logger.info("Min H: %d" % min(H))
    logger.info("Max H: %d" % max(H))
    logger.info("Min K: %d" % min(K))
    logger.info("Max K: %d" % max(K))
    logger.info("Min L: %d" % min(L))
    logger.info("Max L: %d" % max(L))
    colors = cm.plasma.colors
    min_r = 0
    max_r = 1.0
    B = 255.0 / (max_r - min_r)
    A = -B * min_r
    rgba = []
    for r in result:
        i = A + B * r
        if i < 0:
            i = 0
        if i > 255:
            i = 255
        alpha = 0.1 + 0.9 * (1.0 - i / 255.0)
        rgba.append(colors[int(i)] + [alpha])

    # Display the points
    logger.info("Displaying")
    xyz = list(zip(X, Y, Z))
    meanX = sum(X)/len(X)
    meanY = sum(Y)/len(Y)
    meanZ = sum(Z)/len(Z)
    rangeX = max(X)-min(X)
    rangeY = max(Y)-min(Y)
    rangeZ = max(Z)-min(Z)
    distance = 3*max([rangeX, rangeY, rangeZ])
    v = pptk.viewer(xyz, rgba)
    v.set(show_grid=False)
    v.set(bg_color=[0, 0, 0, 0])
    v.set(point_size=point_size)
    v.set(lookat=[meanX, meanY, meanZ])
    v.set(r=distance)
    v.set(theta=0)

    # If record
    if record:
        if not os.path.exists(directory):
            os.mkdir(directory)
        v.wait()
        lookatX, lookatY, lookatZ = v.get("lookat")
        distance = v.get("r")
        theta = v.get("theta")
        poses = []
        poses.append([lookatX, lookatY, lookatZ, 0 * np.pi/2, theta, distance])
        poses.append([lookatX, lookatY, lookatZ, 1 * np.pi/2, theta, distance])
        poses.append([lookatX, lookatY, lookatZ, 2 * np.pi/2, theta, distance])
        poses.append([lookatX, lookatY, lookatZ, 3 * np.pi/2, theta, distance])
        poses.append([lookatX, lookatY, lookatZ, 4 * np.pi/2, theta, distance])
        #v.play("images", poses, 2 * np.arange(5), repeat=True, interp='linear') 
        v.record(directory, poses, 2 * np.arange(5), interp='cubic_natural') 


def match_reflections(reflections, reference):
    '''
    Match two sets of reflections

    '''
    lookup1 = {}
    lookup2 = {}
    for i, h in enumerate(reflections['asu_miller_index']):
        lookup1[h] = i
    for i, h in enumerate(reference['asu_miller_index']):
        if h in lookup1:
            lookup2[h] = i
    H = flex.int()
    K = flex.int()
    L = flex.int()
    I1 = flex.double()
    I2 = flex.double()
    intensity1 = reflections['intensity']
    intensity2 = reference['intensity']
    for h in lookup2:
        H.append(h[0])
        K.append(h[1])
        L.append(h[2])
        I1.append(intensity1[lookup1[h]])
        I2.append(intensity2[lookup2[h]])
    return H, K, L, I1, I2

def compute_local_cc_vs_ref(reflections, reference, kernel_size):
    '''
    Compute the local cc vs a reference set of intrnsities

    '''

    H, K, L, I1, I2 = match_reflections(reflections, reference)

    # Get the range of miller indices
    min_H, max_H = min(H), max(H)
    min_K, max_K = min(K), max(K)
    min_L, max_L = min(L), max(L)

    n_H = max_H - min_H + 1
    n_K = max_K - min_K + 1
    n_L = max_L - min_L + 1

    
    mask = flex.bool(flex.grid(n_L, n_K, n_H))
    data1 = flex.double(mask.accessor())
    data2 = flex.double(mask.accessor())
    H_array = flex.double(mask.accessor())
    K_array = flex.double(mask.accessor())
    L_array = flex.double(mask.accessor())
    for h, k, l, intensity1, intensity2 in zip(H, K, L, I1, I2):
        x = h - min_H
        y = k - min_K
        z = l - min_L
        assert x >= 0 and y >= 0 and z >= 0
        mask[z, y, x] = True
        data1[z, y, x] = intensity1
        data2[z, y, x] = intensity2
        H_array[z, y, x] = h
        K_array[z, y, x] = k
        L_array[z, y, x] = l

    indices = []
    for k in range(mask.all()[0]):
        for j in range(mask.all()[1]):
            for i in range(mask.all()[2]):
                if mask[k, j, i] > 0:
                    indices.append((k, j, i))

    cc_array = flex.double(mask.accessor())
    mask_out = flex.bool(mask.accessor())
    for k, j, i in indices:
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

        X = flex.double()
        Y = flex.double()
        for kk in range(k0, k1):
            for jj in range(j0, j1):
                for ii in range(i0, i1):
                    if mask[kk, jj, ii]:
                        X.append(data1[kk,jj,ii])
                        Y.append(data2[kk,jj,ii])
        if len(X) > 0.25 * (2 * kernel_size + 1) ** 3:
            c = flex.linear_correlation(X, Y)
            cc_array[k,j,i] = c.coefficient()
            print c.coefficient()
            mask_out[k,j,i] = True

    selection = (mask_out == True).as_1d()
    H_sub = H_array.as_1d().select(selection)
    K_sub = K_array.as_1d().select(selection)
    L_sub = L_array.as_1d().select(selection)
    R_sub = cc_array.as_1d().select(selection)

    return H_sub, K_sub, L_sub, R_sub


def plot_local_cc_vs_ref(reflections, reference, A, kernel_size=5, record=False, directory="output", point_size=0.0025):
    '''
    Plot the local CC vs a reference
    '''

    # Compute the local CC
    logger.info("Compute local CC vs reference")
    H, K, L, result = compute_local_cc_vs_ref(reflections, reference, kernel_size)

    # Get the coordinats of points to show
    logger.info("Generating coordinates")
    X, Y, Z = (flex.mat3_double(len(H), A) * flex.vec3_double(H, K, L)).parts()

    # Get the colours
    logger.info("Map CC to colours")
    logger.info("Max CC 1/2: %f" % max(result))
    logger.info("Min CC 1/2: %f" % min(result))
    logger.info("Min H: %d" % min(H))
    logger.info("Max H: %d" % max(H))
    logger.info("Min K: %d" % min(K))
    logger.info("Max K: %d" % max(K))
    logger.info("Min L: %d" % min(L))
    logger.info("Max L: %d" % max(L))
    colors = cm.plasma.colors
    min_r = 0
    max_r = 1.0
    B = 255.0 / (max_r - min_r)
    A = -B * min_r
    rgba = []
    for r in result:
        i = A + B * r
        if i < 0:
            i = 0
        if i > 255:
            i = 255
        alpha = 0.1 + 0.9 * (1.0 - i / 255.0)
        rgba.append(colors[int(i)] + [alpha])

    # Display the points
    logger.info("Displaying")
    xyz = list(zip(X, Y, Z))
    meanX = sum(X)/len(X)
    meanY = sum(Y)/len(Y)
    meanZ = sum(Z)/len(Z)
    rangeX = max(X)-min(X)
    rangeY = max(Y)-min(Y)
    rangeZ = max(Z)-min(Z)
    distance = 3*max([rangeX, rangeY, rangeZ])
    v = pptk.viewer(xyz, rgba)
    v.set(show_grid=False)
    v.set(bg_color=[0, 0, 0, 0])
    v.set(point_size=point_size)
    v.set(lookat=[meanX, meanY, meanZ])
    v.set(r=distance)
    v.set(theta=0)

    # If record
    if record:
        if not os.path.exists(directory):
            os.mkdir(directory)
        v.wait()
        distance = v.get("r")
        theta = v.get("theta")
        lookatX, lookatY, lookatZ = v.get("lookat")
        poses = []
        poses.append([lookatX, lookatY, lookatZ, 0 * np.pi/2, theta, distance])
        poses.append([lookatX, lookatY, lookatZ, 1 * np.pi/2, theta, distance])
        poses.append([lookatX, lookatY, lookatZ, 2 * np.pi/2, theta, distance])
        poses.append([lookatX, lookatY, lookatZ, 3 * np.pi/2, theta, distance])
        poses.append([lookatX, lookatY, lookatZ, 4 * np.pi/2, theta, distance])
        #v.play("images", poses, 2 * np.arange(5), repeat=True, interp='linear') 
        v.record(directory, poses, 2 * np.arange(5), interp='cubic_natural') 
