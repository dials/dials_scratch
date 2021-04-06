import pprofile
import sys
import time

import h5py
import hdf5plugin  # noqa: F401
import numba
import numpy as np
from scipy import ndimage

import matplotlib.pyplot as plt
from matplotlib import colors

def scipy_dispersion(image, kernel_size):
    kernel = np.ones((kernel_size, kernel_size))

    max_valid = 65534
    mask = (image <= max_valid).astype(image.dtype)
    masked_image = image * mask
    masked_image2 = masked_image**2

    sum_image = ndimage.convolve(masked_image, kernel, mode="constant", cval=0)
    sum_sq = ndimage.convolve(masked_image2, kernel, mode="constant", cval=0)
    n = ndimage.convolve(mask, kernel, mode="constant", cval=0)

    mean_image = np.zeros_like(image, dtype=np.float)
    np.divide(sum_image, n, where=(n > 0), out=mean_image)

    inv_count = np.zeros_like(image, dtype=np.float)
    np.divide(1, n, where=(n > 0), out=inv_count)

    variance_image = (sum_sq - inv_count * sum_image ** 2) * inv_count

    dispersion_index = np.zeros_like(image, dtype=np.float)
    np.divide(variance_image, mean_image, where=(mean_image > 0), out=dispersion_index)
    return dispersion_index


#------------------ SAT algorithms --------------
def sat_python(image):
    '''
    SAT with explicit loops
    '''
    sat = np.zeros_like(image)
    for i in range(image.shape[0]):
        for j in range(image.shape[1]):
            sat[i, j] = image[i, j] + sat[i, j - 1] + sat[i - 1, j] - sat[i - 1, j - 1]
    return sat

@numba.njit(numba.int32[:, ::1](numba.int32[:, ::1]))
def sat_python_numba(image):
    '''
    SAT with explicit loops helped by numba
    '''
    sat = np.zeros_like(image)
    for i in range(image.shape[0]):
        for j in range(image.shape[1]):
            sat[i, j] = image[i, j] + sat[i, j - 1] + sat[i - 1, j] - sat[i - 1, j - 1]
    return sat


def sat_numpy(image):
    '''
    SAT is a double cumsum
    '''
    sat = np.zeros_like(image)
    return image.cumsum(axis=0, out=sat).cumsum(axis=1, out=sat)


#-------------- end of SAT algorithms ---------------

def kernel_sum(image, kernel_size, sat_func):
    '''
    sum = satI(x1, y1) + sat(x0, y0) - I(x1, y0) - I(x0, y1)
    '''
    pad = (kernel_size - 1) // 2
    image = np.pad(image, (pad + 1, pad))
    sat = sat_func(image)
    return (
          sat[:-kernel_size, :-kernel_size ]  # top left
        + sat[kernel_size:,   kernel_size: ]  # bottom right
        - sat[kernel_size:,   :-kernel_size]  # top right
        - sat[:-kernel_size,  kernel_size: ]  # bottom left
    )


@numba.njit
def compute_variance(sum_image, sum_image_sq, n):
    return (sum_image_sq - np.square(sum_image) / n) / n


def calculate_dispersion_index(image, kernel_size, sat_func):
    max_valid = 65534
    mask = (image <= max_valid).astype(image.dtype)
    masked_image = image * mask
    masked_image2 = masked_image**2

    sum_image = kernel_sum(masked_image, kernel_size, sat_func)
    sum_sq = kernel_sum(masked_image2, kernel_size, sat_func)
    n = kernel_sum(mask, kernel_size, sat_func)

    mean_image = np.divide(sum_image, n, where=(n > 0))

    variance_image = compute_variance(sum_image, sum_sq, n)

    return np.divide(variance_image, mean_image, where=(mean_image > 0))


def profile(sat_function, image):
    prof = pprofile.Profile()
    #with prof():
    t0 = time.perf_counter()

    calculate_dispersion_index(image, kernel_size=7, sat_func=sat_function)
    #scipy_dispersion(image, kernel_size=7)

    t1 = time.perf_counter()
    print(f"thresholding with {sat_function.__name__} took {t1 -t0:.4f} seconds")
    #prof.dump_stats("profile.txt")


if __name__ == "__main__":
    assert len(sys.argv[1:]) == 1
    with h5py.File(sys.argv[1], mode="r") as f:
        image = f["/entry/data/data"][0]

    #plt.figure()
    #plt.imshow(image, norm=colors.SymLogNorm(1))

    #profile(sat_python, image)
    profile(sat_python_numba, image)
    #profile(sat_numpy, image)




