import pprofile
import sys
import time

import h5py
import hdf5plugin  # noqa: F401
import numba
import numpy as np

#------------------ SAT algorithms --------------
def sat_python(image):
    '''
    SAT with explicit loops
    '''
    sat = np.zeros(image.shape, dtype=image.dtype)
    for i in range(image.shape[0]):
        for j in range(image.shape[1]):
            sat[i, j] = image[i, j] + sat[i, j - 1] + sat[i - 1, j] - sat[i - 1, j - 1]
    return sat

@numba.njit(numba.int32[:, ::1](numba.int32[:, ::1]))
def sat_python_numba(image):
    '''
    SAT with explicit loops helped by numba
    '''
    sat = np.zeros(image.shape, dtype=image.dtype)
    for i in range(image.shape[0]):
        for j in range(image.shape[1]):
            sat[i, j] = image[i, j] + sat[i, j - 1] + sat[i - 1, j] - sat[i - 1, j - 1]
    return sat


def sat_numpy(image):
    '''
    SAT is a double cumsum
    '''
    sat = np.empty_like(image, dtype=image.dtype)
    image.cumsum(axis=0, out=sat)
    return sat.cumsum(axis=1, out=sat)


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
    im = masked_image
    im2 = im ** 2

    sum_image = kernel_sum(masked_image, kernel_size, sat_func)
    sum_sq = kernel_sum(im2, kernel_size, sat_func)
    n = kernel_sum(mask, kernel_size, sat_func)
    mean_image = np.zeros(im.shape)
    np.divide(sum_image, n, where=(n > 0), out=mean_image)
    variance_image = compute_variance(sum_image, sum_sq, n)
    dispersion_index = np.ones(mean_image.shape)
    np.divide(variance_image, mean_image, where=(mean_image > 0), out=dispersion_index)


def profile(sat_function, image):
    prof = pprofile.Profile()
    with prof():
        t0 = time.perf_counter()
        for _ in range(100):
            calculate_dispersion_index(image, kernel_size=7, sat_func=sat_function)
        t1 = time.perf_counter()
        print(f"thresholding 100 times with {sat_function.__name__} took {t1 -t0:.4f} seconds")
    prof.dump_stats("profile.txt")


if __name__ == "__main__":
    assert len(sys.argv[1:]) == 1
    with h5py.File(sys.argv[1], mode="r") as f:
        image = f["/entry/data/data"][0]

    #profile(sat_python, image)
    #profile(sat_python_numba, image)
    profile(sat_numpy, image)




