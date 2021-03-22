import pprofile
import sys
import time

import h5py
import hdf5plugin  # noqa: F401
import numba
import numpy as np


def summed_area_table(image):
    sat = np.empty(image.shape, dtype=image.dtype)
    image.cumsum(axis=0, out=sat)
    sat.cumsum(axis=1, out=sat)
    return sat


@numba.njit(numba.int32[:, ::1](numba.int32[:, ::1]))
def summed_area_table(image):
    sat = np.zeros(image.shape, dtype=image.dtype)
    for i in range(image.shape[0]):
        for j in range(image.shape[1]):
            sat[i, j] = image[i, j] + sat[i, j - 1] + sat[i - 1, j] - sat[i - 1, j - 1]
    return sat


def kernel_sum(image, kernel_size):
    pad = (kernel_size - 1) // 2
    image = np.pad(image, (pad + 1, pad))
    sat = summed_area_table(image)
    return (
        sat[:-kernel_size, :-kernel_size]
        + sat[kernel_size:, kernel_size:]  # top left
        - sat[kernel_size:, :-kernel_size]  # bottom right
        - sat[:-kernel_size, kernel_size:]  # top right  # bottom left
    )


@numba.njit
def compute_variance(sum_image, sum_image_sq, n):
    return (sum_image_sq - np.square(sum_image) / n) / n


def calculate_dispersion_index(image, kernel_size):
    max_valid = 65534
    mask = (image <= max_valid).astype(image.dtype)
    masked_image = image * mask
    im = masked_image
    im2 = im ** 2

    sum_image = kernel_sum(masked_image, kernel_size)
    sum_sq = kernel_sum(im2, kernel_size)
    n = kernel_sum(mask, kernel_size)
    mean_image = np.zeros(im.shape)
    np.divide(sum_image, n, where=(n > 0), out=mean_image)
    variance_image = compute_variance(sum_image, sum_sq, n)
    dispersion_index = np.ones(mean_image.shape)
    np.divide(variance_image, mean_image, where=(mean_image > 0), out=dispersion_index)


if __name__ == "__main__":
    assert len(sys.argv[1:]) == 1
    with h5py.File(sys.argv[1], mode="r") as f:
        image = f["/entry/data/data"][0]

    calculate_dispersion_index(image, kernel_size=7)
    prof = pprofile.Profile()
    with prof():
        t0 = time.perf_counter()
        calculate_dispersion_index(image, kernel_size=7)
        t1 = time.perf_counter()
        print(f"thresholding took {t1 -t0:.4f} seconds")
    prof.dump_stats("profile.txt")
