#!/usr/bin/env dials.python
# coding: utf-8

"""
Plot the magnitude of the refined goniometer rotation from a refined.expt file.

A simple utility for plotting the magnitude of the refined goniometer rotation
versus frame number, as measured from the orientation in the first frame.
"""

import argparse

import numpy as np
from scipy.spatial.transform import Rotation

from dxtbx.model import Crystal, ExperimentList

import matplotlib

matplotlib.use("Agg")

from matplotlib import pyplot as plt
from matplotlib import ticker

try:
    from typing import Iterable
except ImportError:
    pass


def get_angles(crystal):  # type: (Crystal) -> (np.array, np.array)
    """
    Get the magnitude of the refined goniometer rotation from a crystal model.

    Args:
        crystal: A dxtbx.model.Crystal object representing the refined crystal model.

    Returns:
        A tuple of arrays of equal length containing frame numbers and the magnitudes
        of the corresponding rotations.
    """
    # Create an array of frame numbers.
    images = np.arange(crystal.num_scan_points)

    # Extract the refined crystal orientation matrices.
    # We must force int(image) here because DIALS uses old Boost that doesn't play
    # nice with NumPy on Python 3.
    rotations = np.array([crystal.get_U_at_scan_point(int(image)) for image in images])
    # Recast these as rotation operators.
    rotations = Rotation.from_dcm(rotations.reshape((images.size, 3, 3)))
    # Reframe the rotations as being from the initial orientation.
    rotations_from_start = rotations[0].inv() * rotations

    # Extract the magnitude of the sample rotation for each frame.
    angles = np.linalg.norm(rotations_from_start.as_rotvec(), axis=1)

    return images, angles


def plot_angles(images, angles):  # type: (Iterable[int], Iterable[float]) -> None
    """
    Plot tumbling angle versus frame number.

    Args:
        images:  An array of frame numbers (zero-based).
        angles:  An array of corresponding tumbling angles.
    """
    fig, ax = plt.subplots()

    # Double-check that images is a NumPy array so that we can do broadcast arithmetic.
    images = np.array(images)

    # Plot the rotation magnitude (in degrees) versus image number (counted from 1).
    ax.plot(images + 1, np.degrees(angles), ".")
    ax.set_xlabel("Frame number")
    ax.set_ylabel("Sample tumbling angle")

    # Get the frame number axis tick locations.
    locations = ax.xaxis.get_major_locator().tick_values(images.min(), images.max())
    # Replace the tick at the 0th frame (which is meaningless) with one at the 1st.
    locations = np.where(locations == 0, 1, locations)
    ax.xaxis.set_major_locator(ticker.FixedLocator(locations))
    # Format the rotation magnitude (tumbling angle) axis tick labels in degrees.
    ax.yaxis.set_major_formatter(ticker.PercentFormatter(symbol=u"°", decimals=1))

    fig.savefig("tumbling_angle")


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "experiment_list",
    metavar="refined.expt",
    type=str,
    help="A refined experiment list from DIALS, such as 'refined.expt'.",
)


if __name__ == "__main__":
    experiment_list = ExperimentList.from_file(parser.parse_args().experiment_list)

    if len(experiment_list) > 1:
        print(
            "Warning: You have supplied an input file containing multiple "
            "experiment models.  Only the first will be used."
        )

    print("Saving a plot of tumbling angle versus image number as tumbling_angle.png.")
    plot_angles(*get_angles(experiment_list[0].crystal))
