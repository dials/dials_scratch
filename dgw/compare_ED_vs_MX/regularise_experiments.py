#!/usr/bin/env dials.python

# load original electron diffraction experiments_ED.json
from dxtbx.model.experiment_list import ExperimentListFactory

el = ExperimentListFactory.from_json_file("experiments_ED.json", check_format=False)
exp = el[0]
beam = exp.beam
detector = exp.detector
crystal = exp.crystal
goniometer = exp.goniometer
scan = exp.scan

# First fix the beam direction. Sample to source along Z.
beam.set_direction((0, 0, 1))

# Get information from the current detector
pixel_size = detector[0].get_pixel_size()
from scitbx import matrix

diagonal = matrix.col(detector[0].get_pixel_lab_coord((0, 0))) - matrix.col(
    detector[3].get_pixel_lab_coord(detector[3].get_image_size())
)
from math import sqrt

side_length = sqrt(diagonal.length() ** 2 / 2)
npx = int(round(side_length / pixel_size[0]))
beam_centre = [npx / 2 * pixel_size[0]] * 2
distance = detector[0].get_distance()

# Create a new single panel detector orthogonal to the beam, which intersects
# at the centre
from dxtbx.model.detector import DetectorFactory

detector2 = DetectorFactory.simple(
    sensor="PAD",
    distance=distance,
    beam_centre=beam_centre,
    fast_direction="+x",
    slow_direction="-y",
    pixel_size=pixel_size,
    image_size=(npx, npx),
)
exp.detector = detector2

# Make the scan a full turn to ensure no reflections get thrown out during
# refinement for being outside the scan range
image_width_deg = scan.get_oscillation(deg=True)[1]
nimages = int(round(360.0 / image_width_deg))
image_width_deg = 360.0 / nimages
image_range = 1, nimages
epochs = [0] * nimages
exposure_times = 0.0
oscillation = (0, image_width_deg)
from dxtbx.model import ScanFactory

exp.scan = ScanFactory.make_scan(
    image_range, exposure_times, oscillation, epochs, deg=True
)

el.as_file("experiments_ED_regularised.json")

# Now regularize to standard MX geometry
# Set beam energy to 12 keV
# lambda = h*c / E
# lambda = 1.98645e-25 / 1.92261e-15
beam.set_wavelength(1.0332)

# Pilatus 6M-like
pixel_size = (0.172, 0.172)
image_size = (2463, 2527)
beam_centre = (image_size[0] / 2 * pixel_size[0], image_size[1] / 2 * pixel_size[1])
distance = 200

detector3 = DetectorFactory.simple(
    sensor="PAD",
    distance=distance,
    beam_centre=beam_centre,
    fast_direction="+x",
    slow_direction="-y",
    pixel_size=pixel_size,
    image_size=image_size,
)
exp.detector = detector3

el.as_file("experiments_MX_regularised.json")
