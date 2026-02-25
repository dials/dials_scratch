import argparse
from datetime import datetime, timedelta
import h5py
import math
import numpy as np
from pathlib import Path, PurePath
import sys

from dxtbx.format.FormatHDF5EigerNearlyNexus import EigerNXmxFixer
from dxtbx.util import ersatz_uuid4

parser = argparse.ArgumentParser(description="Convert EIGER master files from SPring-8 BL05XU into NXmx files.")

parser.add_argument("input_file", help="EIGER master file name")
parser.add_argument("output_file", help="Output NXmx file name")
parser.add_argument("--frame_time", help="Frame time (sec/frame)", type=float, default=0.1, required=True)
parser.add_argument("--omega_start", help="Starting OMEGA angle (deg)", type=float, default=0.0)
parser.add_argument("--omega_step", help="OMEGA angle per frame (deg/frame)", type=float, default=0.5, required=True)
parser.add_argument("--phi", help="PHI angle (deg)", type=float, default=0)
parser.add_argument("--kappa", help="KAPPA angle (deg). Usually 45.", type=float, default=45)
parser.add_argument("--two_theta", help="Detector TWO THETA angle (deg)", type=float, default=0)
parser.add_argument("--distance", help="Detector distance (m)", type=float, default=0.06,required=True)
parser.add_argument("--wavelength", help="Wavelenegth (angstrom", type=float, default=0.62, required=True)
parser.add_argument("--beam_center_x", help="Beam center (fast axis) at TWO THETA = 0 (px)", type=float, default=384.91)
parser.add_argument("--beam_center_y", help="Beam center (slow axis) at TWO THETA = 0 (px)", type=float, default=554.48)

args = parser.parse_args()

frame_time = args.frame_time
omega_start = args.omega_start
omega_step = args.omega_step
phi = args.phi
kappa = args.kappa
two_theta = args.two_theta
distance = args.distance
wavelength = args.wavelength
beam_center_x = args.beam_center_x
beam_center_y = args.beam_center_y

print("frame_time: %f sec" % frame_time)
print("omega: %f deg, step %f deg/frame" % (omega_start, omega_step))
print("phi: %f deg" % phi)
print("kappa: %f deg" % kappa)
print("two_theta: %f deg" % two_theta)
print("wavelength: %f A" % wavelength)
print("distance: %f m" % distance)
print("beam center at two_theta=0: (%f, %f) px (fast, slow)" % (beam_center_x, beam_center_y))

fn_in = args.input_file
fn_out = args.output_file
orig = h5py.File(fn_in, "r")

basedir = Path(fn_in).parent.absolute()
outputdir = Path(fn_out).parent.absolute()
if basedir != outputdir:
    sys.stderr.write(
        "The output file must be in the same directory as the input files.\n"
    )
    sys.stderr.write("If not, please make symbolic links.\n")
    exit(-1)

if False:
    print("Omega:", orig["entry"]["sample"]["goniometer"]["omega"][()])
    print("Chi:", orig["entry"]["sample"]["goniometer"]["chi"][()])
    print("Kappa:", orig["entry"]["sample"]["goniometer"]["kappa"][()])
    print("Phi:", orig["entry"]["sample"]["goniometer"]["phi"][()])
    print(
        "Two theta:",
        orig["entry"]["instrument"]["detector"]["goniometer"]["two_theta"][()],
    )
    print(
        "beam_center_x:", orig["entry"]["instrument"]["detector"]["beam_center_x"][()]
    )
    print(
        "beam_center_y:", orig["entry"]["instrument"]["detector"]["beam_center_y"][()]
    )
    print(
        "detector_distance:",
        orig["entry"]["instrument"]["detector"]["detector_distance"][()],
    )

temp_file = "tmp_master_%s.nxs" % ersatz_uuid4()
fixed = EigerNXmxFixer(fn_in, temp_file).handle

# BL05XU layout

# looking from the source to the crystal,
#  the OMEGA axis and the TWO_THETA axis are horizontal.
#     OMEGA is clockwise seen from the right.
#     TWO_THETA is the opposite.
#  the PHI axis (at OMEGA 0) is from bottom left to upper right (i.e. fixed at 45 degrees)
#  the detector origin is at the top left
#  the fast axis is horizontal and towards right
#  the slow axis is vertical and downwards
#
# In NeXus/McStats coordinate system,
#  +Z is from the source to the crystal
#  +X is horizontal, leftwards
#  +Y completes the right hand system, so vertical, upwards
#
# Thus,
# OMEGA is (1, 0, 0)
# KAPPA is (0, 0, -1), fixed at 45 deg
# PHI is (1, 0, 0) at KAPPA = 0
# TWO_THETA is (-1, 0, 0)
# fast and slow (at TWO_THETA=0) are (-1, 0, 0), (0, -1, 0)
#
# When importing, DIALS converts this into:
#  +Z is flipped (from the crystal to the source)
#  +X is flipped (horizontal, rightwards)
#  +Y as is (vertical, upwards)
#  Note that this preserves the hand.

number_of_frames = fixed["/entry/sample/goniometer/omega"].shape[0]

# Sample depends on phi, not omega
del fixed["/entry/sample/depends_on"]
fixed["/entry/sample/depends_on"] = np.string_("/entry/sample/transformations/phi")

# Set up two theta
fixed.copy(
    "/entry/instrument/detector/goniometer/two_theta",
    "/entry/instrument/detector/transformations/two_theta",
)
fixed["/entry/instrument/detector/transformations/two_theta"][()] = [two_theta] * number_of_frames
fixed["/entry/instrument/detector/transformations/two_theta"].attrs["vector"] = (
    -1.0,
    0.0,
    0.0,
)
fixed["/entry/instrument/detector/transformations/two_theta"].attrs[
    "units"
] = np.string_("degree")
fixed["/entry/instrument/detector/transformations/two_theta"].attrs[
    "transformation_type"
] = np.string_("rotation")
fixed["/entry/instrument/detector/transformations/two_theta"].attrs["offset"] = (
    0.0,
    0.0,
    0.0,
)
fixed["/entry/instrument/detector/transformations/two_theta"].attrs[
    "depends_on"
] = np.string_(".")
fixed["/entry/instrument/detector/transformations/translation"].attrs[
    "depends_on"
] = np.string_("/entry/instrument/detector/transformations/two_theta")

# Set up phi
fixed.copy("/entry/sample/goniometer/phi", "/entry/sample/transformations/phi")
fixed["/entry/sample/transformations/phi"][()] = [phi] * number_of_frames
fixed["/entry/sample/transformations/phi"].attrs["vector"] = (-1.0, 0.0, 0.0)
fixed["/entry/sample/transformations/phi"].attrs["units"] = np.string_("degree")
fixed["/entry/sample/transformations/phi"].attrs["transformation_type"] = np.string_(
    "rotation"
)
fixed["/entry/sample/transformations/phi"].attrs["offset"] = (0.0, 0.0, 0.0)
fixed["/entry/sample/transformations/phi"].attrs["depends_on"] = np.string_(
    "/entry/sample/transformations/kappa"
)

# Set up kappa
fixed.copy("/entry/sample/goniometer/kappa", "/entry/sample/transformations/kappa")
fixed["/entry/sample/transformations/kappa"][()] = [kappa] * number_of_frames
fixed["/entry/sample/transformations/kappa"].attrs["vector"] = (0.0, 0.0, -1.0)
fixed["/entry/sample/transformations/kappa"].attrs["units"] = np.string_("degree")
fixed["/entry/sample/transformations/kappa"].attrs["transformation_type"] = np.string_(
    "rotation"
)
fixed["/entry/sample/transformations/kappa"].attrs["offset"] = (0.0, 0.0, 0.0)
fixed["/entry/sample/transformations/kappa"].attrs["depends_on"] = np.string_(
    "/entry/sample/transformations/omega"
)

# Update omega
# I don't know why but EigerNXmxFixer recalculates omega from /entry/sample/goniometer/omega_range_average,
# not using per-frame values
fixed["/entry/sample/transformations/omega"][()] = [omega_start + i * omega_step for i in range(number_of_frames)]
fixed["/entry/sample/transformations/omega"].attrs["vector"] = (1.0, 0.0, 0.0)
fixed["/entry/sample/transformations/omega"].attrs["offset"] = (0.0, 0.0, 0.0)

# Set the wavelength
fixed["/entry/instrument/beam/incident_wavelength"][()] = wavelength

# Delete redundant information

# /entry/sample/beam is copied from /entry/instrument/beam by EigerNXmxFixer but actually the old place
# is correct according to NXmx.
del fixed["/entry/sample/beam"]
del fixed["/entry/sample/goniometer"]  # re-written in /entry/sample/transformations
del fixed[
    "/entry/instrument/detector/geometry"
]  # re-written in /entry/instrument/detector/transformations
del fixed[
    "/entry/instrument/detector/goniometer"
]  # re-written in /entry/instrument/detector/transformations

# Make the links in /entry/data relative (again!)
update = {}
for key in fixed["/entry/data"]:
    item = fixed["/entry/data"].get(key, getlink=True)
    if item.__class__ == h5py.ExternalLink:
        relpath = PurePath(item.filename).relative_to(basedir)
        update[key] = h5py.ExternalLink(relpath, item.path)

for key, newitem in update.items():
    del fixed["entry/data"][key]
    fixed["/entry/data"][key] = newitem

# Set the detector distance
print(fixed["/entry/instrument/detector/depends_on"][()])
fixed["/entry/instrument/detector/transformations/translation"][()] = distance
fixed["/entry/instrument/detector/transformations/translation"].attrs["vector"] = (
    0.0,
    0.0,
    1.0, # along the beam
)
fixed["/entry/instrument/detector/distance"] = distance # recommended by NXmx, not used actually
del fixed["/entry/instrument/detector/detector_distance"]

# Set the beam center
fixed["/entry/instrument/detector/transformations/translation"].attrs["offset"] = (
    beam_center_x * fixed["/entry/instrument/detector/module/fast_pixel_direction"][()],
    beam_center_y * fixed["/entry/instrument/detector/module/slow_pixel_direction"][()],
    0.0,
)
# these two are not used actually
fixed["/entry/instrument/detector/beam_center_x"][()] = beam_center_x
fixed["/entry/instrument/detector/beam_center_y"][()] = beam_center_y

# suppress expected "NX_BOOLEAN, got H5T_STD_I32LE"
for item in [
    "/entry/instrument/detector/countrate_correction_applied",
    "/entry/instrument/detector/pixel_mask_applied",
]:
    val = np.int8(fixed[item][()])
    del fixed[item]
    fixed[item] = val
    fixed[item].attrs["NX_CLASS"] = np.string_("NX_BOOLEAN")

# Add mandatory NXmx entries
fixed["/entry/instrument/name"] = np.string_("BL05XU")

fixed["/entry/start_time"] = fixed[
    "/entry/instrument/detector/detectorSpecific/data_collection_date"
][()]

start_time = datetime.fromisoformat(fixed["/entry/start_time"][()].decode("ascii"))

fixed["/entry/instrument/detector/frame_time"][()] = frame_time
collection_time = timedelta(
    seconds=fixed["/entry/instrument/detector/frame_time"][()] * number_of_frames
)
end_time_estimated = start_time + collection_time
fixed["/entry/end_time_estimated"] = np.string_(end_time_estimated.isoformat())

fixed["/entry/sample/name"] = np.string_("Unknown sample")
fixed["/entry/source/name"] = np.string_("SPring-8")
fixed["/entry/source/"].attrs["NX_class"] = np.string_("NXsource")

# /entry/detector/countrate_correction_lookup_table
#  FIXME: Is this really required!? Some detectors don't apply this correction.

out_file = h5py.File(fn_out, "w")
fixed.copy("entry", out_file)
out_file.close()
