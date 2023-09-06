from datetime import datetime, timedelta
import h5py
import math
import numpy as np
from pathlib import Path, PurePath
import sys

from dxtbx.format.FormatHDF5EigerNearlyNexus import EigerNXmxFixer
from dxtbx.util import ersatz_uuid4

if len(sys.argv) != 3:
    sys.stderr.write("Usage: dials.python sp8bl40xu.py input_master.h5 output.nxs\n")
    exit(-1)

fn_in = sys.argv[1]
fn_out = sys.argv[2]
orig = h5py.File(fn_in, "r")

basedir = Path(fn_in).parent.absolute()
outputdir = Path(fn_out).parent.absolute()
if basedir != outputdir:
    sys.stderr.write("The output file must be in the same directory as the input files.\n")
    sys.stderr.write("If not, please make symbolic links.\n")
    exit(-1)

if False:
    print("Omega:", orig["entry"]["sample"]["goniometer"]["omega"][()])
    print("Chi:", orig["entry"]["sample"]["goniometer"]["chi"][()])
    print("Kappa:", orig["entry"]["sample"]["goniometer"]["kappa"][()])
    print("Phi:", orig["entry"]["sample"]["goniometer"]["phi"][()])
    print("Two theta:", orig["entry"]["instrument"]["detector"]["goniometer"]["two_theta"][()])
    print("beam_center_x:", orig["entry"]["instrument"]["detector"]["beam_center_x"][()])
    print("beam_center_y:", orig["entry"]["instrument"]["detector"]["beam_center_y"][()])
    print("detector_distance:", orig["entry"]["instrument"]["detector"]["detector_distance"][()])

temp_file = "tmp_master_%s.nxs" % ersatz_uuid4()
fixed = EigerNXmxFixer(fn_in, temp_file).handle

# BL40XU layout

# looking from the source to the crystal,
#  the OMEGA axis and the TWO_THETA axis are vertical, downwards
#   (i.e. clockwise when looked down from the ceiling)
#  the PHI axis (at OMEGA 0) is towards source, lower (i.e. fixed at 45 degrees)
#   (i.e. clockwise when looked down from the ceiling)
#  the detector origin is at the top left
#  the fast axis is horizontal and towards right
#  the slow axis is vertical and downwards
#
# In NeXus/McStats coordinate system,
#  +Z is from the source to the crystal
#  +X is horizontal, leftwards
#  +Y completes the right hand system, so vertical, upwards
#
# Thus, OMEGA and TWO_THETA is (0, -1, 0)
# PHI is (0, -cos45, -sin45)
# fast and slow (at TWO_THETA=0) are (-1, 0, 0), (0, -1, 0)
#
# When importing, DIALS converts this into:
#  +Z is flipped (from the crystal to the source)
#  +X is flipped (horizontal, rightwards)
#  +Y as is (vertical, upwards)
#  Note that this preserves the hand.

# Sample depends on phi, not omega
del fixed["/entry/sample/depends_on"]
fixed["/entry/sample/depends_on"] = np.string_("/entry/sample/transformations/phi")

# Set up two theta
fixed.copy("/entry/instrument/detector/goniometer/two_theta", "/entry/instrument/detector/transformations/two_theta")
fixed["/entry/instrument/detector/transformations/two_theta"].attrs["vector"] = (0, -1, 0)
fixed["/entry/instrument/detector/transformations/two_theta"].attrs["units"] = np.string_("degree")
fixed["/entry/instrument/detector/transformations/two_theta"].attrs["transformation_type"] = np.string_("rotation")
fixed["/entry/instrument/detector/transformations/two_theta"].attrs["offset"] = 0
fixed["/entry/instrument/detector/transformations/two_theta"].attrs["depends_on"] = np.string_(".")
fixed["/entry/instrument/detector/transformations/translation"].attrs["depends_on"] = np.string_("/entry/instrument/detector/transformations/two_theta")

# Set up phi
fixed.copy("/entry/sample/goniometer/phi", "/entry/sample/transformations/phi")
#sin45 = math.sin(math.pi * 45 / 180.0) # = cos45
#fixed["/entry/sample/transformations/phi"].attrs["vector"] = (0, -sin45, -sin45)
# Refined values based on XRDa-159
fixed["/entry/sample/transformations/phi"].attrs["vector"] = (-0.0106, -0.7094, -0.7047)
fixed["/entry/sample/transformations/phi"].attrs["units"] = np.string_("degree")
fixed["/entry/sample/transformations/phi"].attrs["transformation_type"] = np.string_("rotation")
fixed["/entry/sample/transformations/phi"].attrs["offset"] = 0
fixed["/entry/sample/transformations/phi"].attrs["depends_on"] = np.string_("/entry/sample/transformations/omega")

# Update omega
# I don't know why but EigerNXmxFixer recalculates omega from /entry/sample/goniometer/omega_range_average,
# not using per-frame values
fixed["/entry/sample/transformations/omega"][()] = fixed["/entry/sample/goniometer/omega"][()]
fixed["/entry/sample/transformations/omega"].attrs["vector"] = (0, -1, 0)

# Delete redundant information

# /entry/sample/beam is copied from /entry/instrument/beam by EigerNXmxFixer but actually the old place
# is correct according to NXmx.
del fixed["/entry/sample/beam"]
del fixed["/entry/sample/goniometer"] # re-written in /entry/sample/transformations
del fixed["/entry/instrument/detector/geometry"] # re-written in /entry/instrument/detector/transformations
del fixed["/entry/instrument/detector/goniometer"] # re-written in /entry/instrument/detector/transformations

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

# Mandatory NXmx entries
fixed["/entry/start_time"] = fixed["/entry/instrument/detector/detectorSpecific/data_collection_date"][()]

start_time = datetime.fromisoformat(fixed["/entry/start_time"][()].decode("ascii"))
number_of_frames = fixed["/entry/sample/transformations/omega"].shape[0]
collection_time = timedelta(seconds=fixed["/entry/instrument/detector/frame_time"][()] * number_of_frames)
end_time_estimated = start_time + collection_time
fixed["/entry/end_time_estimated"] = np.string_(end_time_estimated.isoformat())

fixed["/entry/sample/name"] = np.string_("Unknown sample")
fixed["/entry/source/name"] = np.string_("SPring-8")
fixed["/entry/source/"].attrs["NX_class"] = np.string_("NXsource")

# /entry/detector/countrate_correction_lookup_tablea
#  FIXME: Is this really required!? Some detectors don't apply this correction.
#  Is this a bug in the NXmx definition?
# https://manual.nexusformat.org/classes/applications/NXmx.html says required but
# https://github.com/nexusformat/definitions/blob/main/applications/NXmx.nxdl.xml says
# minOccurs="0", i.e. optional.

out_file = h5py.File(fn_out, "w")
fixed.copy("entry", out_file)
out_file.close()
