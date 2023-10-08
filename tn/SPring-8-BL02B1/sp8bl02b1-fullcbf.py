import numpy as np
import sys

from dxtbx.format.FormatCBFMini import FormatCBFMini

# The CBF header does not contain Distance and Start_angle
# so we have to avoid reading them.


def dummy_func(self):
    return None


FormatCBFMini._detector = dummy_func
FormatCBFMini._scan = dummy_func

if len(sys.argv) != 3:
    sys.stderr.write("Usage: sp8bl02b1-fullcbf.py input.cbf output.cbf\n")
    exit(-1)

cbffile = sys.argv[1]
outfile = sys.argv[2]
assert cbffile.endswith(".cbf")
inffile = cbffile[:-3] + "inf"

header_items = {}
with open(inffile, "r") as inf:
    inf.readline()  # skip {
    for line in inf:
        line = line.rstrip()
        if line == "}":
            break

        assert line.endswith(";")
        idx = line.index("=")
        header_items[line[:idx]] = line[(idx + 1) : -1]

reader = FormatCBFMini(cbffile)
# print(reader._cif_header_dictionary)
# print(header_items)

fullcbf_out = """
# Takanori Nakane (Institute for Protein Research, Osaka University)
# constructed this full CBF header using 1191_00005.cbf in CBFLib
# as an example.

_diffrn.id SPRING8_BL02B1
_diffrn.crystal_id xtal001

loop_
_diffrn_source.diffrn_id
_diffrn_source.source
_diffrn_source.type
SPRING8_BL02B1 synchrotron 'SPring-8 Beamline BL02B1'

loop_
_diffrn_radiation.diffrn_id
_diffrn_radiation.wavelength_id
_diffrn_radiation.monochromator
SPRING8_BL02B1 WAVELENGTH1 'Si 111'

# FIXME: I don't know these values.
# diffrn_radiation.polarizn_source_ratio
# diffrn_radiation.polarizn_source_norm
# diffrn_radiation.div_x_source
# diffrn_radiation.div_y_source
# diffrn_radiation.div_x_y_source

loop_
_diffrn_detector.diffrn_id
_diffrn_detector.id
_diffrn_detector.type
_diffrn_detector.number_of_axes
SPRING8_BL02B1 None 'pilatus' 4

loop_
_diffrn_detector_axis.detector_id
_diffrn_detector_axis.axis_id
None DETECTOR_2THETA
None DETECTOR_X
None DETECTOR_Y
None DETECTOR_Z

loop_
_diffrn_detector_element.id
_diffrn_detector_element.detector_id
ELEMENT1 None

loop_
_diffrn_data_frame.id
_diffrn_data_frame.detector_element_id
_diffrn_data_frame.array_id
_diffrn_data_frame.binary_id
FRAME1 ELEMENT1 ARRAY1 1

loop_
_diffrn_scan.id
_diffrn_scan.frame_id_start
_diffrn_scan.frame_id_end
_diffrn_scan.frames
SCAN1 FRAME1 FRAME1 1

loop_
_diffrn_measurement.diffrn_id
_diffrn_measurement.id
_diffrn_measurement.number_of_axes
_diffrn_measurement.method
SPRING8_BL02B1 GONIOMETER 3 rotation

# Removed `diffrn_measurement.sample_detector_distance` because this is not present in the CIF dictionary.

loop_
_diffrn_measurement_axis.measurement_id
_diffrn_measurement_axis.axis_id
GONIOMETER GONIOMETER_OMEGA
GONIOMETER GONIOMETER_KAPPA
GONIOMETER GONIOMETER_PHI
"""

wavelength = float(header_items["SCAN_WAVELENGTH"])

fullcbf_out += (
    """
loop_
_diffrn_radiation_wavelength.id
_diffrn_radiation_wavelength.wavelength
_diffrn_radiation_wavelength.wt
WAVELENGTH1 %f 1
"""
    % wavelength
)

assert header_items["CRYSTAL_GONIO_NAMES"] == " Omega Chi Phi"
assert header_items["ROTATION_AXIS_NAME"] == "Omega"
omega, kappa, phi = [float(x) for x in header_items["CRYSTAL_GONIO_VALUES"].split()]
# I don't know why they didn't use the second item in CCD_GONIO_VALUES...
_, two_theta, distance = [float(x) for x in header_items["SCAN_DET_RELZERO"].split()]
# start, end, increment; for frame ROTATION, end == increment
# See https://www.rigaku.com/downloads/software/free/dTREK%20Image%20Format%20v1.1.pdf.
omega_incr = float(header_items["ROTATION"].split()[2])

fullcbf_out += """
loop_
_diffrn_scan_axis.scan_id
_diffrn_scan_axis.axis_id
_diffrn_scan_axis.angle_start
_diffrn_scan_axis.angle_range
_diffrn_scan_axis.angle_increment
_diffrn_scan_axis.displacement_start
_diffrn_scan_axis.displacement_range
_diffrn_scan_axis.displacement_increment
SCAN1 GONIOMETER_OMEGA %.4f %.4f %.4f 0.0 0.0 0.0
SCAN1 GONIOMETER_KAPPA %.4f 0.0000 0.0000 0.0 0.0 0.0
SCAN1 GONIOMETER_PHI %.4f 0.0000 0.0000 0.0 0.0 0.0
SCAN1 DETECTOR_2THETA %.4f 0.0 0.0 0.0 0.0 0.0
SCAN1 DETECTOR_Z 0.0 0.0 0.0 %.2f 0.0 0.0
SCAN1 DETECTOR_Y 0.0 0.0 0.0 0.0 0.0 0.0
SCAN1 DETECTOR_X 0.0 0.0 0.0 0.0 0.0 0.0

# Because this file contains only a single frame, angle_start == angle_range.
# FIXME: Am I correct?
""" % (
    omega,
    omega_incr,
    omega_incr,
    kappa,
    phi,
    two_theta,
    distance,
)

exp_time = float(reader._cif_header_dictionary["Exposure_time"].split()[0])
exp_period = float(reader._cif_header_dictionary["Exposure_period"].split()[0])
date = reader._cif_header_dictionary["timestamp"]

fullcbf_out += """
loop_
_diffrn_scan_frame.frame_id
_diffrn_scan_frame.frame_number
_diffrn_scan_frame.integration_time
_diffrn_scan_frame.exposure_time
_diffrn_scan_frame.scan_id
_diffrn_scan_frame.date
FRAME1 1 %f %f SCAN1 %s

# `diffrn_scan_frame.date` is not mandatory in imgCIF but DIALS requires it.
# `diffrn_scan_frame.integration_time` was taken from the miniCIF header `Exposure_time`.
# `diffrn_scan_frame.exposure_time` (non-standard in imgCIF) was from the miniCIF header `Exposure_period`.
# See also https://www.iucr.org/__data/iucr/lists/imgcif-archive/msg00349.html
""" % (
    exp_time,
    exp_period,
    date,
)

fullcbf_out += """
loop_
_diffrn_scan_frame_axis.frame_id
_diffrn_scan_frame_axis.axis_id
_diffrn_scan_frame_axis.angle
_diffrn_scan_frame_axis.displacement
FRAME1 GONIOMETER_OMEGA %.4f 0.0
FRAME1 GONIOMETER_KAPPA %.4f 0.0
FRAME1 GONIOMETER_PHI %.4f 0.0
FRAME1 DETECTOR_2THETA %.4f 0.0
FRAME1 DETECTOR_Z 0.0 %.2f
FRAME1 DETECTOR_Y 0.0 0.0
FRAME1 DETECTOR_X 0.0 0.0
""" % (
    omega,
    kappa,
    phi,
    two_theta,
    distance,
)

# Same information in "Pixel_size" (but in meter)
# below are in mm.
_, _, pixelsize_x, pixelsize_y = [
    float(x) for x in header_items["CCD_SPATIAL_DISTORTION_INFO"].split()
]
shift_x = float(header_items["CCD_SPATIAL_BEAM_POSITION"].split()[0]) * pixelsize_x
shift_y = -float(header_items["CCD_SPATIAL_BEAM_POSITION"].split()[1]) * pixelsize_y

fullcbf_out += """
loop_
_axis.id
_axis.type
_axis.equipment
_axis.depends_on
_axis.vector[1] _axis.vector[2] _axis.vector[3]
_axis.offset[1] _axis.offset[2] _axis.offset[3]
GONIOMETER_OMEGA rotation goniometer . 1 0 0 . . .
GONIOMETER_KAPPA rotation goniometer GONIOMETER_OMEGA 0 0 -1 . . .
GONIOMETER_PHI   rotation goniometer GONIOMETER_KAPPA 1 0 0 . . .
SOURCE           general source . 0 0 1 . . .
GRAVITY          general gravity . 0 1 0 . . .
DETECTOR_2THETA  rotation    detector . 1 0 0 . . . 
DETECTOR_Z       translation detector DETECTOR_2THETA 0 0 -1 0 0 0
DETECTOR_Y       translation detector DETECTOR_Z 0 1 0 0 0 0
DETECTOR_X       translation detector DETECTOR_Y 1 0 0 0 0 0
ELEMENT_X        translation detector DETECTOR_X -1 0 0 %.2f %.2f 0
ELEMENT_Y        translation detector ELEMENT_X 0 1 0 0 0 0

# The unit of axis_offset is mm.
""" % (
    shift_x,
    shift_y,
)

# Same information in CCD_DETECTOR_DIMENSIONS or SIZE1/SIZE2
fs = int(reader._cif_header_dictionary["X-Binary-Size-Fastest-Dimension"])
ss = int(reader._cif_header_dictionary["X-Binary-Size-Second-Dimension"])

fullcbf_out += """
# According to the imgCIF definition (https://www.iucr.org/__data/iucr/cifdic_html/2/cif_img.dic/Caxis.html),
# +Z is along the beam, from the sample to the source.
# +X is along the omega axis.
# Y completes the right handed system.
#
# ELEMENT_X/Y: how a movement from a pixel to the neighbour changes the coordinate.
#              So this is multiplication.
# DETECTOR_X/Y: this is addition.

loop_
_array_structure_list.array_id
_array_structure_list.index
_array_structure_list.dimension
_array_structure_list.precedence
_array_structure_list.direction
_array_structure_list.axis_set_id
ARRAY1 1 %d 1 increasing ELEMENT_X
ARRAY1 2 %d 2 increasing ELEMENT_Y

# CHECKME: Is `array_structure_list.direction` ignored in DIALS?
# Changing this does not flip images...
""" % (
    fs,
    ss,
)

fullcbf_out += """
loop_
_array_structure_list_axis.axis_set_id
_array_structure_list_axis.axis_id
_array_structure_list_axis.displacement
_array_structure_list_axis.displacement_increment
ELEMENT_X ELEMENT_X 0.0 %f
ELEMENT_Y ELEMENT_Y 0.0 %f

# We can change the direction by setting a negative `displacement_increment`
# but we should rather use ELEMENT_X/Y.

loop_
_array_element_size.array_id
_array_element_size.index
_array_element_size.size
ARRAY1 1 %f
ARRAY1 2 %f

# `array_element_size.size` is in m, while `array_structure_list_axis.displacement_increment` is in mm!
""" % (
    pixelsize_x,
    pixelsize_y,
    pixelsize_x / 1000.0,
    pixelsize_y / 1000.0,
)

cbf_cutoff = int(reader._cif_header_dictionary["Count_cutoff"].split()[0])
cutoff = int(header_items["SATURATED_VALUE"])

fullcbf_out += """
loop_
_array_intensities.array_id
_array_intensities.binary_id
_array_intensities.linearity
_array_intensities.gain
_array_intensities.gain_esd
_array_intensities.overload
_array_intensities.undefined_value
ARRAY1 1 linear 1.0 . %d -1

# `array_intensities.overload` was taken from `SATURATED_VALUE`;
# the miniCBF header `Count_cutoff` says %d instead.
# FIXME: I am not sure which is correct.
""" % (
    cutoff,
    cbf_cutoff,
)

fullcbf_out += """
loop_
_array_structure.id
_array_structure.encoding_type
_array_structure.compression_type
_array_structure.byte_order
ARRAY1 "signed 32-bit integer" CBF_BYTE_OFFSET little_endian

# The miniCBF header below is preserved from the original input.
# This does not mention goniometer angles.
"""

out = open(outfile, "wb")

# Copy up to the first block
with open(cbffile, "rb") as f:
    while True:
        line = f.readline()
        if line is None:
            sys.stderr.write("ERROR: did not find any datablock")
            exit(-1)
        out.write(line)
        if line.decode("ascii").startswith("data"):
            break

    # Write full CBF header
    out.write(fullcbf_out.encode("ascii"))

    # Copy the rest
    out.write(f.read())
    out.close()
