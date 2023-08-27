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
    sys.stderr.write("Usage: cbf2smv.py input.cbf output.img\n")
    exit(-1)

cbffile = sys.argv[1]
reader = FormatCBFMini(cbffile)
image_data = np.array(reader.get_raw_data())
# print(image_data.shape, image_data.dtype)

assert cbffile.endswith(".cbf")
inffile = cbffile[:-3] + "inf"

header_items = []
two_theta = 0
with open(inffile, "r") as inf:
    inf.readline()  # skip {
    for line in inf:
        line = line.rstrip()
        if line == "}":
            break

        assert line.endswith(";")
        if line.startswith("DETECTOR_NAMES="):
            # This is used as the prefix for items like DETECTOR_DIMENSIONS.
            line = "DETECTOR_NAMES=CCD_;"
        elif line.startswith("CCD_DETECTOR_DESCRIPTION=") or line.startswith(
            "CCD_DETECTOR_IDENTIFICATION"
        ):
            # to satisfy FormatSMVRigakuSaturn.understand()
            line = line[:-1] + " fake saturn;"
        elif line.startswith("CCD_DETECTOR_VECTORS="):
            # the detector origin is at the top left, not lower left corner
            # the fast axis is antiparallel to the OMEGA axis
            # I *guess* this information is encoded in CCD_SPATIAL_DISTORTION_VECTORS
            # but somehow the below is necessary...
            line = "CCD_DETECTOR_VECTORS=0 1 0 1 0 0;"
        elif line.startswith("SCAN_DET_RELZERO="):
            # I don't know why they didn't use the second item in CCD_GONIO_VALUES...
            items = line[17:-1].split()
            two_theta = float(items[1])
        elif line.startswith("SATURATED_VALUE="):
            # The trusted range is [1,65533] (65535 == -1 == panel gaps, 65534 == -2 == bad pixels)
            line = "SATURATED_VALUE=65533;"
        elif line.startswith("CRYSTAL_GONIO_VECTORS="):
            # The CHI axis is along (0, 0, -1), i.e., parallel to the beam towards the detector, not (0, -1, 0)!!
            line = "CRYSTAL_GONIO_VECTORS=1.0000 0.0000 0.0000 0.0000 0.0000 -1.0000 1.0000 0.0000 0.0000;"

        header_items.append(line)

# fudge DTREK_DATE_TIME
header_items.append("DTREK_DATE_TIME=01-Jan-2000 00:00:00")

# Set two_theta
if two_theta != 0:
    for idx, item in enumerate(header_items):
        if item.startswith("CCD_GONIO_VALUES="):
            items = [float(x) for x in item[17:-1].split()]
            items[1] = two_theta
            header_items[idx] = (
                "CCD_GONIO_VALUES=" + " ".join([str(x) for x in items]) + ";"
            )
            break

smvfile = sys.argv[2]
with open(smvfile, "wb") as out:
    out.write(bytes("{\n", "ascii"))
    for item in header_items:
        out.write(bytes(item + "\n", "ascii"))
    out.write(bytes("}\n", "ascii"))

    out.seek(4096, 0)  # hard coded
    # FIXME: will this overflow?
    out.write(image_data.astype(np.uint16).tobytes())
