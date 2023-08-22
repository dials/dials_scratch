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
#print(image_data.shape, image_data.dtype)

assert cbffile.endswith(".cbf")
inffile = cbffile[:-3] + "inf"

header_items = []
with open(inffile, "r") as inf:
    inf.readline() # skip {
    for line in inf:
        line = line.rstrip()
        if line == "}": break

        assert line.endswith(";")
        if line.startswith("DETECTOR_NAMES="):
            # This is used as the prefix for items like DETECTOR_DIMENSIONS.
            line = "DETECTOR_NAMES=CCD_;"
        elif line.startswith("CCD_DETECTOR_DESCRIPTION=") or line.startswith("CCD_DETECTOR_IDENTIFICATION"):
            # to satisfy FormatSMVRigakuSaturn.understand()
            line = line[:-1] + " fake saturn;"
        elif line.startswith("CCD_GONIO_VALUES="):
            items = [float(x) for x in line[17:-1].split()]
            items[0] += 90 # the detector origin is at the top left, not lower left corner
            line = "CCD_GONIO_VALUES=" + " ".join([str(x) for x in items]) + ";"
        elif line.startswith("SATURATED_VALUE="):
            # The trusted range is [1,65533] (65535 == -1 == panel gaps, 65534 == -2 == bad pixels)
            line = "SATURATED_VALUE=65533;"

        header_items.append(line)

# fudge DTREK_DATE_TIME
header_items.append("DTREK_DATE_TIME=01-Jan-2000 00:00:00")

smvfile = sys.argv[2]
with open(smvfile, "wb") as out:
    out.write(bytes("{\n", "ascii"))
    for item in header_items:
        out.write(bytes(item + "\n", "ascii"))
    out.write(bytes("}\n", "ascii"))

    out.seek(4096, 0) # hard coded
    # FIXME: will this overflow?
    out.write(image_data.astype(np.uint16).tobytes())
