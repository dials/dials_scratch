#!/usr/bin/env cctbx.python

"""Load an XDS GAIN.cbf, flatten the values, then write out again"""
from __future__ import print_function

import binascii
import copy
import sys
from cbflib_adaptbx import compress, uncompress
from scitbx.array_family import flex


def squishGain(cbf_file, out_name, force_gain=None):

    start_tag = binascii.unhexlify("0c1a04d5")

    data = cbf_file.read()
    data_offset = data.find(start_tag) + 4
    cbf_header = data[: data_offset - 4]

    fast = 0
    slow = 0
    length = 0

    for record in cbf_header.split("\n"):
        if "X-Binary-Size-Fastest-Dimension" in record:
            fast = int(record.split()[-1])
        elif "X-Binary-Size-Second-Dimension" in record:
            slow = int(record.split()[-1])
        elif "X-Binary-Size:" in record:
            xbsize_record = record
            length = int(record.split()[-1])

    values = uncompress(
        packed=data[data_offset : data_offset + length], fast=fast, slow=slow
    )

    # remainder of the file, contains another CIF-BINARY-FORMAT-SECTION that looks
    # like just zero padding.
    tail = data[data_offset + length :]

    modified = copy.deepcopy(values)

    if force_gain:
        new_val = int(1000 * force_gain)
    else:
        # calculate the mean of values that are greater than zero. This is close
        # to 1000 times the "MEAN GAIN VALUE" reported in INIT.LP
        dval1d = modified.as_1d().as_double()
        mval = flex.mean(dval1d.select(dval1d > 0))
        new_val = int(mval)

    # Set this value everywhere that is not a masked value marked by -3
    print("Setting flat gain of {0}".format(new_val / 1000.0))
    modified.set_selected(modified >= 0, new_val)

    # Compress the data
    compressed = compress(modified)
    nbytes = len(compressed)

    # Update the header
    pre, post = cbf_header.split(xbsize_record)
    new_xbsize_record = "X-Binary-Size:{0:10d}".format(nbytes)
    if xbsize_record.endswith("\r"):
        new_xbsize_record += "\r"
    new_cbf_header = pre + new_xbsize_record + post

    # Write the file out.
    open(out_name, "wb").write(new_cbf_header + start_tag + compressed + tail)


if __name__ == "__main__":

    usage = """Usage:
  * To calculate the mean value and set that:
      squishGainCBF.py /path/to/GAIN.cbf

  * To force the flat gain value to e.g. 2.5:
      squishGainCBF.py /path/to/GAIN.cbf 2.5
"""

    if len(sys.argv) < 2:
        sys.exit(usage)
    if len(sys.argv) > 2:
        force_gain = float(sys.argv[2])
    else:
        force_gain = None

    with open(sys.argv[1], "rb") as cbf:
        squishGain(cbf, "squishedGAIN.cbf", force_gain)
