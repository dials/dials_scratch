#!/usr/bin/env cctbx.python

"""Load a miniCBF, multiply the counts, then write out again"""
from __future__ import print_function

import binascii
import copy
import sys
import os
from cbflib_adaptbx import compress, uncompress
from scitbx.array_family import flex


def multiplyCounts(cbf_file, out_name, multiplier):

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

    # multiply all positive values
    modified = copy.deepcopy(values).as_1d()
    sel = modified > 0
    new_val = modified.select(sel) * multiplier
    modified.set_selected(sel, new_val)

    # reshape
    modified.reshape(values.accessor())

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

  To double the counts in an image:
      scaleCBF.py /path/to/image.cbf 2
"""

    if len(sys.argv) < 3:
        sys.exit(usage)
    multiplier = int(sys.argv[2])
    print("set multiplier to {0}".format(multiplier))

    with open(sys.argv[1], "rb") as cbf:
        out_name = "mod_" + os.path.split(sys.argv[1])[1]
        multiplyCounts(cbf, out_name, multiplier=2)
