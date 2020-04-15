import h5py
import sys
import time
import bitshuffle
import numpy
from dials.array_family import flex

master = sys.argv[1]
try:
    dcid = int(sys.argv[2])
except:
    dcid = 12345

with h5py.File(master, "r") as f:
    idx = 0
    data = f["entry/data"]
    total_read = 0
    t0 = time.time()
    for k in sorted(data):
        if not k.startswith("data_"):
            continue
        d = data[k]
        d_id = d.id
        chunks = d.shape[0]
        ny, nx = d.shape[1:3]
        dtype = d.dtype
        for j in range(chunks):
            offset = (j, 0, 0)
            filter_mask, chunk = d_id.read_direct_chunk(offset)

            # decompress
            blob = numpy.fromstring(chunk[12:], dtype=numpy.uint8)
            image = bitshuffle.decompress_lz4(blob, (ny, nx), dtype)

            image = flex.int(image.astype("int32"))
            print(j, flex.max(image))

            # only need to overwrite values if read and used in 16-bit mode
            if dtype == numpy.uint16:
                bad = 2 ** 16 - 1
                sel = image.as_1d() >= bad
                image.as_1d().set_selected(sel, -1)

            # increment counter - we push this in people numbers
            idx += 1
            total_read += len(chunk)

            # push chunk to object store

    t1 = time.time()

print("{0} GB read in {1}s".format(total_read / (1024.0 ** 3), t1 - t0))
