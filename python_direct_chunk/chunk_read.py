import h5py
import sys
import time

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
        for j in range(chunks):
            offset = (j, 0, 0)
            filter_mask, chunk = d_id.read_direct_chunk(offset)

            # increment counter - we push this in people numbers
            idx += 1
            total_read += len(chunk)

            # push chunk to object store

    t1 = time.time()

print("{0} GB read in {1}s".format(total_read / (1024.0 ** 3), t1 - t0))
