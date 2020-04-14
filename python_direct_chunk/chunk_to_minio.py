import minio
import io
import h5py
import sys
import time

master = sys.argv[1]
try:
    dcid = int(sys.argv[2])
except:
    dcid = 12345

# create a connection to minio service
client = minio.Minio(
    "localhost:9000", access_key="minioadmin", secret_key="minioadmin", secure=False
)

# create a bucket for this DCID - if it already exists, will throw exception
client.make_bucket("%d" % dcid, location="right-here-1")

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
            chunk_size = len(chunk)

            # increment counter - we push this in people numbers
            idx += 1
            total_read += chunk_size

            # push chunk to object store
            etag = client.put_object(
                "%d" % dcid, "%d" % idx, io.BytesIO(chunk), chunk_size
            )
            print(etag, chunk_size)

    t1 = time.time()

print("{0} GB read in {1}s".format(total_read / (1024.0 ** 3), t1 - t0))
