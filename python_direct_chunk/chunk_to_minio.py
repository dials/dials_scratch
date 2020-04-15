import minio
import io
import h5py
import sys
import time

master = sys.argv[1]
meta = sys.argv[1].replace("_master.h5", "_meta.h5")

config = ""
with h5py.File(meta, "r") as f:
    config = f["config"][()]
    print(config)

try:
    dcid = int(sys.argv[2])
except:
    dcid = 12345

key = "hello"
pwd = "k1tty-pass"

# create a connection to minio service
client = minio.Minio("localhost:9000", access_key=key, secret_key=pwd, secure=False)

# create a bucket for this DCID - if it already exists, will throw exception
bucket = "%d" % dcid
client.make_bucket(bucket, location="right-here-1")

etag = client.put_object(bucket, "header", io.BytesIO(config), len(config))

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
            name = "%06d" % idx

            total_read += chunk_size

            # push chunk to object store
            etag = client.put_object(bucket, name, io.BytesIO(chunk), chunk_size)

    t1 = time.time()

print("{0} GB read in {1}s".format(total_read / (1024.0 ** 3), t1 - t0))

# now read back every object from the bucket
total_read = 0
t0 = time.time()
for o in client.list_objects(bucket):
    back = client.get_object(o.bucket_name, o.object_name).read()
    total_read += len(back)
t1 = time.time()

print("{0} GB read back in {1}s".format(total_read / (1024.0 ** 3), t1 - t0))

# now clean up
# for o in client.list_objects(bucket):
#    client.remove_object(o.bucket_name, o.object_name)

# remove the bucket
# client.remove_bucket("%d" % dcid)
