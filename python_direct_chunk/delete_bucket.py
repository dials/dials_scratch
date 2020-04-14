import minio
import sys

try:
    dcid = int(sys.argv[2])
except:
    dcid = 12345

# create a connection to minio service
client = minio.Minio(
    "localhost:9000", access_key="minioadmin", secret_key="minioadmin", secure=False
)

# remove everything from the bucket
for o in client.list_objects("%d" % dcid):
    print("Removing {0}/{1}".format(o.bucket_name, o.object_name))
    client.remove_object(o.bucket_name, o.object_name)

# remove the bucket
client.remove_bucket("%d" % dcid)
