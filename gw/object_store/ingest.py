import io
import time
from argparse import ArgumentParser

import h5py
import minio


class Ingester(object):
    def __init__(self, minio_server, username, password):
        self._client = minio.Minio(
            "localhost:9000", access_key=username, secret_key=password, secure=False
        )

    def ingest(self, master, meta, bucket):
        """Ingest the data from master h5, with the metadata pulled from meta h5
        (which is in the /config string dataset for Diamond HDF5 data) and
        move it into the minio server - the header goes in as "header" with the
        images numbered from 1 as "000001" etc."""

        # FIXME add error trapping etc.
        self._client.make_bucket(bucket)
        # if location is required: , location="right-here-1")

        with h5py.File(meta, "r") as f:
            config = f["config"][()]

        etag = self._client.put_object(
            bucket, "header", io.BytesIO(config), len(config)
        )

        # now iterate through the chunks in the actual data set - assumed to
        # be stored in external data files

        with h5py.File(master, "r") as f:
            frame_index = 0
            data = f["/entry/data"]
            total_read = 0
            t0 = time.time()

            # assume all the "real" data in data_(something) datasets
            for k in sorted(d for d in data if d.startswith("data_")):
                d = data[k]
                d_id = d.id

                # getting the dataset shape should
                chunks = d.shape[0]

                # iterate through the chunks found in here and push each up
                for j in range(chunks):
                    offset = (j, 0, 0)
                    filter_mask, chunk = d_id.read_direct_chunk(offset)
                    chunk_size = len(chunk)

                    # increment counter - we push this in people numbers
                    frame_index += 1

                    # push chunk to object store
                    etag = self._client.put_object(
                        bucket, "%06d" % frame_index, io.BytesIO(chunk), chunk_size
                    )
                    total_read += chunk_size
        t1 = time.time()
        return total_read, t1 - t0


def ingest(master, dcid):
    ingester = Ingester("localhost:9000", "hello", "k1tty-pass")
    r, t = ingester.ingest(master, master.replace("_master.h5", "_meta.h5"), str(dcid))

    print("{0} GB read in {1}s".format(r / (1024.0 ** 3), t))


if __name__ == "__main__":
    parser = ArgumentParser(description="Loads h5 images to object store")
    parser.add_argument("filename", metavar="MASTER_FILE")
    parser.add_argument("dcid", metavar="DCID", type=int)
    args = parser.parse_args()
    ingest(args.filename, args.dcid)
