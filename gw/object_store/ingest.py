import io
import time
import urllib.parse
from argparse import ArgumentParser

import h5py
import minio
from tqdm import tqdm


class Ingester(object):
    def __init__(self, minio_server):
        parts = urllib.parse.urlparse(minio_server)

        host = f"{parts.hostname}:{parts.port}"
        user = f"{parts.username}@" if parts.username else ""
        print(f"Using {user}{host}")
        self._client = minio.Minio(
            host,
            access_key=parts.username or None,
            secret_key=parts.password or None,
            secure=False,
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

            # Add up length of chunks
            n_image_chunks = sum(
                data[d].shape[0] for d in data if d.startswith("data_")
            )

            with tqdm(total=n_image_chunks) as progress:
                # assume all the "real" data in data_(something) datasets
                for k in sorted(d for d in data if d.startswith("data_")):
                    dataset = data[k]
                    d_id = dataset.id

                    # getting the dataset shape should
                    chunks = dataset.shape[0]

                    # iterate through the chunks found in here and push each up
                    for j in range(chunks):
                        offset = (j, 0, 0)
                        filter_mask, chunk = d_id.read_direct_chunk(offset)
                        chunk_size = len(chunk)

                        # increment counter - we push this in people numbers
                        frame_index += 1

                        # push chunk to object store
                        self._client.put_object(
                            bucket, "%06d" % frame_index, io.BytesIO(chunk), chunk_size
                        )
                        total_read += chunk_size
                        progress.update(1)

        t1 = time.time()
        return total_read, t1 - t0


def ingest(host, master, dcid):
    ingester = Ingester(host)
    r, t = ingester.ingest(master, master.replace("_master.h5", "_meta.h5"), str(dcid))

    print("{0} GB read in {1}s".format(r / (1024.0 ** 3), t))


if __name__ == "__main__":
    parser = ArgumentParser(description="Loads h5 images to object store")
    parser.add_argument("filename", metavar="MASTER_FILE")
    parser.add_argument("dcid", metavar="DCID", type=int)
    parser.add_argument(
        "--host",
        metavar="URL",
        type=str,
        help="Endpoint; user:pass@host:port",
        default="minio://hello:k1tty-pass@localhost:9000",
    )
    args = parser.parse_args()
    # Try to be default-sensible
    if not "://" in args.host:
        args.host = f"minio://{args.host}"
    ingest(args.host, args.filename, args.dcid)
