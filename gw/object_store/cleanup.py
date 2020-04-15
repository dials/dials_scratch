import minio
import sys


class Cleaner(object):
    def __init__(self, minio_server, username, password):
        self._client = minio.Minio(
            "localhost:9000", access_key=username, secret_key=password, secure=False
        )

    def cleanup(self, bucket):
        """Remove every item from the bucket before removing the bucket"""

        for item in self._client.list_objects(bucket):
            self._client.remove_object(item.bucket_name, item.object_name)

        self._client.remove_bucket(bucket)


def cleanup(dcid):
    cleaner = Cleaner("localhost:9000", "hello", "k1tty-pass")
    cleaner.cleanup(str(dcid))


if __name__ == "__main__":
    dcid = int(sys.argv[1])
    cleanup(dcid)
