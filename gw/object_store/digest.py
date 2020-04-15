import minio
import sys
import ast
import numpy
import bitshuffle

from dials.array_family import flex
from dials.algorithms.spot_finding.threshold import DispersionThresholdStrategy
from dials.model.data import PixelList, PixelListLabeller


def image_to_shoeboxes(image):
    """For a given image, find spots 2 - 100 pixels im size, assuming a gain
    of 1, return the list of spot shoeboxes. Also assumes valid intensities in
    range 0...N."""

    thresholder = DispersionThresholdStrategy(gain=1)

    mask = image.as_1d() >= 0
    mask.reshape(flex.grid(*image.focus()))

    threshold_mask = thresholder(image, mask=mask)
    plist = PixelList(0, image, threshold_mask)

    pixel_labeller = PixelListLabeller()
    pixel_labeller.add(plist)

    creator = flex.PixelListShoeboxCreator(pixel_labeller, 0, 0, True, 2, 100, False)
    shoeboxes = creator.result()
    return shoeboxes


class Digester(object):
    def __init__(self, minio_server, username, password):
        self._client = minio.Minio(
            "localhost:9000", access_key=username, secret_key=password, secure=False
        )

    def chunk_to_image(self, chunk):
        """Decompress chunk to image."""

        blob = numpy.fromstring(chunk[12:], dtype=numpy.uint8)
        image = bitshuffle.decompress_lz4(blob, self._dim, self._dtype)
        image = flex.int(image.astype("int32"))
        return image

    def digest(self, bucket):
        """Recover the data from bucket, header first then find spots on all
        the image data found therein."""

        # recover header package & unpack

        header_text = self._client.get_object(bucket, "header").read()
        self._header = ast.literal_eval(header_text.decode("utf-8"))
        self._dtype = numpy.dtype("uint%d" % self._header["bit_depth_image"])
        self._dim = (
            self._header["y_pixels_in_detector"],
            self._header["x_pixels_in_detector"],
        )

        result = {}

        for item in self._client.list_objects(bucket):
            if item.object_name == "header":
                continue
            image_number = int(item.object_name)
            chunk = self._client.get_object(item.bucket_name, item.object_name).read()
            image = self.chunk_to_image(chunk)
            shoeboxes = image_to_shoeboxes(image)
            result[image_number] = len(shoeboxes)

        return result


def digest(dcid):
    digester = Digester("localhost:9000", "hello", "k1tty-pass")
    result = digester.digest(str(dcid))

    for image_number in sorted(result):
        print(image_number, result[image_number])


if __name__ == "__main__":
    dcid = int(sys.argv[1])
    digest(dcid)
