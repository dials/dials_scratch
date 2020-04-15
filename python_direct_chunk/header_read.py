import minio
import sys
import ast
import numpy
import bitshuffle
from dials.array_family import flex

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

header_text = client.get_object(bucket, "header").read()
header = ast.literal_eval(header_text.decode("utf-8"))

print(header)

ny = header["y_pixels_in_detector"]
nx = header["x_pixels_in_detector"]

if header["bit_depth_image"] == 32:
    dtype = numpy.dtype("uint32")
elif header["bit_depth_image"] == 16:
    dtype = numpy.dtype("uint16")
elif header["bit_depth_image"] == 8:
    dtype = numpy.dtype("uint8")
else:
    raise RuntimeError("bit depth %d unknown" % header["bit_depth_image"])


def find_spots(image):
    from dials.algorithms.spot_finding.threshold import DispersionThresholdStrategy
    from dials.model.data import PixelList
    from dials.model.data import PixelListLabeller

    thresholder = DispersionThresholdStrategy(gain=1)
    mask = image.as_1d() >= 0  # flex.bool(image.size(), True)
    mask.reshape(flex.grid(*image.focus()))

    threshold_mask = thresholder(image, mask=mask)
    plist = PixelList(0, image, threshold_mask)

    pixel_labeller = PixelListLabeller()
    pixel_labeller.add(plist)

    creator = flex.PixelListShoeboxCreator(pixel_labeller, 0, 0, True, 2, 100, False)
    shoeboxes = creator.result()
    return shoeboxes


for o in client.list_objects(bucket):
    if o.object_name is "header":
        continue
    image_number = int(o.object_name)
    chunk = client.get_object(o.bucket_name, o.object_name).read()
    blob = numpy.fromstring(chunk[12:], dtype=numpy.uint8)
    image = bitshuffle.decompress_lz4(blob, (ny, nx), dtype)
    image = flex.int(image.astype("int32"))
    shoeboxes = find_spots(image)
    print(image_number, flex.max(image), len(shoeboxes))
