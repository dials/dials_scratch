#!/usr/bin/env python3
import json
import logging
import pathlib
import sys
import urllib.parse
import uuid
from argparse import ArgumentParser
from queue import Queue
from typing import Dict, List, Union

import workflows.recipe
from minio import Minio
from tqdm import tqdm
from workflows.transport.stomp_transport import StompTransport

logger = logging.getLogger(__name__)


UNIQUE_ID = str(uuid.uuid4())
PIA_UNIQUE_QUEUE = "transient.pia_feedback_{}".format(UNIQUE_ID)


def get_client(minio_server: str) -> Minio:
    """Parse a server URL and return a minio client"""
    parts = urllib.parse.urlparse(minio_server)

    host = f"{parts.hostname}:{parts.port}"
    user = f"{parts.username}@" if parts.username else ""
    logger.info(f"Using {user}{host}")
    return Minio(
        host,
        access_key=parts.username or None,
        secret_key=parts.password or None,
        secure=False,
    )


def ensure_url(url: str) -> str:
    """Make sure a string URL has a schema, for urllib.parse.urlparse consumption"""
    if "://" not in url:
        return f"minio://{args.host}"
    return url


logger = logging.getLogger()

parser = ArgumentParser(description="Submit an S3 bucket for PIA")

StompTransport.load_configuration_file(
    "/dls_sw/apps/zocalo/secrets/credentials-testing.cfg"
)
StompTransport.add_command_line_options(parser)

parser.add_argument(
    "s3_url",
    metavar="S3_URL",
    help="The access URL for the S3 bucket",
    type=ensure_url,
)
# parser.add_argument("images", metavar="IMAGES", help="Image numbers to submit. Index or ranges '1,10' '1-10'. Defauts to all")
parser.add_argument(
    "-v", "--verbose", help="increase output verbosity", action="store_true"
)

args = parser.parse_args()
logging.basicConfig(
    format="%(message)s", level=logging.DEBUG if args.verbose else logging.INFO
)
logging.getLogger("stomp.py").setLevel(logging.WARNING)

# Make sure we were passed a bucket
url = urllib.parse.urlparse(args.s3_url)
if not url.path:
    parser.print_usage()
    sys.exit("Error: Must specify bucket as path s3://host/bucket")

# Extract the bucket name from the URL
path = pathlib.PurePosixPath(url.path)
args.bucket = path.parts[1]
logger.info(f"Bucket: {args.bucket}")

client = get_client(args.s3_url)
objects = [
    x.object_name for x in client.list_objects(args.bucket) if x.object_name.isdigit()
]

stomp = StompTransport()
try:
    stomp.connect()

    # for image in tqdm(sorted(objects), desc=f"Submitting {len(objects)} images..."):
    # image_urls = [
    #     urllib.parse.urlunparse(url._replace(path=str(path / image)))
    #     for image in objects
    # ]
    def submit_images(transport, s3_url):
        recipe = {
            "1": {
                "service": "DLS Cephwatcher",
                "queue": "cephwatcher",
                "output": {"every": 2},
                "parameters": {"s3_url": args.s3_url,},
            },
            "2": {
                "service": "DLS Per-Image-Analysis",
                "queue": "per_image_analysis",
                "output": 3,
            },
            "3": {
                "service": "Return results Queue",
                "queue": PIA_UNIQUE_QUEUE,
                "id": UNIQUE_ID,
            },
            "start": [[1, {}]],
        }
        message = {"custom_recipe": recipe, "parameters": {}, "recipes": []}
        logger.info(f"Message length: {len(json.dumps(recipe))}")
        transport.send("processing_recipe", message)

    submit_images(stomp, args.s3_url)
    # print(message)
    # logger.info(f"Submitted {len(objects)} images.")

    q: "Queue[Union[Dict, List]]" = Queue()
    # Now wait for the reply messages
    logger.info(f"Waiting for image-analysis messages on {PIA_UNIQUE_QUEUE}")
    workflows.recipe.wrap_subscribe(
        stomp,
        PIA_UNIQUE_QUEUE,
        lambda _, _2, message: q.put(message),
        acknowledgement=False,
    )

    with tqdm(total=len(objects)) as progress:
        # def receive(rw, header, message):
        #     """Process a hitfinding message"""
        #     q.put(message)

        # Wait for results to filter in
        result_count = 0
        while result_count < len(objects):
            message = q.get()
            result_count += 1
            progress.write(
                f"Got results for image {message['file-number']} = {message['n_spots_total']}"
            )
            progress.update()


finally:
    stomp.disconnect()
