# LIBTBX_SET_DISPATCHER_NAME dev.dials.event_counter

from __future__ import absolute_import, division, print_function

import iotbx.phil
import libtbx.load_env
from scitbx.array_family import flex
from dials.util import Sorry

help_message = (
    """

Examples::

  %s data_master.h5

"""
    % libtbx.env.dispatcher_name
)

phil_scope = iotbx.phil.parse(
    """
images = None
  .type = ints
  .help = "Images on which to perform the analysis (otherwise use all images)"
image_range = None
  .type = ints
  .help = "Image range to consider"
"""
)


def run(args):

    from dials.util.options import OptionParser
    from dials.util.options import flatten_experiments

    usage = "%s [options] data_master.h5" % (libtbx.env.dispatcher_name)

    parser = OptionParser(
        usage=usage,
        phil=phil_scope,
        read_experiments=True,
        read_experiments_from_images=True,
        epilog=help_message,
    )

    params, options = parser.parse_args(show_diff_phil=True)

    experiments = flatten_experiments(params.input.experiments)
    if len(experiments) != 1:
        parser.print_help()
        print("Please pass an experiment list\n")
        return

    imagesets = experiments.imagesets()

    if len(imagesets) != 1:
        raise Sorry("Please pass an experiment list that contains one imageset")

    imageset = imagesets[0]

    first, last = imageset.get_scan().get_image_range()
    images = range(first, last + 1)

    if not params.images and params.image_range:
        params.images = list(range(params.image_range[0], params.image_range[1] + 1))

    if params.images:
        if min(params.images) < first or max(params.images) > last:
            raise Sorry("image outside of scan range")
        images = params.images

    detectors = imageset.get_detector()
    assert len(detectors) == 1
    detector = detectors[0]
    trusted = detector.get_trusted_range()

    # construct an integer array same shape as image; accumulate number of
    # "signal" pixels in each pixel across data

    total = None

    from dials.util.command_line import ProgressBar

    p = ProgressBar(title="Counting events")

    events_per_image = {}

    for idx in images:

        p.update(idx * 100.0 / len(images))

        pixels = imageset.get_raw_data(idx - 1)
        assert len(pixels) == 1
        mask = ~imageset.get_mask(idx)[0].as_1d()
        data = pixels[0].as_1d()

        data.set_selected(mask, 0)
        events_per_image[idx] = flex.sum(data)

    for idx in images:
        print(idx, events_per_image[idx])


if __name__ == "__main__":
    import sys

    run(sys.argv[1:])
