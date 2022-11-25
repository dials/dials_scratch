# LIBTBX_SET_DISPATCHER_NAME dev.dials.find_singla_direct_beam

from __future__ import annotations
from dials.array_family import flex
import libtbx
from dials.util import Sorry, show_mail_handle_errors
from dials.util.options import ArgumentParser, flatten_experiments

from dials.algorithms.image.filter import convolve
import numpy as np

# from matplotlib import pyplot
# fig = pyplot.figure()
# ax = fig.add_subplot(111)

help_message = """

Examples::

  dev.dials.find_singla_direct_beam image.nxs

  dev.dials.find_singla_direct_beam imported.expt
"""

phil_scope = libtbx.phil.parse(
    """\
blur = narrow wide
    .type = choice
    .help = "Optional preprocessing of the image by a convolution with"
            "a simple Gaussian kernel of size either 3×3 (narrow) or"
            "5×5 (wide)."

""",
    process_includes=True,
)


def find_beam_centre(image):
    """Find the centre of gravity of the maximum pixels"""

    max_val = flex.max(image)
    positions = (image == max_val).iselection()
    ncol, nrow = image.all()
    x, y = zip(*[(pos % ncol, pos // ncol) for pos in positions])
    mean_x = sum(x) / len(x)
    mean_y = sum(y) / len(y)

    return mean_x, mean_y


@show_mail_handle_errors()
def run(args=None):
    usage = "dev.dials.find_singla_direct_beam image.nxs"

    parser = ArgumentParser(
        usage=usage,
        phil=phil_scope,
        read_experiments=True,
        read_experiments_from_images=True,
        epilog=help_message,
    )

    params, options = parser.parse_args(args, show_diff_phil=False)

    # Ensure we have either a data block or an experiment list
    experiments = flatten_experiments(params.input.experiments)
    imagesets = experiments.imagesets()
    if len(imagesets) != 1:
        raise Sorry("Multiple imagesets are not supported")

    # Set approximate Gaussian kernel for blurring
    if params.blur == "narrow":
        # fmt: off
        kernel = flex.double(
               (0.0625, 0.125, 0.0625,
                0.125,  0.25,  0.125,
                0.0625, 0.125, 0.0625)
        )
        # fmt: on
        kernel.reshape(flex.grid((3, 3)))
    elif params.blur == "wide":
        # fmt: off
        kernel = (
            flex.double(
                (
                    1,  4,  7,  4,  1,
                    4, 16, 26, 16,  4,
                    7, 26, 41, 26,  7,
                    4, 16, 26, 16,  4,
                    1,  4,  7,  4,  1,
                )
            ) / 273
        )
        # fmt: on
        kernel.reshape(flex.grid((5, 5)))
    else:
        kernel = None

    # Set the ROI to be +/- 100 pixels around the image centre
    xc, yc = (e // 2 for e in experiments[0].detector[0].get_image_size())
    x0 = xc - 100
    x1 = xc + 100
    y0 = yc - 100
    y1 = yc + 100

    saturation = experiments[0].detector[0].get_trusted_range()[1]

    imageset = imagesets[0]
    first, last = imageset.get_scan().get_image_range()
    step = imageset.get_scan().get_num_images() // 10
    images = []
    for i in range(0, last - first, step):
        image = imageset.get_raw_data(i)[0].as_double()
        mask = imageset.get_mask(i)[0]

        image = image[y0:y1, x0:x1]
        mask = mask[y0:y1, x0:x1]

        # Set masked pixels to zero to remove the module gap
        image = image.set_selected(~mask, 0)

        # Set overloaded pixels to the saturation value
        image = image.set_selected(image == -1, saturation)

        # Gaussian blur
        if kernel:
            image = convolve(image, kernel)

        images.append(image)

        # pyplot.imshow((image).as_numpy_array(), cmap=pyplot.cm.gray)
        # pyplot.show()

    beam_centres = [find_beam_centre(im) for im in images]
    x, y = zip(*beam_centres)

    # Need to be robust against blank images. Remove any value more than 5 px
    # from the median
    med_x = np.median(x)
    med_y = np.median(y)
    x = [e for e in x if abs(e - med_x) < 5]
    y = [e for e in y if abs(e - med_y) < 5]

    mean_x = sum(x) / len(x)
    mean_y = sum(y) / len(y)

    # Offset of the ROI
    fast = xc - 100 + mean_x
    slow = yc - 100 + mean_y

    print(f"fast_slow_beam_centre={fast,slow}")


if __name__ == "__main__":
    run()
