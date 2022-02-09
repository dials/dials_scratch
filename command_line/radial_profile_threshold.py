# LIBTBX_SET_DISPATCHER_NAME dials.radial_profile_threshold
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1


from __future__ import annotations

import logging
import math

import iotbx.phil
from libtbx.phil import parse
from scitbx import matrix

import dials.util.masking
from dials.algorithms.spot_finding.factory import SpotFinderFactory
from dials.algorithms.spot_finding.factory import phil_scope as spot_phil
from dials.array_family import flex
from dials.util import Sorry, log, show_mail_handle_errors, tabulate
from dials.util.options import ArgumentParser, flatten_experiments
from dials.util.version import dials_version

logger = logging.getLogger("dials")

help_message = """

Examples::

  dials.radial_profile_threshold image_*.cbf

  dials.radial_profile_threshold imported.expt
"""

phil_scope = iotbx.phil.parse(
    """\
n_bins = 100
  .type = int

image = 0
  .type = int
  .help = "Image on which to perform the analysis"

n_sigma = 3
  .type = int
  .help = "Sigma multiplier for determining the threshold value"

masking {
  include scope dials.util.masking.phil_scope
}

output {
    plot = None
      .type = path
      .help = "Save background plot to file"
    size_inches = None
      .type = floats(value_min=0, size=2)
    log = dials.radial_profile_threshold.log
      .type = str
}

""",
    process_includes=True,
)


@show_mail_handle_errors()
def run(args=None):
    usage = "dials.radial_profile_threshold [options] image_*.cbf"

    parser = ArgumentParser(
        usage=usage,
        phil=phil_scope,
        read_experiments=True,
        read_experiments_from_images=True,
        epilog=help_message,
    )

    params, options = parser.parse_args(args, show_diff_phil=True)

    # Ensure we have either a data block or an experiment list
    experiments = flatten_experiments(params.input.experiments)
    imagesets = experiments.imagesets()
    if len(imagesets) != 1:
        raise Sorry("Multiple imagesets are not supported")

    # Configure the logging
    log.config(logfile=params.output.log)
    logger.info(dials_version())

    if params.output.plot:
        import matplotlib

        matplotlib.use("agg")

        from matplotlib import pyplot

        fig = pyplot.figure(figsize=params.output.size_inches)
        ax = fig.add_subplot(111)

    imageset = imagesets[0]
    first, last = imageset.get_scan().get_image_range()
    if params.image < first or params.image > last:
        raise Sorry("Image outside of scan range")
    image = params.image

    logger.info(f"For image {image}:")
    peaks = radial_profile_threshold(
        imageset,
        image - first,  # indices passed to imageset.get_raw_data start from zero
        n_bins=params.n_bins,
        mask_params=params.masking,
        show_summary=True,
        n_sigma=params.n_sigma,
    )
    logger.info(f"Found {peaks.count(True)} strong pixels")

    if params.output.plot:
        pyplot.imshow(~peaks.as_numpy_array(), cmap=pyplot.cm.gray)
        pyplot.savefig(params.output.plot)


def radial_profile_threshold(
    imageset, indx, n_bins, mask_params=None, show_summary=True, n_sigma=3
):
    if mask_params is None:
        # Default mask params for trusted range
        mask_params = phil_scope.fetch(parse("")).extract().masking

    detector = imageset.get_detector()
    beam = imageset.get_beam()

    # Only working with single panel detector for now
    assert len(detector) == 1
    panel = detector[0]
    imageset_mask = imageset.get_mask(indx)[0]
    mask = dials.util.masking.generate_mask(imageset, mask_params)[0]
    mask = imageset_mask & mask

    n = matrix.col(panel.get_normal()).normalize()
    b = matrix.col(beam.get_s0()).normalize()
    wavelength = beam.get_wavelength()

    if math.fabs(b.dot(n)) < 0.95:
        raise Sorry("Detector not perpendicular to beam")

    # Get the corrected data
    data = imageset.get_corrected_data(indx)
    assert len(data) == 1
    data = data[0].as_double()

    spot_params = spot_phil.fetch(source=parse("")).extract()
    threshold_function = SpotFinderFactory.configure_threshold(spot_params)
    peak_pixels = threshold_function.compute_threshold(data, mask)
    signal = data.select(peak_pixels.iselection())
    background_pixels = mask & ~peak_pixels
    background = data.select(background_pixels.iselection())

    # print some summary information
    if show_summary:
        logger.info(f"Mean background: {flex.sum(background) / background.size():.3f}")
        if len(signal) > 0:
            logger.info(
                f"Max/total signal pixels: {flex.max(signal):.0f} / {flex.sum(signal):.0f}"
            )
        else:
            logger.info("No signal pixels on this image")
        logger.info(
            "Peak/background/masked pixels: %d / %d / %d"
            % (peak_pixels.count(True), background.size(), mask.count(False))
        )

    # compute histogram of two-theta values, then same weighted
    # by pixel values, finally divide latter by former to get
    # the radial profile out, need to set the number of bins
    # sensibly; inspired by method in PyFAI

    full_two_theta_array = panel.get_two_theta_array(beam.get_s0())
    two_theta_array = full_two_theta_array.as_1d().select(
        background_pixels.iselection()
    )

    # Use flex.weighted_histogram
    h0 = flex.weighted_histogram(two_theta_array, n_slots=n_bins)
    h1 = flex.weighted_histogram(two_theta_array, background, n_slots=n_bins)
    h2 = flex.weighted_histogram(
        two_theta_array, background * background, n_slots=n_bins
    )

    d0 = h0.slots()
    d1 = h1.slots()
    d2 = h2.slots()

    I = d1 / d0
    I2 = d2 / d0
    sig = flex.sqrt(I2 - flex.pow2(I))

    tt = h0.slot_centers()
    d_spacings = wavelength / (2.0 * flex.sin(0.5 * tt))

    # Determine the threshold value for each bin
    threshold = I + n_sigma * sig

    # Shift the full 2θ array to the lower bound and truncate
    infos = list(h0.slot_infos())
    lower_bound = infos[0].low_cutoff
    lookup = full_two_theta_array - lower_bound
    lookup.set_selected(lookup < 0, 0)

    # Truncate just under the shifted upper bound and rescale
    upper_bound = infos[-1].high_cutoff - lower_bound
    lookup.set_selected(lookup >= upper_bound, upper_bound - 1e-10)
    lookup /= upper_bound  # values now in range [0,1)

    # Convert to a size_t lookup into the threshold array
    lookup *= n_bins
    lookup = flex.floor(lookup).iround().as_size_t()

    # Now construct a threshold image
    thresh_im = threshold.select(lookup.as_1d())
    thresh_im.reshape(data.accessor())

    peaks = data > thresh_im
    peaks.reshape(data.accessor())

    if show_summary:
        logger.info("Threshold values at resolution bin centres")
        header = ["d min (Å)", "Threshold"]
        rows = []
        for d, t in zip(d_spacings, threshold):
            rows.append([f"{d:.4f}", f"{t:.4f}"])
        logger.info(tabulate(rows, header))

    return peaks


if __name__ == "__main__":
    run()
