from __future__ import division
from __future__ import print_function

import os

from cctbx import uctbx
import iotbx.phil
from dials.util.options import OptionParser

help_message = """
"""

phil_scope = iotbx.phil.parse(
    """
space_group = None
  .type = space_group
n_bins = 20
  .type = int(value_min=1)
anomalous = False
  .type = bool
use_internal_variance = False
  .type = bool
eliminate_sys_absent = False
  .type = bool
size_inches = None
  .type = floats(size=2, value_min=0)
image_dir = None
  .type = path
labels = None
  .type = str
prefix = None
  .type = str
""",
    process_includes=True,
)


def run(args):
    import matplotlib

    matplotlib.use("Agg")
    import libtbx.load_env

    usage = "%s [options]" % libtbx.env.dispatcher_name

    parser = OptionParser(
        usage=usage, phil=phil_scope, check_format=False, epilog=help_message
    )

    params, options, args = parser.parse_args(
        show_diff_phil=True, return_unhandled=True
    )

    for mtz in args:
        print(mtz)
        assert os.path.isfile(mtz), mtz
        import iotbx.merging_statistics

        i_obs = iotbx.merging_statistics.select_data(mtz, data_labels=params.labels)
        if params.space_group is not None:
            i_obs = i_obs.customized_copy(space_group_info=params.space_group)

        from scitbx.array_family import flex

        # set the sigmas to 1, and calculate the mean intensities and internal variances
        intensities_copy = i_obs.customized_copy(sigmas=flex.double(i_obs.size(), 1))
        merging_internal = intensities_copy.merge_equivalents(
            use_internal_variance=True
        )
        merged = merging_internal.array()

        merging_external = i_obs.merge_equivalents(use_internal_variance=False)

        sigmas_internal = merging_internal.array().sigmas()
        sigmas_external = merging_external.array().sigmas()

        variances_internal = flex.pow2(sigmas_internal)
        variances_external = flex.pow2(sigmas_external)

        n_bins = 100
        i_obs.setup_binner_counting_sorted(n_bins=n_bins)
        sigmas_ratio = sigmas_external / sigmas_internal
        variance_ratio = variances_external / variances_internal

        array_sr = merging_external.array().customized_copy(
            data=sigmas_ratio, sigmas=None
        )
        array_sr.use_binning_of(i_obs)
        mean_sr = array_sr.mean(use_binning=True)

        ds2 = mean_sr.binner.bin_centers(2)
        sr = mean_sr.data[1:-1]

        array_vr = merging_external.array().customized_copy(
            data=variance_ratio, sigmas=None
        )
        array_vr.use_binning_of(i_obs)
        mean_vr = array_vr.mean(use_binning=True)

        d_star_sq = mean_vr.binner.bin_centers(2)
        vr = mean_vr.data[1:-1]

        prefix = params.prefix
        if prefix is None:
            prefix = ""

        from matplotlib import pyplot

        pyplot.style.use("ggplot")

        pyplot.plot(d_star_sq, sr)
        ax = pyplot.gca()
        xticks = ax.get_xticks()
        xticks_d = [
            "%.2f" % uctbx.d_star_sq_as_d(ds2) if ds2 > 0 else 0 for ds2 in xticks
        ]
        ax.set_xticklabels(xticks_d)
        pyplot.xlabel("d spacing (A)")
        pyplot.ylabel("<sigI_ext/sigI_int>")
        pyplot.savefig("%ssigmas_ratio.png" % prefix)
        pyplot.clf()

        pyplot.plot(d_star_sq, vr)
        ax = pyplot.gca()
        xticks = ax.get_xticks()
        xticks_d = [
            "%.2f" % uctbx.d_star_sq_as_d(ds2) if ds2 > 0 else 0 for ds2 in xticks
        ]
        ax.set_xticklabels(xticks_d)
        pyplot.xlabel("d spacing (A)")
        pyplot.ylabel("<varI_ext/varI_int>")
        pyplot.savefig("%svariances_ratio.png" % prefix)
        pyplot.clf()


if __name__ == "__main__":
    import sys
    from libtbx.utils import show_times_at_exit

    show_times_at_exit()
    run(sys.argv[1:])
