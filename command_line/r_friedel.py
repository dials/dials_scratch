# LIBTBX_SET_DISPATCHER_NAME dev.dials.r_friedel
"""
Calculate R_{Friedel} on both structure factors (as per the original definition)
and on intensities. Definition as given in Glaeser & Downing (1993)
https://www.doi.org/10.1016/0304-3991(93)90064-5

    .. math::
      R_{Friedel} = \\dfrac{\\sum_{hkl}{|F_{hkl} - F_{-h,-k,-l}|}}{\\sum_{hkl}{\\left \\langle F_{hkl} \\right \\rangle}}

Usage: dev.dials.r_friedel hklin=scaled.mtz
"""

from __future__ import absolute_import, division, print_function

import logging
from iotbx import mtz
from cctbx import sgtbx
from scitbx.array_family import flex
import dials.util.log

from dials.util import Sorry, show_mail_handle_errors
from dials.util.options import OptionParser
from dials.algorithms.merging.merge import truncate

# Define a logger. __name__ cannot be used as this script is called directly.
# Omit the dev prefix to use the DIALS logger
logger = logging.getLogger("dials.r_friedel")


def merge_in_P1(intensities):
    return (
        intensities.customized_copy(
            space_group_info=sgtbx.space_group_info("P1"), anomalous_flag=True
        )
        .merge_equivalents()
        .array()
        .customized_copy(
            space_group_info=intensities.space_group_info(), anomalous_flag=True
        )
    )


def r_friedel(data_array):
    d_ano = data_array.anomalous_differences()
    friedel_mean = data_array.average_bijvoet_mates().common_set(other=d_ano)
    abs_anom_diff = flex.abs(d_ano.data())
    assert friedel_mean.size() == abs_anom_diff.size()
    logger.debug("R_{Friedel} using " + f"{friedel_mean.size()} pairs")
    numerator = flex.sum(abs_anom_diff)
    denominator = flex.sum(flex.abs(friedel_mean.data()))
    assert denominator > 0
    return numerator / denominator


class Script(object):
    """A class for running the script."""

    def __init__(self):
        """Initialise the script."""
        from libtbx.phil import parse

        # The phil scope
        phil_scope = parse(
            """
            output {
                log = dev.dials.r_friedel.log
                    .type = path
            }

            hklin = None
                .type = path
                .help = "MTZ file (containing observed and calculated structure "
                        "factors)"

            Io = I
                .type = str
                .help = "MTZ column name for Iobs"
            """,
            process_includes=True,
        )

        # The script usage
        usage = "usage: dev.dials.r_friedel hklin=scaled.mtz"

        # Create the parser
        self.parser = OptionParser(usage=usage, phil=phil_scope, epilog=__doc__)

        return

    def _extract_data_from_mtz(self):
        try:
            m = mtz.object(self.params.hklin)
        except RuntimeError:
            raise Sorry("Could not read {0}".format(self.params.hklin))

        mad = m.as_miller_arrays_dict(merge_equivalents=False)
        mad = {k[-1]: v for (k, v) in mad.items()}
        iobs = mad.get(self.params.Io)

        original_indices = m.extract_original_index_miller_indices()
        iobs = iobs.customized_copy(indices=original_indices)

        if iobs is None:
            raise Sorry(
                "Intensity column {0} not found in available labels: {1}".format(
                    self.params.Io,
                    ", ".join(m.column_labels()),
                )
            )

        if not iobs.is_unmerged_intensity_array():
            raise Sorry("Please check that the file contains unmerged intensities")

        return iobs

    def run(self, args=None):
        """Execute the script."""

        # Parse the command line
        self.params, options = self.parser.parse_args(args, show_diff_phil=True)

        # Configure the logging
        dials.util.log.config(verbosity=options.verbose, logfile=self.params.output.log)

        if self.params.hklin is None:
            self.parser.print_help()
            sys.exit()

        iobs = self._extract_data_from_mtz()

        i_p1 = merge_in_P1(iobs)
        f_p1 = truncate(i_p1)[1]
        fsq_p1 = f_p1.customized_copy(data=flex.pow2(f_p1.data()))

        logger.info("R_friedel(F) = {0:.5f}".format(r_friedel(f_p1)))
        logger.info("R_friedel(F^2) = {0:.5f}".format(r_friedel(fsq_p1)))
        logger.info("R_friedel(I) = {0:.5f}".format(r_friedel(i_p1)))

        return


@show_mail_handle_errors()
def run(args=None):
    script = Script()
    script.run(args)


if __name__ == "__main__":
    run()
