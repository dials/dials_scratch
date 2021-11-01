# LIBTBX_SET_DISPATCHER_NAME dev.dials.convert_to_XDS_frame
"""
Convert an imported.expt to the XDS coordinate system by comparison with an
XDS.INP by aligning the detector axes.

Usage: dev.dials.convert_to_XDS_frame xds_inp=XDS.INP
"""

import logging
import sys
import numpy as np
from scitbx import matrix
import libtbx.phil
from dxtbx.model import ExperimentList
import dials.util
import dials.util.log
from dials.util.options import OptionParser, flatten_experiments
from scitbx.math import r3_rotation_axis_and_angle_from_matrix
from dials.util.version import dials_version

try:
    from typing import List
except ImportError:
    pass

from rstbx.cftbx.coordinate_frame_helpers import align_reference_frame
import iotbx.xds

# Define a logger. __name__ cannot be used as this script is called directly.
# Omit the dev prefix to use the DIALS logger
logger = logging.getLogger("dials.convert_to_XDS_frame")

# Define the master PHIL scope for this program.
phil_scope = libtbx.phil.parse(
    """
    input {
        xds_inp = None
            .type = path
            .help = "The input XDS.INP file."
    }
    output {
        experiments = converted.expt
            .type = path
            .help = "File name for the converted DIALS experiments file."
        log = dev.dials.convert_to_XDS_frame.log
            .type = path
    }
    """
)


def read_xds_inp(xds_inp_file):
    """Extract the detector axes from an XDS.INP file"""
    handle = iotbx.xds.xds_inp.reader()
    handle.read_file(xds_inp_file)
    x, y = handle.direction_of_detector_x_axis, handle.direction_of_detector_y_axis

    x = matrix.col(x).normalize()
    y = matrix.col(y).normalize()

    return x, y


def align_experiments(
    experiments: ExperimentList,
    params: libtbx.phil.scope_extract,
) -> ExperimentList:

    if len(experiments) > 1:
        logger.info(
            "Only the first experiment will be used to determine the detector axes"
        )
    expt = experiments[0]
    detector = expt.detector
    if len(detector) > 1:
        logger.info("Only the first panel will be used to determine the detector axes")
    panel = detector[0]

    xds_x, xds_y = read_xds_inp(params.input.xds_inp)
    R = align_reference_frame(
        panel.get_fast_axis(), xds_x, panel.get_slow_axis(), xds_y
    )

    axis_angle = r3_rotation_axis_and_angle_from_matrix(R)
    axis = axis_angle.axis
    angle = axis_angle.angle()
    logger.info(
        f"Rotating experiment{'s' if len(experiments) else ''} about axis {axis} by {np.degrees(angle)}°"
    )

    for expt in experiments:
        expt.detector.rotate_around_origin(axis, angle, deg=False)
        expt.beam.rotate_around_origin(axis, angle, deg=False)
        expt.goniometer.rotate_around_origin(axis, angle, deg=False)

        if expt.crystal is not None:
            expt.crystal = rotate_crystal(expt.crystal, R, axis, angle)

    return experiments


def rotate_crystal(crystal, Rmat, axis, angle):

    Amats = []
    if crystal.num_scan_points > 0:
        scan_pts = list(range(crystal.num_scan_points))
        Amats = [
            Rmat
            * matrix.sqr(crystal.get_U_at_scan_point(t))
            * matrix.sqr(crystal.get_B_at_scan_point(t))
            for t in scan_pts
        ]

    crystal.rotate_around_origin(axis, angle, deg=False)
    if Amats:
        crystal.set_A_at_scan_points(Amats)

    return crystal


@dials.util.show_mail_on_error()
def run(args: List[str] = None, phil: libtbx.phil.scope = phil_scope) -> None:
    """
    Args:
        args: The arguments supplied by the user (default: sys.argv[1:])
        phil: The PHIL scope definition (default: phil_scope, the master PHIL scope
        for this program).
    """
    usage = "dev.dials.convert_to_XDS_frame [options] imported.expt xds_inp=XDS.INP"

    parser = OptionParser(
        usage=usage,
        phil=phil,
        read_experiments=True,
        check_format=False,
        epilog=__doc__,
    )

    params, options, args = parser.parse_args(
        args=args, show_diff_phil=False, return_unhandled=True
    )

    # Configure the logging.
    dials.util.log.config(options.verbose, logfile=params.output.log)

    # Log the dials version
    logger.info(dials_version())

    # Log the difference between the PHIL scope definition and the active PHIL scope,
    # which will include the parsed user inputs.
    diff_phil = parser.diff_phil.as_str()
    if diff_phil:
        logger.info("The following parameters have been modified:\n%s", diff_phil)

    experiments = flatten_experiments(params.input.experiments)

    if not params.input.xds_inp:
        if args:
            params.input.xds_inp = args[0]
        else:
            parser.print_help()
            return

    # Check the models and data
    if len(experiments) == 0:
        parser.print_help()
        return

    if not iotbx.xds.xds_inp.reader.is_xds_inp_file(params.input.xds_inp):
        sys.exit(f"Cannot interpret {params.input.xds_inp} as an XDS.INP file")

    # Do whatever this program is supposed to do.
    aligned_experiments = align_experiments(experiments, params)

    logger.info(f"Saving optimised experiments to {params.output.experiments}")
    aligned_experiments.as_file(params.output.experiments)


# Keep this minimal. Calling run() should do exactly the same thing as running this
if __name__ == "__main__":
    run()
