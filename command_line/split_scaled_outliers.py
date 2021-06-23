# LIBTBX_SET_DISPATCHER_NAME dev.dials.split_scaled_outliers

"""
Separate scaled reflections from intensity outliers.  Save them in separate files.

The scaling component of the DIALS pipeline, dials.scale, rejects reflections that
have intensities that are deemed to be outliers according to a Wilson intensity
distribution.  They are not used to calculate scaling statistics, but are labelled
and retained in the output table of reflections.  This utility separates the outlier
reflections and the scaled reflections in the reflection table output from
dials.scale and saves them in separate reflection table files.
"""

import argparse
import pathlib
import sys
from typing import Union, List, Optional

from dials.array_family import flex

File = Union[str, pathlib.Path]


parser = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter
)
parser.add_argument(
    "reflection-table",
    help=(
        "A DIALS reflection table (.refl) file containing scaled reflections, "
        "as produced by dials.scale."
    ),
    type=pathlib.Path,
)
parser.add_argument(
    "--scaled-output",
    "-s",
    default="scaled.refl",
    help="File name for output table of scaled reflections.",
    type=pathlib.Path,
)
parser.add_argument(
    "--outliers-output",
    "-o",
    default="outliers.refl",
    help="File name for output table of reflections deemed to be intensity outliers.",
    type=pathlib.Path,
)
parser.add_argument(
    "--force-overwrite",
    "-f",
    help="Overwrite pre-existing output files.",
    action="store_true",
)


def split_scaled_and_outliers(
    reflection_table: File, scaled_output: File, outliers_output: File
) -> None:
    """
    Save scaled reflections and intensity outliers to separate reflection tables.

    Args:
        reflection_table:  A DIALS reflection table containing reflections after
                           scaling.
        scaled_output:  The file path for the scaled reflections.
        outliers_output:  The file path for the intensity outliers.
    """
    reflections = flex.reflection_table.from_file(reflection_table)

    scaled = reflections.get_flags(reflections.flags.scaled)
    outliers = reflections.get_flags(reflections.flags.outlier_in_scaling)

    reflections.select(scaled).as_file(scaled_output)
    reflections.select(outliers).as_file(outliers_output)


def command_line_interface(args: Optional[List[str]] = None):
    """
    Call the split_scaled_and_outliers routine as if from the command line.

    Args:
        args:  List of argument words passed as if from the command line input
               (defaults to sys.argv[1:]).
    """
    args = parser.parse_args(args)

    if not args.reflection_table.is_file():
        sys.exit(f"Could not find specified input file {args.reflection_table}.")

    split_scaled_and_outliers(
        args.reflection_table, args.scaled_output, args.outliers_output
    )


if __name__ == "__main__":
    command_line_interface()
