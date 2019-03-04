#!/usr/bin/env python
#
#  Copyright (C) 2016 Diamond Light Source and STFC Rutherford Appleton
#  Laboratory, UK.
#
#  Author: David Waterman
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

"""Work through a list of indexed.pickles and do refinement of a static model
and scan-varying models with various interval widths (i.e. smoothness). Perform
centroid analysis on each to produce periodograms of the refined residuals.
Save the plots to see if the varying model accounts for the main features of
the residuals.

This will be run on datasets from the Metrix database. The idea is to explore
datasets with many different features, to investigate in which situations it is
possible to automatically determine a sensible interval width for scan-varying
refinement.

Usage: dials.python analyse_periodograms.py filelist

where filelist consists of paths to indexed.pickles. Associated
experiments.jsons will be determined by matching expected filenames."""

from __future__ import division
from __future__ import print_function
import os, shutil
import datetime
import glob
from libtbx.utils import Sorry
from libtbx import easy_run
import libtbx.load_env
from libtbx.test_utils import open_tmp_directory
from dials.command_line.analyse_output import ensure_directory


class Script(object):
    """A class for running the script."""

    def __init__(self):
        """Initialise the script."""
        from dials.util.options import OptionParser
        from libtbx.phil import parse
        import libtbx.load_env

        # The phil scope
        phil_scope = parse(
            """
      force_analysis = False
        .type = bool
        .help = "Set True to force overwriting of previous analysis results."
                "Otherwise subdirectories for jobs that already exist will be"
                "skipped."

      remove_working_dir = True
        .type = bool
        .help = "Set False to keep working directories containing refinement"
                "results, otherwise these will be deleted."

      output {
        directory = analyse_periodograms
          .type = path
          .help = "The root directory in which to store results"
      }
    """,
            process_includes=True,
        )

        # The script usage
        import __main__

        usage = "usage: dials.python {0} [options] [param.phil] filelist ".format(
            __main__.__file__
        )

        # Create the parser
        self.parser = OptionParser(usage=usage, phil=phil_scope, epilog=__doc__)

        return

    def run(self):
        """Execute the script."""

        # Parse the command line
        self.params, options, filelist = self.parser.parse_args(
            show_diff_phil=True, return_unhandled=True
        )

        if len(filelist) != 1:
            raise Sorry(
                "Please provide input in the form:\n  {0}".format(self.parser.usage)
            )

        # Read the filelist and tidy it up
        try:
            with open(filelist[0], "r") as f:
                self.filelist = f.readlines()
        except IOError:
            raise Sorry("Cannot read the specified file list")
        self.filelist = [l.rstrip() for l in self.filelist]
        self.filelist = [l for l in self.filelist if l != ""]

        # Create output directory if it does not already exist
        self._directory = os.path.abspath(self.params.output.directory)
        ensure_directory(self._directory)
        print("Analysis will be performed in {0}".format(self._directory))

        # locate centroid_analysis script
        dials_dir = libtbx.env.find_in_repositories("dials")
        self._ca = os.path.join(
            dials_dir, "algorithms", "refinement", "analysis", "centroid_analysis.py"
        )

        # Do analysis for the files in filelist
        for i, f in enumerate(self.filelist):
            # strip newline from the filenames
            s = "{:%Y-%m-%d %H:%M:%S} ".format(datetime.datetime.now())
            s += "Processing file {0}: ".format(i) + f
            print(s)
            status = self.process(f, i)
            if status:
                print("Incomplete:", status)

        return

    def process(self, idx_path, i):
        """Perform processing tasks for one indexed.pickle"""

        cwd = os.path.abspath(os.curdir)
        self._working_dir = None
        status = self._process_core(idx_path, i)
        os.chdir(cwd)
        # clean up working directory
        if self._working_dir is not None and self.params.remove_working_dir:
            shutil.rmtree(self._working_dir)
        return status

    def _process_core(self, idx_path, i):

        directory = os.path.join(self._directory, "run%05d" % i)
        try:
            os.mkdir(directory)
        except OSError:
            if not self.params.force_analysis:
                return "Skipping " + directory + " that already exists"
        os.chdir(directory)

        # try to find an associated experiment.json
        exp_path = self.find_experiments_json(idx_path)
        if exp_path is None:
            return "Cannot find an experiments.json for refinement"

        # do refinement in a temporary directory
        tmp_dir = open_tmp_directory(suffix="_working_dir")
        self._working_dir = os.path.abspath(tmp_dir)
        os.chdir(tmp_dir)

        # refine static model first
        cmd = (
            "dials.refine {0} {1} "
            "output.experiments=refined_static.json "
            "output.reflections=refined_static.pickle"
        ).format(exp_path, idx_path)
        result = easy_run.fully_buffered(command=cmd)
        tst = [os.path.exists("refined_static" + e) for e in [".json", ".pickle"]]
        if tst.count(True) != 2:
            return "Static refinement output was not found"

        # refine scan-varying, 54 degrees interval width
        args = [
            "dials.refine",
            "refined_static.pickle",
            "refined_static.json",
            "scan_varying=true",
            "unit_cell.smoother.interval_width_degrees=54",
            "orientation.smoother.interval_width_degrees=54",
            "output.experiments=sv_refined_54deg.json",
            "output.reflections=sv_refined_54deg.pickle",
        ]
        cmd = " ".join(args)
        result = easy_run.fully_buffered(command=cmd)
        tst = [os.path.exists("sv_refined_54deg" + e) for e in [".json", ".pickle"]]
        if tst.count(True) != 2:
            return "54 deg interval width scan-varying refinement output was not found"

        # refine scan-varying, 36 degrees interval width
        args = [
            "dials.refine",
            "refined_static.pickle",
            "refined_static.json",
            "scan_varying=true",
            "unit_cell.smoother.interval_width_degrees=36",
            "orientation.smoother.interval_width_degrees=36",
            "output.experiments=sv_refined_36deg.json",
            "output.reflections=sv_refined_36deg.pickle",
        ]
        cmd = " ".join(args)
        result = easy_run.fully_buffered(command=cmd)
        tst = [os.path.exists("sv_refined_36deg" + e) for e in [".json", ".pickle"]]
        if tst.count(True) != 2:
            return "36 deg interval width scan-varying refinement output was not found"

        # refine scan-varying, 18 degrees interval width
        args = [
            "dials.refine",
            "refined_static.pickle",
            "refined_static.json",
            "scan_varying=true",
            "unit_cell.smoother.interval_width_degrees=18",
            "orientation.smoother.interval_width_degrees=18",
            "output.experiments=sv_refined_18deg.json",
            "output.reflections=sv_refined_18deg.pickle",
        ]
        cmd = " ".join(args)
        result = easy_run.fully_buffered(command=cmd)
        tst = [os.path.exists("sv_refined_18deg" + e) for e in [".json", ".pickle"]]
        if tst.count(True) != 2:
            return "18 deg interval width scan-varying refinement output was not found"

        # run analysis for static refinement job
        os.chdir(directory)
        if not self._centroid_analysis("refined_static.pickle"):
            return "Failed to run analysis on static refinement job"

        # run analysis for 54 deg scan-varying refinement job
        if not self._centroid_analysis("sv_refined_54deg.pickle", "_02"):
            return "Failed to run analysis on sv_refined_54deg.pickle"

        # run analysis for 36 deg scan-varying refinement job
        if not self._centroid_analysis("sv_refined_36deg.pickle", "_03"):
            return "Failed to run analysis on sv_refined_36deg.pickle"

        # run analysis for 18 deg scan-varying refinement job
        if not self._centroid_analysis("sv_refined_18deg.pickle", "_04"):
            return "Failed to run analysis on sv_refined_18deg.pickle"

        # Return empty status for success
        return ""

    def _centroid_analysis(self, filename, tag="_01"):
        existing_pngs = set(glob.glob("*.png"))
        cmd = "dials.python {0} {1}".format(
            self._ca, os.path.join(self._working_dir, filename)
        )
        result = easy_run.fully_buffered(command=cmd)
        if result.return_code != 0:
            return False
        # append tag to the plot filenames
        pngs = set(glob.glob("*.png"))
        new_pngs = pngs.difference(existing_pngs)
        for f in new_pngs:
            root, ext = os.path.splitext(f)
            os.rename(f, root + tag + ext)
        return True

    def find_experiments_json(self, idx_path):
        """Given the path to an indexed.pickle, try to identify an associated
    experiments.json in the same directory assuming xia2's naming convention"""

        head, tail = os.path.split(idx_path)
        root, ext = os.path.splitext(tail)
        exp_path = os.path.join(head, root.replace("indexed", "experiments") + ".json")
        if os.path.exists(exp_path):
            return exp_path
        else:
            return None


if __name__ == "__main__":
    from dials.util import halraiser

    try:
        script = Script()
        script.run()
    except Exception as e:
        halraiser(e)
