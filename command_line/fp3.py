import procrunner
import sys
import os
import glob
import copy
import concurrent.futures

from iotbx import phil
from dxtbx.model.experiment_list import ExperimentList, ExperimentListFactory
from dials_scratch.fp3 import even_blocks, index_blocks, format_phil_include
from dials_scratch.fp3 import nproc, combine_reflections, combine_experiments
from dials_scratch.fp3 import find_setup_script

scope = phil.parse(
    """
find_spots {
  include scope dials.command_line.find_spots.phil_scope
}
index {
  include scope dials.command_line.index.phil_scope
}
refine {
  include scope dials.command_line.refine.phil_scope
}
integrate {
  include scope dials.command_line.integrate.phil_scope
}
symmetry {
  include scope dials.command_line.symmetry.phil_scope
}
scale {
  include scope dials.command_line.scale.phil_scope
}
resolution {
  include scope dials.command_line.estimate_resolution.phil_scope
}
parallelism = *process drmaa
  .type = choice
max_workers = 1
  .type = int
  .help = "Maximum number of concurrent workers"
worker_nproc = 1
  .type = int
  .help = "Number of processes for each worker"
""",
    process_includes=True,
)


class FP3:
    def __init__(self, filenames, params):
        self._experiment = ExperimentListFactory.from_filenames(filenames)

        # parse PHIL parameters
        clai = scope.command_line_argument_interpreter()
        self._working = scope.fetch(clai.process_and_fetch(params))
        self._params = self._working.extract()

        print(scope.fetch_diff(self._working).as_str())

        self._crystal = None
        self._root = os.getcwd()
        self._n = nproc()

        # quick checks...
        scan = self._experiment[0].scan
        osc = scan.get_oscillation()
        rng = scan.get_image_range()

        self._osc = osc
        self._rng = rng

        self._integrated = []
        self._combined = None
        self._symmetry = None
        self._scaled = None

    def _write_phil(self, what, where):
        """helper function: extract (what) from global scope and write as
        (where)/(what).phil iff it contains anything. If file written,
        returns [filename] else [] to allow straight command-line adding."""

        phil = format_phil_include(scope, self._working, what).strip()
        if not phil:
            return []

        out = os.path.join(where, f"{what}.phil")
        with open(out, "w") as f:
            f.write(phil)
        return [out]

    def index(self):
        """Initial find spots and indexing: if it has been run already will
        just reload the results from the previous run."""

        work = os.path.join(self._root, "index")
        if os.path.exists(work):
            indexed = ExperimentList.from_file(os.path.join(work, "indexed.expt"))

            self._experiment[0].crystal = indexed[0].crystal
            return

        if not os.path.exists(work):
            os.mkdir(work)

        self._experiment.as_file(os.path.join(work, "input.expt"))

        five = int(round(5 / self._osc[1]))
        i0, i1 = self._rng
        blocks = [(b[0] + 1, b[1]) for b in index_blocks(i0 - 1, i1, self._osc[1])]

        phil = self._write_phil("find_spots", work)
        result = procrunner.run(
            ["dials.find_spots", "input.expt", "nproc=%d" % self._n]
            + phil
            + ["scan_range=%d,%d" % block for block in blocks],
            working_directory=work,
        )

        # let's just assume that was fine - so index

        phil = self._write_phil("index", work)
        result = procrunner.run(
            ["dials.index", "input.expt", "strong.refl"] + phil, working_directory=work
        )

        indexed = ExperimentList.from_file(os.path.join(work, "indexed.expt"))

        self._experiment[0].crystal = indexed[0].crystal

    def integrate(self):
        """Integration of the complete scan: will split the data into 5 deg
        chunks and spin the integration of each chunk off separately"""

        rng = self._rng

        nblocks = int(round(self._osc[1] * (rng[1] - rng[0] + 1) / 5.0))
        blocks = even_blocks(rng[0] - 1, rng[1], nblocks)

        with concurrent.futures.ProcessPoolExecutor(
            max_workers=self._params.max_workers
        ) as pool:
            jobs = []
            for j, block in enumerate(blocks):
                if self._params.parallelism == "process":
                    jobs.append(pool.submit(self.integrate_chunk, j, block))
                else:
                    jobs.append(pool.submit(self.integrate_chunk_script, j, block))
            for job in concurrent.futures.as_completed(jobs):
                self._integrated.append(job.result())

    def integrate_chunk(self, no, chunk):
        """Integrate a chunk of data: performs -
        - spot finding
        - indexing by using the UB matrix determined above (quick)
        - scan varying refinement
        - integration
        And works in the usual way for DIALS of using spots from every
        image for modelling. This is designed to be data-local e.g. could
        somehow stash the data as read for spot finding and not need to
        read more times in the integration."""

        # FIXME this should probably be logging
        print(f"Processing block {no} for images {chunk[0]} to {chunk[1]}")

        work = os.path.join(self._root, "integrate%03d" % no)
        if os.path.exists(work):
            if all(
                os.path.exists(os.path.join(work, f"integrated.{exten}"))
                for exten in ["refl", "expt"]
            ):
                return work

        if not os.path.exists(work):
            os.mkdir(work)

        # fix up the scan to correspond to input chunk
        expt = copy.deepcopy(self._experiment)
        expt[0].scan = expt[0].scan[chunk[0] : chunk[1]]

        expt.as_file(os.path.join(work, "input.expt"))

        phil = self._write_phil("find_spots", work)
        # FIXME nproc here depends on max_workers and parallelism mode
        np = self._params.worker_nproc

        result = procrunner.run(
            ["dials.find_spots", "input.expt", f"nproc={np}"] + phil,
            working_directory=work,
            print_stdout=False,
            print_stderr=False,
        )

        phil = self._write_phil("index", work)
        result = procrunner.run(
            [
                "dials.index",
                "input.expt",
                "strong.refl",
                "index_assignment.method=local",
            ]
            + phil,
            working_directory=work,
            print_stdout=False,
            print_stderr=False,
        )

        phil = self._write_phil("refine", work)
        result = procrunner.run(
            ["dials.refine", "indexed.expt", "indexed.refl"] + phil,
            working_directory=work,
            print_stdout=False,
            print_stderr=False,
        )

        phil = self._write_phil("integrate", work)
        # FIXME nproc here depends on max_workers and parallelism mode
        result = procrunner.run(
            ["dials.integrate", "refined.expt", "refined.refl", f"nproc={np}"] + phil,
            working_directory=work,
            print_stdout=False,
            print_stderr=False,
        )

        return work

    def integrate_chunk_script(self, no, chunk):
        """Integrate a chunk of data: performs -
        - spot finding
        - indexing by using the UB matrix determined above (quick)
        - scan varying refinement
        - integration
        And works in the usual way for DIALS of using spots from every
        image for modelling. This is designed to be data-local e.g. could
        somehow stash the data as read for spot finding and not need to
        read more times in the integration. This writes one script for
        submission to a cluster to do the processing."""

        # FIXME this should probably be logging
        print(f"Writing script {no} for images {chunk[0]} to {chunk[1]}")

        # first check if there is nothing to be done here
        work = os.path.join(self._root, "integrate%03d" % no)
        if os.path.exists(work):
            if all(
                os.path.exists(os.path.join(work, f"integrated.{exten}"))
                for exten in ["refl", "expt"]
            ):
                return work

        if not os.path.exists(work):
            os.mkdir(work)

        # fix up the scan to correspond to input chunk, save to working area
        expt = copy.deepcopy(self._experiment)
        expt[0].scan = expt[0].scan[chunk[0] : chunk[1]]

        expt.as_file(os.path.join(work, "input.expt"))

        # this is used to (i) write script and (ii) submit to cluster
        np = self._params.worker_nproc

        fout = open(os.path.join(work, "integrate.sh"), "w")

        fout.write("\n".join(["#!/bin/bash", f". {find_setup_script()}", "", ""]))

        phil = self._write_phil("find_spots", work)
        fout.write(f"dials.find_spots input.expt nproc={np}" + " ".join(phil) + "\n")
        fout.write("if [ $? -ne 0 ] ; then exit $? ; fi\n")

        phil = self._write_phil("index", work)
        fout.write(
            " ".join(
                [
                    "dials.index",
                    "input.expt",
                    "strong.refl",
                    "index_assignment.method=local",
                ]
            )
            + " ".join(phil)
            + "\n"
        )
        fout.write("if [ $? -ne 0 ] ; then exit $? ; fi\n")

        phil = self._write_phil("refine", work)
        fout.write(
            " ".join(["dials.refine", "indexed.expt", "indexed.refl"])
            + " ".join(phil)
            + "\n"
        )
        fout.write("if [ $? -ne 0 ] ; then exit $? ; fi\n")

        phil = self._write_phil("integrate", work)
        fout.write(
            " ".join(["dials.integrate", "refined.expt", "refined.refl", f"nproc={np}"])
            + " ".join(phil)
            + "\n"
        )
        fout.write("if [ $? -ne 0 ] ; then exit $? ; fi\n")

        fout.close()

        # FIXME submit the script for processing
        print(f"Executing script {no} for images {chunk[0]} to {chunk[1]}")
        result = procrunner.run(
            ["bash", "integrate.sh"],
            working_directory=work,
            print_stdout=False,
            print_stderr=False,
        )

        return work

    def combine(self):
        """Collect together the data so far integrated."""

        work = os.path.join(self._root, "combine")

        if os.path.exists(work):
            if all(
                os.path.exists(os.path.join(work, f"combined.{exten}"))
                for exten in ["refl", "expt"]
            ):
                self._combined = work
                return

        if not os.path.exists(work):
            os.mkdir(work)

        assert self._integrated is not []

        integrated_refl = [
            os.path.join(directory, "integrated.refl")
            for directory in sorted(self._integrated)
        ]
        integrated_expt = [
            os.path.join(directory, "integrated.expt")
            for directory in sorted(self._integrated)
        ]

        combine_reflections(integrated_refl, os.path.join(work, "combined.refl"))
        combine_experiments(integrated_expt, os.path.join(work, "combined.expt"))

        self._combined = work

    def symmetry(self):
        """Take the combined data, determine the symmetry."""

        assert self._combined is not None

        work = os.path.join(self._root, "symmetry")

        if os.path.exists(work):
            if all(
                os.path.exists(os.path.join(work, f"symmetrized.{exten}"))
                for exten in ["refl", "expt"]
            ):
                self._symmetry = work
                return work

        if not os.path.exists(work):
            os.mkdir(work)

        phil = self._write_phil("symmetry", work)
        result = procrunner.run(
            ["dials.symmetry"]
            + phil
            + [
                os.path.join(self._combined, f"combined.{exten}")
                for exten in ["refl", "expt"]
            ],
            working_directory=work,
        )

        self._symmetry = work

    def scale(self, d_min=None):
        """Scale the data, to the resolution input or everything if unspecified.
        If the data have been previously scaled, start from the scaled data
        for sake of efficiency."""

        work = os.path.join(self._root, "scale")

        if os.path.exists(work):
            if d_min is None and all(
                os.path.exists(os.path.join(work, f"scaled.{exten}"))
                for exten in ["refl", "expt"]
            ):
                self._scaled = work
                return work

        if not os.path.exists(work):
            os.mkdir(work)

        if self._scaled is not None:
            source = os.path.join(self._scaled, "scaled")
        else:
            source = os.path.join(self._symmetry, "symmetrized")

        phil = self._write_phil("scale", work)
        command = ["dials.scale"] + phil
        if d_min is not None:
            command.append(f"d_min={d_min}")

        result = procrunner.run(
            command + [f"{source}.{exten}" for exten in ["refl", "expt"]],
            working_directory=work,
        )

        self._scaled = work

    def resolution(self):
        """Determine the resolution of the data after scaling."""

        assert self._scaled

        source = os.path.join(self._scaled, "scaled")

        work = os.path.join(self._root, "scale")
        phil = self._write_phil("resolution", work)

        result = procrunner.run(
            ["dials.estimate_resolution"]
            + phil
            + [f"{source}.{exten}" for exten in ["refl", "expt"]],
            working_directory=self._scaled,
        )

        d_min = None
        for record in result["stdout"].split(b"\n"):
            if record.startswith(b"Resolution cc_half:"):
                d_min = float(record.split()[-1])

        return d_min


if __name__ == "__main__":
    args = [arg for arg in sys.argv[1:] if not "=" in arg]
    params = [arg for arg in sys.argv[1:] if "=" in arg]
    filenames = sum(map(glob.glob, args), [])

    fp3 = FP3(filenames, params)
    fp3.index()
    fp3.integrate()
    fp3.combine()
    fp3.symmetry()
    fp3.scale()
    d_min = fp3.resolution()
    fp3.scale(d_min=d_min)
