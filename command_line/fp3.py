# LIBTBX_SET_DISPATCHER_NAME dials.fp3

import procrunner
import sys
import os
import glob
import copy
import concurrent.futures
import logging
import time
import functools

logger = logging.getLogger("fp3")

from iotbx import phil
from dxtbx.model.experiment_list import ExperimentList, ExperimentListFactory
from dials_scratch.fp3 import even_blocks, index_blocks, format_phil_include
from dials_scratch.fp3 import nproc, combine_reflections, combine_experiments
from dials_scratch.fp3 import find_setup_script, prime_factors

# FIXME add block size for processing to the input parameters

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
log_level = CRITICAL ERROR WARNING *INFO DEBUG
  .type = choice
""",
    process_includes=True,
)


class FP3:
    def __init__(self, filenames, params):
        self._experiments = ExperimentListFactory.from_filenames(filenames)

        # parse PHIL parameters
        clai = scope.command_line_argument_interpreter()
        self._working = scope.fetch(clai.process_and_fetch(params))
        self._params = self._working.extract()
        self._debug = self._params.log_level == "DEBUG"

        # configure logging
        logging.basicConfig(
            filename="dials.fp3.log", filemode="w", format="%(message)s"
        )
        logger.addHandler(logging.StreamHandler(sys.stdout))
        logger.setLevel(getattr(logging, self._params.log_level))
        logger.info(scope.fetch_diff(self._working).as_str())

        logger.info(f"Found {len(self._experiments)} scans to process")

        if len(self._experiments) > 1:
            sys.exit("Multi-sweep data not supported at this time")

        self._crystal = None
        self._root = os.getcwd()
        self._n = nproc()

        # quick checks...
        # FIXME MULTI-SWEEEP will require attention
        scan = self._experiments[0].scan
        osc = scan.get_oscillation()

        self._osc = osc
        self._image_range = scan.get_image_range()

        # default nproc here depends on max_workers and parallelism mode -
        # if both are set to 1 have a guess for something reasonable - factor
        # nproc to primes then share these out between workers and cores...

        np = self._params.worker_nproc
        nw = self._params.max_workers
        if np == 1 and nw == 1:
            f = prime_factors(self._n)
            nw = functools.reduce((lambda x, y: x * y), f[0::2])
            np = functools.reduce((lambda x, y: x * y), f[1::2])
            self._params.worker_nproc = np
            self._params.max_workers = nw

        logger.info(f"Using {self._n} processors to for images {self._image_range}")
        logger.info(
            f"with up to {self._params.max_workers} x {self._params.worker_nproc} core workers"
        )

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
        if os.path.exists(os.path.join(work, "indexed.expt")):
            indexed = ExperimentList.from_file(os.path.join(work, "indexed.expt"))

            self._experiments[0].crystal = indexed[0].crystal
            logger.info("Picked up pre-existing crystal:")
            logger.info(indexed[0].crystal)
            return

        if not os.path.exists(work):
            os.mkdir(work)

        self._experiments.as_file(os.path.join(work, "input.expt"))

        # FIXME MULTI-SWEEEP will require attention - for 5° blocks for each
        # sweep will need to ensure that we compute the right blocks for each
        five = int(round(5 / self._osc[1]))
        i0, i1 = self._image_range
        blocks = [(b[0] + 1, b[1]) for b in index_blocks(i0 - 1, i1, self._osc[1])]

        logger.info("Finding spots...")
        phil = self._write_phil("find_spots", work)
        result = procrunner.run(
            ["dials.find_spots", "input.expt", "nproc=%d" % self._n]
            + phil
            + ["scan_range=%d,%d" % block for block in blocks],
            working_directory=work,
            print_stdout=self._debug,
            print_stderr=self._debug,
        )

        logger.info("Indexing...")
        # let's just assume that was fine - so index

        phil = self._write_phil("index", work)
        result = procrunner.run(
            ["dials.index", "input.expt", "strong.refl"] + phil,
            working_directory=work,
            print_stdout=self._debug,
            print_stderr=self._debug,
        )

        # FIXME MULTI-LATTICE quite possible at this stage that we have more
        # than one lattice, on more than one experiment, so probably want to
        # collate the unique crystals at this stage

        indexed = ExperimentList.from_file(os.path.join(work, "indexed.expt"))
        logger.info(indexed[0].crystal)
        self._experiments[0].crystal = indexed[0].crystal

    def integrate(self):
        """Integration of the complete scan: will split the data into 5 deg
        chunks and spin the integration of each chunk off separately"""

        # FIXME MULTI-SWEEP - in here if we have multiple sweeps (i.e. more than
        # one imageset) need to integrate these separately, with the appropriate
        # block sizes set. Probably want to form a list of tasks which includes
        # all blocks for all image sets for parallel processing rather than
        # iterating through the imagesets doing parallel processing on each.

        irange = self._image_range

        nblocks = int(round(self._osc[1] * (irange[1] - irange[0] + 1) / 5.0))
        blocks = even_blocks(irange[0] - 1, irange[1], nblocks)

        logger.info("Integrating...")

        # FIXME in here define the integrate directory template based on the
        # maximum number of blocks => don't hard code this. Should also have
        # the imageset in the name e.g. integrate 0.01 for the second block of
        # imageset 0.

        # API wise, if we make blocks a list of tuples here this could be pretty
        # simple to modify i.e. (imageset, block)

        if self._params.parallelism == "process":
            with concurrent.futures.ProcessPoolExecutor(
                max_workers=self._params.max_workers
            ) as pool:
                jobs = []
                for j, block in enumerate(blocks):
                    jobs.append(pool.submit(self.integrate_chunk, j, block))
                for job in concurrent.futures.as_completed(jobs):
                    self._integrated.append(job.result())
        else:
            self.integrate_drmaa_array(blocks)

        # if debug, print the log files from each of the integration stages
        if self._debug:
            for j, block in enumerate(blocks):
                work = os.path.join(self._root, "integrate%04d" % no)
                for log in (
                    "dials.find_spots.log",
                    "dials.index.log",
                    "dials.refine.log",
                    "dials.integrate.log",
                ):
                    for record in open(os.path.join(work, log)):
                        logger.debug(record[:-1])

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

        logger.debug(f"Processing block {no} for images {chunk[0]} to {chunk[1]}")

        # FIXME MULTI-SWEEEP in here dig out the right imageset etc.

        work = os.path.join(self._root, "integrate%04d" % no)
        if os.path.exists(work):
            if all(
                os.path.exists(os.path.join(work, f"integrated.{exten}"))
                for exten in ["refl", "expt"]
            ):
                return work

        if not os.path.exists(work):
            os.mkdir(work)

        # FIXME MULTI-SWEEEP pick the right scan etc. - also need to pick the
        # _right_ experiments to copy i.e. those for which the image set is the
        # _right_ image set.

        # fix up the scan to correspond to input chunk
        expt = copy.copy(self._experiments)
        expt[0].scan = expt[0].scan[chunk[0] : chunk[1]]

        expt.as_file(os.path.join(work, "input.expt"))

        phil = self._write_phil("find_spots", work)

        np = self._params.worker_nproc

        result = procrunner.run(
            ["dials.find_spots", "input.expt", f"nproc={np}"] + phil,
            working_directory=work,
            print_stdout=False,
            print_stderr=False,
        )

        # FIXME check the return status, check for errors, verify that the
        # files we are expecting to exist actually ... exist

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

        # FIXME check the return status, check for errors, verify that the
        # files we are expecting to exist actually ... exist

        phil = self._write_phil("refine", work)
        result = procrunner.run(
            ["dials.refine", "indexed.expt", "indexed.refl"] + phil,
            working_directory=work,
            print_stdout=False,
            print_stderr=False,
        )

        # FIXME check the return status, check for errors, verify that the
        # files we are expecting to exist actually ... exist

        phil = self._write_phil("integrate", work)

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

        logger.debug(f"Writing script {no} for images {chunk[0]} to {chunk[1]}")

        # first check if there is nothing to be done here
        work = os.path.join(self._root, f"integrate{no:04}")
        if os.path.exists(work):
            if all(
                os.path.exists(os.path.join(work, f"integrated.{exten}"))
                for exten in ["refl", "expt"]
            ):
                return work

        if not os.path.exists(work):
            os.mkdir(work)

        # fix up the scan to correspond to input chunk, save to working area
        expt = copy.deepcopy(self._experiments)
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

        return work

    def integrate_drmaa_array(self, blocks):
        script_file = os.path.join(self._root, "run_integrate_array.sh")
        with open(script_file, "w") as script:

            script.write("#!/bin/bash\n")
            nblocks = 0
            for idx, chunk in enumerate(blocks, start=1):

                working_dir = self.integrate_chunk_script(idx, chunk)
                script.write(f"WORKING_DIR_{idx}={working_dir}\n")
                self._integrated.append(working_dir)
                nblocks += 1

            script.write("TASK_WORKING_DIR=WORKING_DIR_${SGE_TASK_ID}\n")
            script.write("cd ${!TASK_WORKING_DIR}\n")
            script.write(
                "sh ./integrate.sh > ${!TASK_WORKING_DIR}/integrate.out  2> ${!TASK_WORKING_DIR}/integrate.err"
            )

        import drmaa

        # TODO in here check if processing already performed before submission
        # of task...

        with drmaa.Session() as session:
            job = session.createJobTemplate()
            job.jobName = "fp3_integrate"
            job.workingDirectory = self._root
            job.remoteCommand = "sh"
            args = [script_file]
            job.args = args
            job.jobCategory = "medium"
            if self._params.worker_nproc > 1:
                smp = f"-pe smp {self._params.worker_nproc}"
            else:
                smp = ""
            job.nativeSpecification = f"-V {smp} -l h_rt=12:00:00 -l mfree=4G -tc {self._params.max_workers}  -o /dev/null -e /dev/null"

            job_ids = session.runBulkJobs(job, 1, nblocks, 1)
            session.synchronize(job_ids, drmaa.Session.TIMEOUT_WAIT_FOREVER, True)
            session.deleteJobTemplate(job)

        # FIXME check the return status, check for errors, verify that the
        # files we are expecting to exist actually ... exist, and if not, do
        # something sensible.

    def combine(self):
        """Collect together the data so far integrated."""

        logger.info("Combining...")
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

        # FIXME this should not be an assertion - if [] then sys.exit()
        # with helpful message
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

        logger.info("Determining symmetry...")
        work = os.path.join(self._root, "symmetry")

        if os.path.exists(work):
            if all(
                os.path.exists(os.path.join(work, f"symmetrized.{exten}"))
                for exten in ["refl", "expt"]
            ):
                symm = ExperimentList.from_file(os.path.join(work, "symmetrized.expt"))
                logger.info(symm[0].crystal)
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
            print_stdout=self._debug,
            print_stderr=self._debug,
        )

        symm = ExperimentList.from_file(os.path.join(work, "symmetrized.expt"))
        logger.info(symm[0].crystal)
        self._symmetry = work

    def scale(self, d_min=None):
        """Scale the data, to the resolution input or everything if unspecified.
        If the data have been previously scaled, start from the scaled data
        for sake of efficiency."""

        logger.info("Scaling...")
        work = os.path.join(self._root, "scale")

        if os.path.exists(work):
            if d_min is None and all(
                os.path.exists(os.path.join(work, f"scaled.{exten}"))
                for exten in ["refl", "expt"]
            ):
                self._scaled = work
                stdout = (
                    open(os.path.join(work, "dials.scale.log"), "rb")
                    .read()
                    .split(b"\n")
                )

                # find summary table
                for j, line in enumerate(stdout):
                    if b"--Summary of merging statistics--" in line:
                        break

                # save table of values
                self._stats = []
                for line in stdout[j + 2 :]:
                    if not line.strip():
                        break
                    self._stats.append(line.decode())

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
            print_stdout=self._debug,
            print_stderr=self._debug,
        )

        stdout = result.stdout.split(b"\n")

        # find summary table
        for j, line in enumerate(stdout):
            if b"----------Summary of merging statistics-----------" in line:
                break

        # save table of values
        self._stats = []
        for line in stdout[j + 2 :]:
            if not line.strip():
                break
            self._stats.append(line.decode())

        self._scaled = work

    def resolution(self):
        """Determine the resolution of the data after scaling."""

        assert self._scaled

        # FIXME check to see if we have done this before and if we have,
        # pick up the outcome from that calculation rather than re-running

        logger.info("Estimating resolution...")
        source = os.path.join(self._scaled, "scaled")

        work = os.path.join(self._root, "scale")
        phil = self._write_phil("resolution", work)

        result = procrunner.run(
            ["dials.estimate_resolution"]
            + phil
            + [f"{source}.{exten}" for exten in ["refl", "expt"]],
            working_directory=self._scaled,
            print_stdout=self._debug,
            print_stderr=self._debug,
        )

        d_min = None
        for record in result.stdout.split(b"\n"):
            if record.startswith(b"Resolution cc_half:"):
                d_min = float(record.split()[-1])

        if d_min is not None:
            logger.info(f"Resolution determined as {d_min:.2f}")
        else:
            logger.info("Recommendation is use all data")
        return d_min


if __name__ == "__main__":
    args = [arg for arg in sys.argv[1:] if not "=" in arg]
    params = [arg for arg in sys.argv[1:] if "=" in arg]
    filenames = sum(map(glob.glob, args), [])

    # FIXME move all this to a main() function

    fp3 = FP3(filenames, params)
    t0 = time.time()
    fp3.index()
    t1 = time.time()
    fp3.integrate()
    t2 = time.time()
    fp3.combine()
    t3 = time.time()
    fp3.symmetry()
    t4 = time.time()
    fp3.scale()
    t5 = time.time()
    d_min = fp3.resolution()
    if d_min:
        fp3.scale(d_min=d_min)
    t6 = time.time()
    logger.info("\n".join(fp3._stats))
    logger.info(f"Total processing time: {t6 - t0:0.2f}s")
