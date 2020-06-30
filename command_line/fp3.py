import procrunner
import sys
import os
import glob
import copy

from dxtbx.model.experiment_list import ExperimentList, ExperimentListFactory
from dials_scratch.fp3 import even_blocks, index_blocks

class fp3:
    def __init__(self, filenames):
        self._experiment = ExperimentListFactory.from_filenames(filenames)
        self._crystal = None
        self._root = os.getcwd()

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
        blocks = [(b[0]+1, b[1]) for b in index_blocks(i0-1, i1, self._osc[1])]

        result = procrunner.run(
            ["dials.find_spots", "input.expt", "nproc=8"]
            + ["scan_range=%d,%d" % block for block in blocks],
            working_directory=work,
        )

        # let's just assume that was fine - so index

        result = procrunner.run(
            ["dials.index", "input.expt", "strong.refl"], working_directory=work
        )

        indexed = ExperimentList.from_file(os.path.join(work, "indexed.expt"))

        self._experiment[0].crystal = indexed[0].crystal

    def integrate(self):
        """Integration of the complete scan: will split the data into 5 deg
        chunks and spin the integration of each chunk off separately"""

        rng = self._rng
        
        nblocks = int(round(self._osc[1] * (rng[1] - rng[0] + 1) / 5.0))
        blocks = even_blocks(rng[0] - 1, rng[1], nblocks)
        
        # need to figure out how to spin this off to somehing running on a
        # cluster node... ideally want this called on many blocks at once

        for j, block in enumerate(blocks):
            self._integrated.append(self.integrate_chunk(j, block))

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

        work = os.path.join(self._root, "integrate%03d" % no)
        if os.path.exists(work):
            if all(
                os.path.exists(os.path.join(work, f"integrated.{exten}"))
                for exten in ["refl", "expt"]
            ):
                return work

        if not os.path.exists(work):
            os.mkdir(work)

        expt = copy.deepcopy(self._experiment)

        # fix up the scan to correspond to input chunk
        scan = expt[0].scan
        epochs = scan.get_epochs()[chunk[0] : chunk[1]]
        exp_times = scan.get_exposure_times()[chunk[0] : chunk[1]]
        scan.set_image_range((chunk[0] + 1, chunk[1]))
        scan.set_oscillation(
            (self._osc[0] + (chunk[0] - self._rng[0]) * self._osc[1], self._osc[1])
        )
        scan.set_epochs(epochs)
        scan.set_exposure_times(exp_times)
        expt[0].scan = scan

        expt.as_file(os.path.join(work, "input.expt"))
        result = procrunner.run(
            ["dials.find_spots", "input.expt", "nproc=8"], working_directory=work
        )
        result = procrunner.run(
            [
                "dials.index",
                "input.expt",
                "strong.refl",
                "index_assignment.method=local",
            ],
            working_directory=work,
        )
        result = procrunner.run(
            ["dials.refine", "indexed.expt", "indexed.refl"], working_directory=work
        )

        result = procrunner.run(
            ["dials.integrate", "refined.expt", "refined.refl", "nproc=8"],
            working_directory=work,
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

        integrated = sum(
            [
                [
                    os.path.join(directory, f"integrated.{exten}")
                    for exten in ["refl", "expt"]
                ]
                for directory in sorted(self._integrated)
            ],
            [],
        )

        result = procrunner.run(
            ["dials.combine_experiments"] + integrated, working_directory=work
        )

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

        result = procrunner.run(
            ["dials.symmetry"]
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

        command = ["dials.scale"]
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

        result = procrunner.run(
            ["dials.resolutionizer"]
            + [f"{source}.{exten}" for exten in ["refl", "expt"]],
            working_directory=self._scaled,
        )

        d_min = None
        for record in result["stdout"].split(b"\n"):
            if record.startswith(b"Resolution cc_half:"):
                d_min = float(record.split()[-1])

        return d_min


if __name__ == "__main__":
    filenames = sum(map(glob.glob, sys.argv[1:]), [])

    derp = fp3(filenames)
    derp.index()
    derp.integrate()
    derp.combine()
    derp.symmetry()
    derp.scale()
    d_min = derp.resolution()
    derp.scale(d_min=d_min)
