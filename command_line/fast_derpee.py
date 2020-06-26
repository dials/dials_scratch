import procrunner
import sys
import os
import glob
import copy

from dxtbx.model.experiment_list import ExperimentList, ExperimentListFactory


class derpee:
    def __init__(self, filenames):
        self._experiment = ExperimentListFactory.from_filenames(filenames)
        self._crystal = None
        self._root = os.getcwd()

        # quick checks...
        scan = self._experiment[0].scan
        osc = scan.get_oscillation()
        rng = scan.get_image_range()

        wedge = osc[1] * (rng[1] - rng[0] + 1)

        assert wedge >= 95

        self._osc = osc
        self._rng = rng

    def index(self):
        """Initial find spots and indexing: if it has been run already will
        just reload the results from the previous run."""

        work = os.path.join(self._root, "index")
        if os.path.exists(work):
            indexed = ExperimentList.from_file(os.path.join(work, "indexed.expt"))

            self._experiment[0].crystal = indexed[0].crystal
            return

        os.mkdir(work)

        # index from 0-5, 45-50 degree blocks - hard code on 0.1 degree frames
        # to start, do properly laters
        self._experiment.as_file(os.path.join(work, "input.expt"))

        blocks = [
            (start + 1, start + int(round(5 / self._osc[1])))
            for start in (
                0,
                int(round(45 / self._osc[1])),
                int(round(90 / self._osc[1])),
            )
        ]

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

        size = int(round(5 / self._osc[1]))

        # need to figure out how to spin this off to somehing running on a
        # cluster node... ideally want this called on many chunks at once

        # count chunks in computer numbers
        for j, start in enumerate(range(self._rng[0] - 1, self._rng[1], size)):
            self.integrate_chunk(j, (start, start + size))

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

        work = os.path.join(self._root, "integrate%02d" % no)
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

    def symmetry_scale(self):
        """Collect together the data so far integrated, use e.g. multiplex to
        combine, determine the symmetry and scale, or combine experiments
        followed by dials.symmetry and dials.scale #TODO. Since the UB matrix
        is in principle the same for each scan should be fine."""
        pass


if __name__ == "__main__":
    filenames = sum(map(glob.glob, sys.argv[1:]), [])

    derp = derpee(filenames)
    derp.index()
    derp.integrate()
