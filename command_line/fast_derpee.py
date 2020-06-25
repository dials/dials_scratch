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

        self._osc = osc[1]
        self._rng = rng

    def index(self):
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
            (start + 1, start + int(round(5 / self._osc)))
            for start in (0, int(round(45 / self._osc)), int(round(90 / self._osc)))
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
        size = int(round(5 / self._osc))
        # count chunks in computer numbers
        for j, start in enumerate(range(self._rng[0] - 1, self._rng[1], size)):
            self.integrate_chunk(j, (start, start + size))

    def integrate_chunk(self, no, chunk):
        # need to figure out how to spin this off to somehing running on a
        # cluster node...
        work = os.path.join(self._root, "integrate%02d" % no)
        if not os.path.exists(work):
            os.mkdir(work)

        expt = copy.deepcopy(self._experiment)

        # fix up the scan to correspond to input chunk
        scan = expt[0].scan
        epochs = scan.get_epochs()[chunk[0] : chunk[1]]
        exp_times = scan.get_exposure_times()[chunk[0] : chunk[1]]
        scan.set_image_range((chunk[0] + 1, chunk[1]))
        scan.set_epochs(epochs)
        scan.set_exposure_times(exp_times)
        expt[0].scan = scan

        expt.as_file(os.path.join(work, "input.expt"))
        result = procrunner.run(
            ["dials.find_spots", "input.expt", "nproc=8"], working_directory=work
        )
        result = procrunner.run(
            ["dials.index", "input.expt", "strong.refl"], working_directory=work
        )
        result = procrunner.run(
            ["dials.refine", "indexed.expt", "indexed.refl"], working_directory=work
        )

        result = procrunner.run(
            ["dials.integrate", "refined.expt", "refined.refl", "nproc=8"],
            working_directory=work,
        )


if __name__ == "__main__":
    filenames = sum(map(glob.glob, sys.argv[1:]), [])

    derp = derpee(filenames)
    derp.index()
    derp.integrate()
