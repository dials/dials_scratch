"""
Script to compare various methods of indexing ED stills

"""

import os
import glob
from libtbx import easy_run
from dials.array_family import flex
from libtbx.table_utils import simple_table
from dxtbx.model.experiment_list import ExperimentListFactory
from dials.algorithms.indexing.compare_orientation_matrices import (
    difference_rotation_matrix_axis_angle,
)


class Script(object):
    def __init__(self):

        self.cmds = [
            (
                "dials.index stills.indexer=sweeps imported.expt "
                "strong.refl unit_cell=78.840,78.840,38.290,90.000,90.000,90.000 "
                "space_group=P43212 indexing.method=real_space_grid_search "
                "n_macro_cycles=2 detector.fix_list=Dist"
            ),  # RSGS
            (
                "dev.dials.indigo imported.expt strong.refl "
                "unit_cell=78.840,78.840,38.290,90.000,90.000,90.000 space_group=P43212 "
                "stills.indexer=sweeps n_macro_cycles=2 detector.fix_list=Dist "
            ),  # indigo
            (
                "dev.dials.indigo imported.expt strong.refl "
                "unit_cell=78.840,78.840,38.290,90.000,90.000,90.000 space_group=P43212 "
                "stills.indexer=sweeps n_macro_cycles=2 detector.fix_list=Dist "
                "bootstrap_crystal=True "
            ),  # indigo + bootstrap crystal
        ]

        self.images = sorted(glob.glob("noiseimage_*.img"))
        print(
            "Found {0} images matching noiseimage_*.img in the current working directory".format(
                len(self.images)
            )
        )

    def run(self):
        results = [self.process(image) for image in self.images]

        header = (
            ["Image", "Num spots"]
            + [str(i + 1) for i, _ in enumerate(self.cmds)]
            + ["Best"]
        )
        rows = []
        for im, res in zip(self.images, results):
            offset = []
            for v in res["offset_deg"]:
                if v is None:
                    offset.append("fail")
                else:
                    offset.append("{:.3f}".format(v))

            row = [str(e1) + ":" + e2 for e1, e2 in zip(res["nindexed"], offset)]
            row = [im, str(res["nspots"])] + row

            nindexed = res["nindexed"]
            top = max(nindexed)
            best = " "
            for i, _ in enumerate(self.cmds):
                if nindexed[i] == top:
                    best += str(i + 1) + " "
            if res["nspots"] == top:
                best += "*"
            rows.append(row + [best])
        st = simple_table(rows, header)
        print(st.format())

    @staticmethod
    def compare_orientation_matrices(image):
        serialno = image.lstrip("noiseimage_").rstrip(".img")
        exp1 = ExperimentListFactory.from_json_file("experiments_" + serialno + ".json")
        exp2 = ExperimentListFactory.from_json_file("indexed.expt")
        R_ij, axis, angle, cb_op = difference_rotation_matrix_axis_angle(
            exp1[0].crystal, exp2[0].crystal
        )
        return abs(angle)

    def process(self, image):
        print("processing {0}".format(image))
        cmd = "dials.import {0}".format(image)
        result = easy_run.fully_buffered(command=cmd)
        cmd = "dials.find_spots imported.expt min_spot_size=3"
        result = easy_run.fully_buffered(command=cmd)
        strong = flex.reflection_table.from_pickle("strong.refl")
        d = {"nspots": len(strong)}

        nindexed = []
        offset_deg = []
        for cmd in self.cmds:
            success = True
            try:
                result = easy_run.fully_buffered(command=cmd).raise_if_errors()
                assert result.return_code == 0
            except RuntimeError:
                success = False

            if success:
                idx = flex.reflection_table.from_file("indexed.refl")
                nindexed.append(idx.get_flags(idx.flags.indexed).count(True))
                offset_deg.append(self.compare_orientation_matrices(image))
            else:
                nindexed.append(0)
                offset_deg.append(None)

            try:
                os.remove("indexed.refl")
                os.remove("indexed.expt")
            except OSError:
                pass

        d["nindexed"] = nindexed
        d["offset_deg"] = offset_deg

        print(
            "#Spots:{0} ".format(d["nspots"])
            + " ".join([str(i + 1) + ":" + str(e) for i, e in enumerate(nindexed)])
        )

        return d


if __name__ == "__main__":
    script = Script()
    script.run()
