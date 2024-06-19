"""
Script to compare various methods of indexing ED stills

"""

import os
import glob
import subprocess
import shutil
from dials.array_family import flex
from libtbx.table_utils import simple_table
from dxtbx.model.experiment_list import ExperimentListFactory
from dials.algorithms.indexing.compare_orientation_matrices import (
    difference_rotation_matrix_axis_angle,
)
from time import perf_counter_ns
import json


class Script(object):
    def __init__(
        self,
        prefix="noiseimage_",
        extension="img",
        spotfind_cmds=[
            "min_spot_size=3",
        ],
    ):

        self.prefix = prefix
        self.extension = extension
        self.spotfind_cmds = spotfind_cmds

        self.cmds = [
            (
                shutil.which("dials.index"),
                "imported.expt",
                "strong.refl",
                "stills.indexer=sequences",
                "unit_cell=78.840,78.840,38.290,90.000,90.000,90.000",
                "space_group=P43212",
                "n_macro_cycles=2",
                "detector.fix=distance",
                "indexing.method=real_space_grid_search",
            ),  # RSGS
            (
                shutil.which("dials.index"),
                "imported.expt",
                "strong.refl",
                "stills.indexer=sequences",
                "unit_cell=78.840,78.840,38.290,90.000,90.000,90.000",
                "space_group=P43212",
                "n_macro_cycles=2",
                "detector.fix=distance",
                "indexing.method=low_res_spot_match",
            ),  # lrsm
            (
                shutil.which("dials.index"),
                "imported.expt",
                "strong.refl",
                "stills.indexer=sequences",
                "unit_cell=78.840,78.840,38.290,90.000,90.000,90.000",
                "space_group=P43212",
                "n_macro_cycles=2",
                "detector.fix=distance",
                "bootstrap_crystal=True",
                "indexing.method=low_res_spot_match",
            ),  # lrsm + bootstrap crystal
            (
                shutil.which("dials.index"),
                "imported.expt",
                "strong.refl",
                "stills.indexer=sequences",
                "unit_cell=78.840,78.840,38.290,90.000,90.000,90.000",
                "space_group=P43212",
                "n_macro_cycles=2",
                "detector.fix=distance",
                "indexing.method=pink_indexer",
                "percent_bandwidth=2",
                "min_lattices=5",
            ),  # pink_indexer
        ]

        self.images = sorted(glob.glob(self.prefix + "*." + self.extension))
        print(
            f"Found {len(self.images)} images matching {self.prefix}*.{self.extension} in the current working directory"
        )

    def run(self):
        results = [self.process(image) for image in self.images]
        with open("indigo_results.json", "w") as f:
            json.dump(results, f)

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
        serialno = image.lstrip(self.prefix).rstrip("." + self.extension)
        exp1 = ExperimentListFactory.from_json_file("experiments_" + serialno + ".json")
        exp2 = ExperimentListFactory.from_json_file("indexed.expt")
        R_ij, axis, angle, cb_op = difference_rotation_matrix_axis_angle(
            exp1[0].crystal, exp2[0].crystal
        )
        return abs(angle)

    def process(self, image):
        print(f"processing {image}")
        cmd = (shutil.which("dials.import"), image)
        result = subprocess.run(cmd)
        cmd = (shutil.which("dials.find_spots"), "imported.expt") + self.spotfind_cmds
        result = subprocess.run(cmd)
        strong = flex.reflection_table.from_file("strong.refl")
        d = {"nspots": len(strong)}

        nindexed = []
        offset_deg = []
        elapsed = []
        for cmd in self.cmds:
            start = perf_counter_ns()
            result = subprocess.run(cmd)
            elapsed.append(perf_counter_ns() - start)
            if os.path.isfile("indexed.refl"):
                idx = flex.reflection_table.from_file("indexed.refl")
                nindexed.append(idx.get_flags(idx.flags.indexed).count(True))
                offset_deg.append(self.compare_orientation_matrices(image))
                os.remove("indexed.refl")
                os.remove("indexed.expt")
            else:
                nindexed.append(0)
                offset_deg.append(None)

        d["nindexed"] = nindexed
        d["offset_deg"] = offset_deg
        d["elapsed"] = elapsed

        print(
            "#Spots:{0} ".format(d["nspots"])
            + " ".join([str(i + 1) + ":" + str(e) for i, e in enumerate(nindexed)])
        )

        return d


if __name__ == "__main__":
    # Edit prefix, extension and spotfinding commands as needed
    script = Script(
        prefix="noiseimage_",
        extension="img",
        spotfind_cmds=[
            "min_spot_size=3",
        ],
    )
    script.run()
