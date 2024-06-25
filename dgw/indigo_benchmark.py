"""
Script to compare various methods of indexing ED stills

"""

import os
import sys
import glob
import subprocess
import shutil
from dials.array_family import flex
from libtbx.table_utils import simple_table
from dxtbx.model.experiment_list import ExperimentListFactory
from dials.algorithms.indexing.compare_orientation_matrices import (
    difference_rotation_matrix_axis_angle,
)
from time import perf_counter
import json


class Script(object):
    def __init__(
        self,
        prefix="noiseimage_",
        extension="img",
        unit_cell="78.840,78.840,38.290,90.000,90.000,90.000",
        space_group="P43212",
        import_cmds=[],
        spotfind_cmds=[
            "min_spot_size=3",
        ],
    ):

        self.prefix = prefix
        self.extension = extension
        self.import_cmds = import_cmds
        self.spotfind_cmds = spotfind_cmds

        self.cmds = [
            (
                shutil.which("dials.index"),
                "imported.expt",
                "strong.refl",
                "stills.indexer=sequences",
                f"unit_cell={unit_cell}",
                f"space_group={space_group}",
                "n_macro_cycles=2",
                "detector.fix=distance",
                "indexing.method=real_space_grid_search",
            ),  # RSGS
            (
                shutil.which("dials.index"),
                "imported.expt",
                "strong.refl",
                "stills.indexer=sequences",
                f"unit_cell={unit_cell}",
                f"space_group={space_group}",
                "n_macro_cycles=2",
                "detector.fix=distance",
                "indexing.method=low_res_spot_match",
            ),  # lrsm
            (
                shutil.which("dials.index"),
                "imported.expt",
                "strong.refl",
                "stills.indexer=sequences",
                f"unit_cell={unit_cell}",
                f"space_group={space_group}",
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
                f"unit_cell={unit_cell}",
                f"space_group={space_group}",
                "n_macro_cycles=2",
                "detector.fix=distance",
                "indexing.method=pink_indexer",
                "percent_bandwidth=2",
                "min_lattices=5",
            ),  # pink_indexer
        ]
        if os.path.isfile("restraint.phil"):
            self.cmds = [cmd + ("restraint.phil",) for cmd in self.cmds]

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

    def reindex_to_ground_truth(self, image):
        """Calculate reindexing operation required to take the indexed image to
        the ground truth model, reindex the experiments and reflections,
        and return the offset angle in degrees."""

        # Calculate the cb_op and rotation required
        serialno = image.replace(self.prefix, "").replace("." + self.extension, "")
        exp1 = ExperimentListFactory.from_json_file("experiments_" + serialno + ".json")
        exp2 = ExperimentListFactory.from_json_file("indexed.expt")
        R_ij, axis, angle, cb_op = difference_rotation_matrix_axis_angle(
            exp1[0].crystal, exp2[0].crystal
        )

        # Do the reindexing
        # exp2 = reindex_experiments(exp2, cb_op)
        # exp2.as_file("indexed.expt")
        # indexed = flex.reflection_table.from_file("indexed.refl")
        # indexed = reindex_reflections(indexed, cb_op)
        # indexed.as_file("indexed.refl")
        cmd = (
            shutil.which("dials.reindex"),
            "indexed.expt",
            "indexed.refl",
            f"change_of_basis_op={str(cb_op)}",
            "output.experiments=indexed.expt",
            "output.reflections=indexed.refl",
        )
        result = subprocess.run(cmd)
        return abs(angle)

    def check_hkl(self, idx):
        ground_truth = flex.reflection_table.from_file("ground_truth.refl")

        # Match reflections by position
        id1, id2, _ = ground_truth.match(idx)
        ground_truth = ground_truth.select(id1)
        idx = idx.select(id2)

        # Select only reflections indexed in both
        sel = idx.get_flags(idx.flags.indexed) & ground_truth.get_flags(
            ground_truth.flags.indexed
        )
        idx = idx.select(sel)
        ground_truth = ground_truth.select(sel)

        # Compare hkl
        idx_hkl = idx["miller_index"].as_vec3_double()
        gt_hkl = ground_truth["miller_index"].as_vec3_double()

        # Count number of reflections indexed identically, and the number with
        # inverted indices (indexed as the Friedel pairs).
        n_direct = (idx_hkl - gt_hkl).norms().count(0.0)
        n_inverse = (idx_hkl + gt_hkl).norms().count(0.0)

        # Return the fraction of correctly indexed reflections
        return max(n_direct, n_inverse) / len(ground_truth)

    def process(self, image):
        print(f"processing {image}")
        cmd = (shutil.which("dials.import"), image) + tuple(self.import_cmds)
        result = subprocess.run(cmd)
        cmd = (shutil.which("dials.find_spots"), "imported.expt") + tuple(
            self.spotfind_cmds
        )
        result = subprocess.run(cmd)
        serialno = image.replace(self.prefix, "").replace("." + self.extension, "")
        cmd = (
            shutil.which("dials.index"),
            f"experiments_{serialno}.json",
            "strong.refl",
            "detector.fix=all",
            "output.experiments=ground_truth.expt",
            "output.reflections=ground_truth.refl",
        )
        if os.path.isfile("restraint.phil"):
            cmd = cmd + ("restraint.phil",)
        result = subprocess.run(cmd)
        strong = flex.reflection_table.from_file("strong.refl")
        d = {"image": image, "nspots": len(strong)}

        nindexed = []
        offset_deg = []
        elapsed = []
        correct_hkl = []
        for cmd in self.cmds:
            start = perf_counter()
            result = subprocess.run(cmd)
            elapsed.append(perf_counter() - start)
            if os.path.isfile("indexed.refl"):
                offset_deg.append(self.reindex_to_ground_truth(image))
                idx = flex.reflection_table.from_file("indexed.refl")
                nindexed.append(idx.get_flags(idx.flags.indexed).count(True))
                correct_hkl.append(self.check_hkl(idx))
                os.remove("indexed.refl")
                os.remove("indexed.expt")
            else:
                nindexed.append(0)
                offset_deg.append(None)
                correct_hkl.append(0)
        try:
            os.remove("ground_truth.refl")
            os.remove("ground_truth.expt")
        except FileNotFoundError:
            pass

        d["nindexed"] = nindexed
        d["offset_deg"] = offset_deg
        d["elapsed"] = elapsed
        d["correct_hkl"] = correct_hkl

        print(
            "#Spots:{0} ".format(d["nspots"])
            + " ".join([str(i + 1) + ":" + str(e) for i, e in enumerate(nindexed)])
        )

        return d


if __name__ == "__main__":

    # Edit prefix, extension and spotfinding commands as needed
    if sys.argv[1] == "TPB":
        script = Script(
            prefix="TPB_25032024_07_",
            extension="cbf",
            unit_cell="7.5840,11.225,19.711,90.0,90.0,90.0",
            space_group="P222",
            import_cmds=["panel.gain=1.7", "convert_sequences_to_stills=True"],
            spotfind_cmds=[
                "d_max=10",
            ],
        )
    elif sys.argv[1] == "lyso":
        script = Script(
            prefix="as_cbf_",
            extension="cbf",
            unit_cell="25.909,30.9595,33.3054,88.229,71.351,67.848",
            space_group="P1",
            import_cmds=["convert_sequences_to_stills=True"],
            spotfind_cmds=[
                "d_max=10",
            ],
        )
    elif sys.argv[1] == "simED":
        script = Script(
            prefix="noiseimage_",
            extension="img",
            unit_cell="78.840,78.840,38.290,90.000,90.000,90.000",
            space_group="P43212",
            import_cmds=[],
            spotfind_cmds=[
                "d_max=10",
            ],
        )
    script.run()
