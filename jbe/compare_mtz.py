from iotbx import mtz
import sys
from dials.array_family import flex
from annlib_ext import AnnAdaptor as ann_adaptor
import math


def data_from_mtz(filename):

    miller_arrays = mtz.object(file_name=filename).as_miller_arrays(
        merge_equivalents=False
    )

    # Select the desired columns
    intensities = None
    batches = None
    x = None
    y = None
    z = None
    g = None
    for array in miller_arrays:
        if array.info().labels == ["I", "SIGI"]:
            intensities = array
        if array.info().labels == ["BATCH"]:
            batches = array
        if array.info().labels == ["XDET"]:
            x = array
        if array.info().labels == ["YDET"]:
            y = array
        if array.info().labels == ["ROT"]:
            z = array
        if array.info().labels == ["SCALEUSED"]:
            g = array
    if not intensities:
        raise KeyError("Intensities not found in mtz file, expected labels I, SIGI")
    if not batches:
        raise KeyError("Batch values not found")
    if not x or not y:
        raise KeyError("Did not find xdet and ydet")
    if not z:
        raise KeyError("Did not find ROT in table")
    if batches.data().size() != intensities.data().size():
        raise ValueError("Batch and intensity array sizes do not match")

    # Get the unit cell and space group
    unit_cell = intensities.unit_cell()
    space_group = intensities.crystal_symmetry().space_group()

    # The reflection data
    table = flex.reflection_table()
    table["miller_index"] = intensities.indices()
    table["d"] = intensities.d_spacings().data()
    table["intensity"] = intensities.data()
    table["variance"] = intensities.sigmas() ** 2
    table["x"] = x.data()
    table["y"] = y.data()
    # scale z to match number of batches
    table["z"] = z.data()
    table["batch"] = batches.data()
    if g:
        table["scale_factor"] = g.data()
    return table


max_allowable_distance = 15
max_z_diff = 7
dmin = 12


def report_on_non_matches(table1, table2, dmin=5.0):
    # First set up binning
    sel = table1["match_idx"] == -1
    dsel = table1["d"] > dmin
    no_matches = table1.select(sel & dsel)
    print("Reflections in table1 and not table2, with d > %s:" % dmin)
    print("Miller index, d, intensity, sigma, xyz")
    for i in range(0, no_matches.size()):
        print(
            "%s, %.3f, %.4f, %.4f, (%.2f, %.2f, %.2f)"
            % (
                no_matches["miller_index"][i],
                no_matches["d"][i],
                no_matches["intensity"][i],
                no_matches["variance"][i] ** 0.5,
                no_matches["x"][i],
                no_matches["y"][i],
                no_matches["z"][i],
            )
        )
    matching_indices = set(table1["match_idx"])
    all_indices = set(range(0, table2.size()))
    no_matching_indices = all_indices.difference(matching_indices)
    sel = flex.size_t(list(no_matching_indices))
    no_match = table2.select(sel)
    dsel = no_match["d"] > dmin
    no_match = no_match.select(dsel)
    print("Reflections in table2 and not table1, with d > %s:" % dmin)
    print("Miller index, d, intensity, sigma, xyz")
    for i in range(0, no_match.size()):
        print(
            "%s, %.3f, %.4f, %.4f, (%.2f, %.2f, %.2f)"
            % (
                no_match["miller_index"][i],
                no_match["d"][i],
                no_match["intensity"][i],
                no_match["variance"][i] ** 0.5,
                no_match["x"][i],
                no_match["y"][i],
                no_match["z"][i],
            )
        )


def plot_scales(table1, table2):
    import matplotlib.pyplot as plt

    sel = table1["match_idx"] != -1
    intensities_1 = table1["intensity"].select(sel)
    indices = table1["match_idx"].select(sel)
    isel = flex.size_t()
    for i in indices:
        isel.append(i)
    intensities_2 = table2["intensity"].select(isel)
    sel = (intensities_1 < 1) | (intensities_2 < 1)
    intensities_1 = intensities_1.select(~sel)
    intensities_2 = intensities_2.select(~sel)

    plt.figure()
    ax = plt.gca()
    ax.scatter(intensities_1, intensities_2, s=0.5)
    ax.grid()
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Intensity, dataset 1")
    ax.set_ylabel("Intensity, dataset 2")
    plt.show()


def run(args):
    table_1 = data_from_mtz(args[0])
    table_2 = data_from_mtz(args[1])

    zr1 = max(table_1["z"]) - min(table_1["z"])
    zr2 = max(table_2["z"]) - min(table_2["z"])
    if zr1 / zr2 > 1.2:
        table_2["z"] *= int(zr1 / zr2)
        print(
            "scaled ROT values of dataset 2 by %s to match dataset 1" % int(zr1 / zr2)
        )
    elif zr2 / zr1 > 1.2:
        table_1["z"] *= int(zr2 / zr1)
        print(
            "scaled ROT values of dataset 1 by %s to match dataset 2" % int(zr2 / zr1)
        )

    # print(br1, br2)
    # print(max(br1/br2, br2/br1))

    table_1["match_idx"] = flex.int(table_1.size(), -1)

    xy0 = flex.vec2_double(table_1["x"], table_1["y"])
    # add 1e-6 as if the values are exactly the same then not found as nearest
    # neighbour! e.g. trying to compare dials vs aimless
    xy1 = flex.vec2_double(table_2["x"] + 1e-2, table_2["y"])

    z0 = table_1["z"]
    z1 = table_2["z"]

    zs = range(int(flex.min(z0)), int(flex.max(z0)) + 1)
    matches = 0
    for z in zs:
        sel = (z0 < z + max_z_diff) & (z0 > z - max_z_diff)
        xy = xy0.select(sel)
        table1_indices = sel.iselection()
        ann = ann_adaptor(xy.as_double().as_1d(), 2)
        sel2 = (z1 < z + max_z_diff) & (z1 > z - max_z_diff)
        xy = xy1.select(sel2)
        table2_indices = sel2.iselection()
        xy1d = xy.as_double().as_1d()
        ann.query(xy1d)
        for i, nn_idx in enumerate(
            ann.nn
        ):  # nn_idx is index of table1, i index of table 2
            if math.sqrt(ann.distances[i]) < max_allowable_distance:
                table_1["match_idx"][table1_indices[nn_idx]] = table2_indices[i]
                matches += 1

    print("Table 1 size:%s" % table_1.size())
    print("Table 2 size:%s" % table_2.size())
    print(
        "Found %s matches (searching for matches to table1 in table2)"
        % str((table_1["match_idx"] != -1).count(True))
    )
    indices = set(table_1["match_idx"])
    n_unique = len(indices)
    if -1 in indices:
        n_unique -= 1
    print("%s unique matches" % str(n_unique))

    correctly_matched = 0
    incorrectly_matched = 0
    for i in range(len(table_1)):
        if table_1["match_idx"][i] != -1:
            if (
                table_1["miller_index"][i]
                == table_2["miller_index"][table_1["match_idx"][i]]
            ):
                correctly_matched += 1
            else:
                # pass
                incorrectly_matched += 1
                # print("Incorrectly matched %s, %s" % (table_1["miller_index"][i], table_2["miller_index"][table_1["match_idx"][i]]))
                table_1["match_idx"][i] = -1
    print("N correctly matched %s" % correctly_matched)
    print("N incorrectly matched %s" % incorrectly_matched)
    report_on_non_matches(table_1, table_2, dmin)
    plot_scales(table_1, table_2)


if __name__ == "__main__":
    run(sys.argv[1:])
