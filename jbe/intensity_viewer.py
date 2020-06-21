"""Experimental simple viewer to plot intensities of symmetry equivalents.

Requires a scaled reflection table and experiments file (can be multi-dataset)

Plots as a function of image number/"dose".
"""
import sys
import os
import collections
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, CheckButtons, TextBox
from dials.util.options import OptionParser, reflections_and_experiments_from_files
from dials.array_family import flex
from dials.util.filter_reflections import (
    SumAndPrfIntensityReducer,
    ScaleIntensityReducer,
)
from orderedset import OrderedSet
from cctbx import miller, crystal, uctbx, sgtbx
from iotbx.reflection_file_reader import any_reflection_file
from iotbx import mtz


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
    ipr = None
    for array in miller_arrays:
        if array.info().labels == ["I", "SIGI"]:
            intensities = array.as_anomalous_array()
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
        if array.info().labels == ["IPR", "SIGIPR"]:
            ipr = array
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
    table["sigma"] = intensities.sigmas()
    if ipr:
        table["intensity.prf.value"] = ipr.data()
        table["intensity.prf.sigma"] = ipr.sigmas()
    table["x"] = x.data()
    table["y"] = y.data()
    # scale z to match number of batches
    table["z"] = z.data()
    table["batch"] = batches.data()
    if g:
        table["scale_factor"] = g.data()
    return table, unit_cell, space_group


def read_xds_ascii(filename):

    hkl = flex.miller_index()
    I = flex.double()
    sigma = flex.double()
    xyz = flex.vec3_double()
    uc = None
    sg = None

    with open(filename, "r") as f:
        for line in f.readlines():
            vals = line.split()
            if line[0] != "!" and len(vals) == 12:
                hkl.append((int(vals[0]), int(vals[1]), int(vals[2])))
                I.append(float(vals[3]))
                sigma.append(float(vals[4]))
                xyz.append((float(vals[5]), float(vals[6]), float(vals[7])))
                # assert
            elif line.startswith("!UNIT_CELL_CONSTANTS="):
                vals = line.split()
                uc = uctbx.unit_cell(
                    parameters=(
                        float(vals[1]),
                        float(vals[2]),
                        float(vals[3]),
                        float(vals[4]),
                        float(vals[5]),
                        float(vals[6]),
                    )
                )
            elif line.startswith("!SPACE_GROUP_NUMBER="):
                sg = sgtbx.space_group_info(number=line.split()[1]).group()

    data = flex.reflection_table()
    data["intensity"] = I
    data["sigma"] = sigma
    data["miller_index"] = hkl
    x, y, z = xyz.parts()
    data["x"] = x
    data["y"] = y
    data["z"] = z
    return data, uc, sg


def map_indices_to_asu(miller_indices, space_group, uc, anom=False):
    """Map the indices to the asymmetric unit."""
    assert anom in (False, True)
    crystal_symmetry = crystal.symmetry(space_group=space_group, unit_cell=uc,)
    miller_set = miller.set(
        crystal_symmetry=crystal_symmetry, indices=miller_indices, anomalous_flag=anom
    )
    miller_set_in_asu = miller_set.map_to_asu()
    return miller_set_in_asu.indices(), miller_set_in_asu.d_spacings().data()


def setup_dials_models(refls, expts):

    ### Initial setting up of data
    dose = flex.ceil(refls["xyzobs.px.value"].parts()[2])
    refls["dose"] = dose.iround()

    # want to get scaled data with outliers still present
    excluded = refls.get_flags(refls.flags.excluded_for_scaling) | refls.get_flags(
        refls.flags.user_excluded_in_scaling
    )
    good = ~excluded
    refls = refls.select(good)  # still has outliers
    refls = ScaleIntensityReducer.apply_scaling_factors(refls)
    scaled_refls = SumAndPrfIntensityReducer.apply_scaling_factors(refls)

    sg = expts[0].crystal.get_space_group()
    uc = expts[0].crystal.get_unit_cell()

    asu_index_s, d_s = map_indices_to_asu(scaled_refls["miller_index"], sg, uc)
    anom_index_s, _ = map_indices_to_asu(
        scaled_refls["miller_index"], sg, uc, anom=True
    )
    ### finished initial setting up, now know indices

    ### Get relevant data and sort by asu index
    intensity_s = scaled_refls["intensity.scale.value"]
    sigma_s = scaled_refls["intensity.scale.variance"] ** 0.5
    dose_s = scaled_refls["dose"]

    integrated_refls = scaled_refls
    intensity_u = integrated_refls["intensity.prf.value"]
    sigma_u = flex.sqrt(integrated_refls["intensity.prf.variance"])

    isel_s = flex.sort_permutation(d_s, reverse=True)
    sorted_asu_s = asu_index_s.select(isel_s)
    sorted_anom_s = anom_index_s.select(isel_s)
    scaled_groups = list(OrderedSet(sorted_asu_s))
    outliers = scaled_refls.get_flags(scaled_refls.flags.outlier_in_scaling)

    Data = collections.namedtuple(
        "Data", ["intensity", "sigma", "dose", "asu_index", "anom_index", "outliers"]
    )

    data = Data(
        intensity=intensity_s.select(isel_s),
        sigma=sigma_s.select(isel_s),
        dose=dose_s.select(isel_s),
        asu_index=sorted_asu_s,
        anom_index=sorted_anom_s,
        outliers=outliers.select(isel_s),
    )
    dials_scaled = DialsScaledModel(data)

    unscaled_data = Data(
        intensity=intensity_u.select(isel_s),
        sigma=sigma_u.select(isel_s),
        dose=dose_s.select(isel_s),
        asu_index=sorted_asu_s,
        anom_index=sorted_anom_s,
        outliers=None,
    )
    dials_unscaled = DialsUnScaledModel(unscaled_data)
    return [dials_scaled, dials_unscaled], scaled_groups, uc, sg


class DialsScaledModel(object):
    def __init__(self, data, label="DIALS\nscaled"):
        self.scaled_data = data
        self.label = label

    def get_data(self, miller_idx):
        sel = self.scaled_data.asu_index == miller_idx
        I = self.scaled_data.intensity.select(sel)
        if not I:
            return None
        s = self.scaled_data.sigma.select(sel)
        d = self.scaled_data.dose.select(sel)
        if self.scaled_data.outliers:
            outliers = self.scaled_data.outliers.select(sel)
        else:
            outliers = None
        pairs = self.scaled_data.anom_index.select(sel)
        anom = pairs == pairs[0]
        return d, I, s, anom, outliers

    def add_to_plot(self, ax, miller_idx):
        result = self.get_data(miller_idx)
        if not result:
            return
        x, y, err, anom, outliers = result

        n_outliers = 0
        if outliers:
            n_outliers = outliers.count(True)
        if n_outliers:
            sel1 = ~outliers & anom
            sel2 = ~outliers & ~anom
        else:
            sel1 = anom
            sel2 = ~anom
        ax.errorbar(
            x.select(sel1),
            y.select(sel1),
            yerr=err.select(sel1),
            fmt="o",
            visible=True,
            color="k",
            label=self.label,
        )
        ax.errorbar(
            x.select(sel2),
            y.select(sel2),
            yerr=err.select(sel2),
            fmt="v",
            visible=True,
            color="k",
        )
        if n_outliers:
            sel3 = outliers & anom
            sel4 = outliers & ~anom
            ax.errorbar(
                x.select(sel3),
                y.select(sel3),
                yerr=err.select(sel3),
                fmt="o",
                visible=True,
                color="r",
                label=self.label + "\noutlier",
            )
            ax.errorbar(
                x.select(sel4),
                y.select(sel4),
                yerr=err.select(sel4),
                fmt="v",
                visible=True,
                color="r",
            )


class DialsUnScaledModel(object):
    def __init__(self, data, label="DIALS\nunscaled"):
        self.unscaled_data = data
        self.label = label

    def get_data(self, miller_idx):
        sel = self.unscaled_data.asu_index == miller_idx
        I = self.unscaled_data.intensity.select(sel)
        if not I:
            return None
        s = self.unscaled_data.sigma.select(sel)
        d = self.unscaled_data.dose.select(sel)
        pairs = self.unscaled_data.anom_index.select(sel)
        anom = pairs == pairs[0]
        return d, I, s, anom

    def add_to_plot(self, ax, miller_idx):
        result = self.get_data(miller_idx)
        if not result:
            return
        x1, y1, err1, anom = result
        neg = ~anom
        ax.errorbar(
            x1.select(anom),
            y1.select(anom),
            yerr=err1.select(anom),
            fmt="o",
            visible=True,
            color="b",
            label=self.label,
        )
        ax.errorbar(
            x1.select(neg),
            y1.select(neg),
            yerr=err1.select(neg),
            fmt="v",
            visible=True,
            color="b",
        )


class XDSModel(object):
    def __init__(self, data, label="XDS\nscaled"):
        self.scaled_data = data
        self.label = label

    def get_data(self, miller_idx):
        sel = self.scaled_data.asu_index == miller_idx
        I = self.scaled_data.intensity.select(sel)
        if not I:
            return None
        s = self.scaled_data.sigma.select(sel)
        d = self.scaled_data.dose.select(sel)
        pairs = self.scaled_data.anom_index.select(sel)
        anom = pairs == pairs[0]
        return d, I, s, anom

    def add_to_plot(self, ax, miller_idx):
        result = self.get_data(miller_idx)
        if not result:
            return
        x, y, err, anom = result
        sel1 = anom
        sel2 = ~anom
        ax.errorbar(
            x.select(sel1),
            y.select(sel1),
            yerr=err.select(sel1),
            fmt="o",
            visible=True,
            color="g",
            label=self.label,
        )
        ax.errorbar(
            x.select(sel2),
            y.select(sel2),
            yerr=err.select(sel2),
            fmt="v",
            visible=True,
            color="g",
        )


def setup_xds_models(xdsasciifile):
    table, uc, sg = read_xds_ascii(xdsasciifile)

    asu_index_s, d_s = map_indices_to_asu(table["miller_index"], sg, uc)
    anom_index_s, _ = map_indices_to_asu(table["miller_index"], sg, uc, anom=True)

    intensity_s = table["intensity"]
    sigma_s = table["sigma"]

    isel_s = flex.sort_permutation(d_s, reverse=True)
    sorted_asu_s = asu_index_s.select(isel_s)
    sorted_anom_s = anom_index_s.select(isel_s)
    scaled_groups = list(OrderedSet(sorted_asu_s))

    Data = collections.namedtuple(
        "Data", ["intensity", "sigma", "dose", "asu_index", "anom_index"]
    )

    data = Data(
        intensity=intensity_s.select(isel_s),
        sigma=sigma_s.select(isel_s),
        dose=table["z"].select(isel_s),
        asu_index=sorted_asu_s,
        anom_index=sorted_anom_s,
    )

    return [XDSModel(data)], scaled_groups, uc, sg


class MTZModel(object):
    def __init__(self, data, label="MTZ\nscaled"):
        self.scaled_data = data
        self.label = label

    def get_data(self, miller_idx):
        sel = self.scaled_data.asu_index == miller_idx
        I = self.scaled_data.intensity.select(sel)
        if not I:
            return None
        s = self.scaled_data.sigma.select(sel)
        d = self.scaled_data.dose.select(sel)
        pairs = self.scaled_data.anom_index.select(sel)
        anom = pairs == pairs[0]
        return d, I, s, anom

    def add_to_plot(self, ax, miller_idx):
        result = self.get_data(miller_idx)
        if not result:
            return

        x, y, err, anom = result

        sel1 = anom
        sel2 = ~anom
        ax.errorbar(
            x.select(sel1),
            y.select(sel1),
            yerr=err.select(sel1),
            fmt="o",
            visible=True,
            color="darkorange",
            label=self.label,
        )
        ax.errorbar(
            x.select(sel2),
            y.select(sel2),
            yerr=err.select(sel2),
            fmt="v",
            visible=True,
            color="darkorange",
        )


def setup_mtz_models(mtzfile):

    table, uc, sg = data_from_mtz(mtzfile)
    asu_index_s, d_s = map_indices_to_asu(table["miller_index"], sg, uc)
    anom_index_s, _ = map_indices_to_asu(table["miller_index"], sg, uc, anom=True)

    asu_index_s, d_s = map_indices_to_asu(table["miller_index"], sg, uc)
    anom_index_s, _ = map_indices_to_asu(table["miller_index"], sg, uc, anom=True)

    intensity_s = table["intensity"]
    sigma_s = table["sigma"]

    isel_s = flex.sort_permutation(d_s, reverse=True)
    sorted_asu_s = asu_index_s.select(isel_s)
    sorted_anom_s = anom_index_s.select(isel_s)
    scaled_groups = list(OrderedSet(sorted_asu_s))

    Data = collections.namedtuple(
        "Data", ["intensity", "sigma", "dose", "asu_index", "anom_index"]
    )

    data = Data(
        intensity=intensity_s.select(isel_s),
        sigma=sigma_s.select(isel_s),
        dose=table["batch"].select(isel_s),
        asu_index=sorted_asu_s,
        anom_index=sorted_anom_s,
    )

    return [MTZModel(data)], scaled_groups, uc, sg


class Model(object):
    def __init__(self, dataseries, scaled_groups, d_spacings):

        self.data_series = dataseries
        self.scaled_groups = flex.miller_index(scaled_groups)
        self.d_spacings = d_spacings
        self.current_miller_index = self.scaled_groups[0]


class Viewer(object):
    def __init__(self, data):
        self.data = data
        if len(self.data.data_series) > 1:
            fig, self.ax = plt.subplots(len(self.data.data_series), 1, sharex=True)
        else:
            fig, self.ax = plt.subplots()
            self.ax = [self.ax]

        self.visibilities = []
        self.labels = None
        if len(self.data.data_series[0]) > 1:
            self.labels = [d.label for d in self.data.data_series[0]]
        for d in self.data.data_series:
            self.visibilities.append([True] + ([False] * (len(d) - 1)))

        plt.subplots_adjust(left=0.3, bottom=0.3)

        self.fig = fig

        self.plot_main_chart()
        self.messages = []

        axcolor = "blanchedalmond"
        axgroupno = plt.axes(
            [0.05, 0.05, 0.70, 0.06], facecolor=axcolor, title="Resolution slider"
        )

        msg = axgroupno.text(
            (len(self.data.scaled_groups) // 2) * 0.85, 0.25, "(click me)"
        )
        self.messages.append(msg)
        self.gno = Slider(
            axgroupno, "", 0, len(data.scaled_groups) - 1, valinit=0, valstep=1
        )

        self.gno.valtext.set_text("")

        self.gno.on_changed(self.update)
        if len(self.data.data_series[0]) > 1:  # i.e. if we have DIALS data.
            rax = plt.axes([0.05, 0.7, 0.15, 0.2], facecolor=axcolor)
            self.check = CheckButtons(
                rax,
                list(c.label for c in self.data.data_series[0]),
                self.visibilities[0],
            )
            self.check.on_clicked(self.change_type)

        axbox = plt.axes(
            [0.8, 0.05, 0.15, 0.06], title="Go to Miller index", facecolor=axcolor
        )
        self.text_box = TextBox(
            axbox, "", initial=str(self.data.scaled_groups[0]), color=axcolor
        )
        self.text_box.on_submit(self.submit)

        plt.show()

    def set_visible(self, label):
        idx = self.labels.index(label)
        self.visibilities[0][idx] = not self.visibilities[0][idx]

    def plot_main_chart(self, n=0):
        # self.ax.set_title("Scaled intensity explorer")

        for vis, d, ax in zip(self.visibilities, self.data.data_series, self.ax):
            miller_idx = self.data.scaled_groups[n]
            for v, subitem in zip(vis, d):
                if v:
                    subitem.add_to_plot(ax, miller_idx)
            if True in vis:
                ax.legend(loc="best")

        self.ax[-1].set_xlabel("Image number/dose")
        for ax in self.ax:
            ax.margins(x=0.03)
            ax.set_ylabel("Intensity")

    def change_type(self, label):
        self.set_visible(label)
        idx = (self.data.scaled_groups == self.data.current_miller_index).iselection()[
            0
        ]
        self.update_via_textbox(idx)
        self.fig.canvas.draw_idle()

    def clear(self):
        for ax in self.ax:
            ax.cla()

    def update(self, val):
        n = int(val)
        self.clear()
        self.plot_main_chart(n)
        self.fig.canvas.draw_idle()
        self.data.current_miller_index = self.data.scaled_groups[n]
        self.gno.ax.set_xlabel(
            str(self.data.scaled_groups[n]) + ", d = %.2f" % self.data.d_spacings[n],
            color="g",
        )
        self.gno.valtext.set_text("")
        self.text_box.text_disp.set_color("k")
        # remove any help messages
        for _ in range(0, len(self.messages)):
            m = self.messages.pop()
            m.remove()

    def update_via_textbox(self, val):
        n = int(val)
        self.clear()
        self.plot_main_chart(n)
        self.data.current_miller_index = self.data.scaled_groups[n]
        self.gno.ax.xaxis.label.set_color("k")
        self.fig.canvas.draw_idle()

    def submit(self, text):
        miller_index = eval(text)
        try:
            idx = (self.data.scaled_groups == miller_index).iselection()[0]
        except IndexError:
            self.text_box.text_disp.set_color("r")
        else:
            self.text_box.text_disp.set_color("g")
            self.update_via_textbox(idx)


def run_viewer():
    """Run the intensity explorer"""

    parser = OptionParser(
        read_experiments=True,
        read_reflections=True,
        check_format=False,
        epilog=__doc__,
    )
    params, _, args = parser.parse_args(show_diff_phil=False, return_unhandled=True)

    data_series = []
    uclist = []
    sglist = []
    allgroups = flex.miller_index()

    reflections, experiments = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )
    if reflections and experiments:
        dials_dataseries, groups, uc, sg = setup_dials_models(
            reflections[0], experiments
        )
        data_series.append(dials_dataseries)
        uclist.append(uc)
        sglist.append(sg)
        allgroups.extend(flex.miller_index(groups))
    if args:
        for arg in args:
            print("trying to read %s as an MTZ or XDS ascii file" % arg)
            assert os.path.isfile(arg), arg
            reader = any_reflection_file(arg)
            if reader.file_type() == "ccp4_mtz":
                mtz_dataseries, groups, uc, sg = setup_mtz_models(arg)
                data_series.append(mtz_dataseries)
                uclist.append(uc)
                sglist.append(sg)
                allgroups.extend(flex.miller_index(groups))
            else:
                try:
                    xds_dataseries, groups, uc, sg = setup_xds_models(arg)
                except Exception as e:
                    pass
                else:
                    data_series.append(xds_dataseries)
                    uclist.append(uc)
                    sglist.append(sg)
                    allgroups.extend(flex.miller_index(groups))

    if not data_series:
        parser.print_help()
        sys.exit()
    # sort all groups and calc d_spacings
    # generate final dpsacings from overall miller set and uc, sg
    sgs = [s.type().number() for s in sglist]
    if len(set(sgs)) > 1:
        raise ValueError(
            "Spacegroups not equal, found space group numbers: %s"
            % ",".join(str(s) for s in set(sgs))
        )
    for uc in uclist[1:]:
        if not uclist[0].is_similar_to(uc):
            raise ValueError(
                """
Unit cells must be similar to compare multiple datasets.
Unit cells not similar:
%s
%s
"""
                % (str(uclist[0]), str(uc))
            )

    crystal_symmetry = crystal.symmetry(space_group=sglist[0], unit_cell=uclist[0])
    miller_set = miller.set(
        crystal_symmetry=crystal_symmetry, indices=allgroups, anomalous_flag=False
    )
    ds = miller_set.d_spacings().data()
    indices = miller_set.indices()

    isel_s = flex.sort_permutation(ds, reverse=True)
    sorted_indices = indices.select(isel_s)
    d_spacings = ds.select(isel_s)
    scaled_groups = list(OrderedSet(sorted_indices))

    data = Model(data_series, scaled_groups, d_spacings)
    _ = Viewer(data)


if __name__ == "__main__":
    run_viewer()
