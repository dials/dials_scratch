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
from dials.util.filter_reflections import SumAndPrfIntensityReducer, ScaleIntensityReducer
from orderedset import OrderedSet
from cctbx import miller, crystal
from iotbx.reflection_file_reader import any_reflection_file
from dials_scratch.jbe.compare_dials_xds import read_xds_ascii, data_from_mtz

def map_indices_to_asu(miller_indices, space_group, uc, anom=False):
    """Map the indices to the asymmetric unit."""
    assert anom in (False, True)
    crystal_symmetry = crystal.symmetry(space_group=space_group, unit_cell=uc,)
    miller_set = miller.set(
        crystal_symmetry=crystal_symmetry, indices=miller_indices, anomalous_flag=anom
    )
    miller_set_in_asu = miller_set.map_to_asu()
    return miller_set_in_asu.indices(), miller_set_in_asu.d_spacings().data()

class Model(object):
    def __init__(self, refls, expts):

        dose = flex.ceil(refls["xyzobs.px.value"].parts()[2])
        refls["dose"] = dose.iround()

        # want to get scaled data with outliers
        excluded = refls.get_flags(refls.flags.excluded_for_scaling) | refls.get_flags(refls.flags.user_excluded_in_scaling)
        good = ~excluded
        refls = refls.select(good) #still has outliers
        refls = ScaleIntensityReducer.apply_scaling_factors(refls)
        scaled_refls = SumAndPrfIntensityReducer.apply_scaling_factors(refls)
        outliers = scaled_refls.get_flags(scaled_refls.flags.outlier_in_scaling)

        intensity_s = scaled_refls["intensity.scale.value"]
        sigma_s = scaled_refls["intensity.scale.variance"] ** 0.5
        dose_s = scaled_refls["dose"]

        sg = expts[0].crystal.get_space_group()
        uc = expts[0].crystal.get_unit_cell()

        asu_index_s, d_s = map_indices_to_asu(scaled_refls["miller_index"], sg, uc)
        anom_index_s, _ = map_indices_to_asu(scaled_refls["miller_index"], sg, uc, anom=True)
        integrated_refls = scaled_refls
        intensity_u = integrated_refls["intensity.prf.value"]
        sigma_u = flex.sqrt(integrated_refls["intensity.prf.variance"])

        isel_s = flex.sort_permutation(d_s, reverse=True)
        sorted_asu_s = asu_index_s.select(isel_s)
        sorted_anom_s = anom_index_s.select(isel_s)
        scaled_groups = list(OrderedSet(sorted_asu_s))
        outliers_s = outliers.select(isel_s)

        # use sorted indices
        crystal_symmetry = crystal.symmetry(space_group=sg, unit_cell=uc)
        miller_set = miller.set(
            crystal_symmetry=crystal_symmetry,
            indices=flex.miller_index(scaled_groups),
            anomalous_flag=False,
        )
        d_spacings = miller_set.d_spacings().data()

        Data = collections.namedtuple(
            "Data", ["intensity", "sigma", "dose", "asu_index", "anom_index", "outliers"]
        )

        self.scaled_data = Data(
            intensity=intensity_s.select(isel_s),
            sigma=sigma_s.select(isel_s),
            dose=dose_s.select(isel_s),
            asu_index=sorted_asu_s,
            anom_index=sorted_anom_s,
            outliers=outliers_s,
        )
        self.unscaled_data = Data(
            intensity=intensity_u.select(isel_s),
            sigma=sigma_u.select(isel_s),
            dose=dose_s.select(isel_s),
            asu_index=sorted_asu_s,
            anom_index=sorted_anom_s,
            outliers=outliers_s,
        )

        self.scaled_groups = flex.miller_index(scaled_groups)
        self.d_spacings = d_spacings
        self.visibility = [True, False]
        self.current_miller_index = self.scaled_groups[0]

    def get_scaled_data(self, idx=0):
        miller_idx = self.scaled_groups[idx]
        sel = self.scaled_data.asu_index == miller_idx
        I = self.scaled_data.intensity.select(sel)
        s = self.scaled_data.sigma.select(sel)
        d = self.scaled_data.dose.select(sel)
        if self.scaled_data.outliers:
            outliers = self.scaled_data.outliers.select(sel)
        else:
            outliers = None
        pairs = self.scaled_data.anom_index.select(sel)
        anom = (pairs == pairs[0])
        return d, I, s, self.visibility[0], anom, outliers

    def get_unscaled_data(self, idx=0):
        miller_idx = self.scaled_groups[idx]
        sel = self.unscaled_data.asu_index == miller_idx
        I = self.unscaled_data.intensity.select(sel)
        s = self.unscaled_data.sigma.select(sel)
        d = self.unscaled_data.dose.select(sel)
        pairs = self.unscaled_data.anom_index.select(sel)
        anom = (pairs == pairs[0])
        return d, I, s, self.visibility[1], anom

    def set_visible(self, label):
        if label == "unscaled":
            self.visibility[1] = not self.visibility[1]
        elif label == "scaled":
            self.visibility[0] = not self.visibility[0]

class XDSModel(Model):
    def __init__(self, xdsasciifile):

        table, uc, sg = read_xds_ascii(xdsasciifile)
        self.unscaled_data = None

        asu_index_s, d_s = map_indices_to_asu(table["miller_index"], sg, uc)
        anom_index_s, _ = map_indices_to_asu(table["miller_index"], sg, uc, anom=True)

        intensity_s = table["intensity"]
        sigma_s = table["sigma"]

        isel_s = flex.sort_permutation(d_s, reverse=True)
        sorted_asu_s = asu_index_s.select(isel_s)
        sorted_anom_s = anom_index_s.select(isel_s)
        scaled_groups = list(OrderedSet(sorted_asu_s))

        # use sorted indices
        crystal_symmetry = crystal.symmetry(space_group=sg, unit_cell=uc)
        miller_set = miller.set(
            crystal_symmetry=crystal_symmetry,
            indices=flex.miller_index(scaled_groups),
            anomalous_flag=False,
        )
        d_spacings = miller_set.d_spacings().data()

        Data = collections.namedtuple(
            "Data", ["intensity", "sigma", "dose", "asu_index", "anom_index", "outliers"]
        )

        self.scaled_data = Data(
            intensity=intensity_s.select(isel_s),
            sigma=sigma_s.select(isel_s),
            dose=table["z"].select(isel_s),
            asu_index=sorted_asu_s,
            anom_index=sorted_anom_s,
            outliers=None,
        )

        self.scaled_groups = flex.miller_index(scaled_groups)
        self.d_spacings = d_spacings
        self.visibility = [True, False]
        self.current_miller_index = self.scaled_groups[0]

class MTZModel(Model):
    def __init__(self, mtzfile):
        self.unscaled_data = None
        table, uc, sg = data_from_mtz(mtzfile)
        asu_index_s, d_s = map_indices_to_asu(table["miller_index"], sg, uc)
        anom_index_s, _ = map_indices_to_asu(table["miller_index"], sg, uc, anom=True)

        intensity_s = table["intensity"]
        sigma_s = table["sigma"]
        if "intensity.prf.value" in table:
            intensity_u = table["intensity.prf.value"]
            sigma_u = table["intensity.prf.sigma"]

        isel_s = flex.sort_permutation(d_s, reverse=True)
        sorted_asu_s = asu_index_s.select(isel_s)
        sorted_anom_s = anom_index_s.select(isel_s)
        scaled_groups = list(OrderedSet(sorted_asu_s))

        # use sorted indices
        crystal_symmetry = crystal.symmetry(space_group=sg, unit_cell=uc)
        miller_set = miller.set(
            crystal_symmetry=crystal_symmetry,
            indices=flex.miller_index(scaled_groups),
            anomalous_flag=False,
        )
        d_spacings = miller_set.d_spacings().data()

        Data = collections.namedtuple(
            "Data", ["intensity", "sigma", "dose", "asu_index", "anom_index", "outliers"]
        )

        self.scaled_data = Data(
            intensity=intensity_s.select(isel_s),
            sigma=sigma_s.select(isel_s),
            dose=table["batch"].select(isel_s),
            asu_index=sorted_asu_s,
            anom_index=sorted_anom_s,
            outliers=None,
        )
        if "intensity.prf.value" in table:
            self.unscaled_data = Data(
                intensity=intensity_u.select(isel_s),
                sigma=sigma_u.select(isel_s),
                dose=table["batch"].select(isel_s),
                asu_index=sorted_asu_s,
                anom_index=sorted_anom_s,
                outliers=None,
            )

        self.scaled_groups = flex.miller_index(scaled_groups)
        self.d_spacings = d_spacings
        self.visibility = [True, False]
        self.current_miller_index = self.scaled_groups[0]

class Viewer(object):
    def __init__(self, data):
        fig, self.ax = plt.subplots()
        plt.subplots_adjust(left=0.3, bottom=0.40)

        self.fig = fig
        self.data = data
        self.plot_main_chart()
        self.messages = []

        axcolor = "blanchedalmond"
        axgroupno = plt.axes(
            [0.05, 0.1, 0.70, 0.06], facecolor=axcolor, title="Resolution slider"
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

        rax = plt.axes([0.05, 0.75, 0.15, 0.10], facecolor=axcolor)
        self.check = CheckButtons(rax, ("scaled", "unscaled"), self.data.visibility)
        self.check.on_clicked(self.change_type)

        axbox = plt.axes(
            [0.8, 0.1, 0.15, 0.06], title="Go to Miller index", facecolor=axcolor
        )
        self.text_box = TextBox(
            axbox, "", initial=str(self.data.scaled_groups[0]), color=axcolor
        )
        self.text_box.on_submit(self.submit)

        plt.show()

    def plot_main_chart(self, n=0):
        self.ax.set_title("Scaled intensity explorer")

        # For scaled data, plot outliers separately and also separate anomalous pairs
        x, y, err, vis, anom, outliers = self.data.get_scaled_data(n)
        neg = ~anom
        if vis:
            n_outliers = 0
            if outliers:
                n_outliers = outliers.count(True)
            if n_outliers:
                sel1 = ~outliers & anom
                sel2 = ~outliers & ~anom
            else:
                sel1 = anom
                sel2 = ~anom
            self.ax.errorbar(
                x.select(sel1), y.select(sel1), yerr=err.select(sel1), fmt="o", visible=vis, color="k", label="scaled"
            )
            self.ax.errorbar(
                x.select(sel2), y.select(sel2), yerr=err.select(sel2), fmt="v", visible=vis, color="k",
            )
            if n_outliers:
                sel3 = outliers & anom
                sel4 = outliers & ~anom
                self.ax.errorbar(
                    x.select(sel3), y.select(sel3), yerr=err.select(sel3), fmt="o", visible=vis, color="r", label="outlier"
                )
                self.ax.errorbar(
                    x.select(sel4), y.select(sel4), yerr=err.select(sel4), fmt="v", visible=vis, color="r",
                )
        if self.data.unscaled_data:
            x1, y1, err1, vis, anom = self.data.get_unscaled_data(n)
            neg = ~anom
            if vis:
                self.ax.errorbar(
                    x1.select(anom), y1.select(anom), yerr=err1.select(anom), fmt="o", visible=vis, color="b", label="unscaled",
                )
                self.ax.errorbar(
                    x1.select(neg), y1.select(neg), yerr=err1.select(neg), fmt="v", visible=vis, color="b",
                )
        self.ax.margins(x=0.1)
        self.ax.set_xlabel("Image number/dose")
        self.ax.set_ylabel("Intensity")
        self.ax.legend(loc="best")

    def change_type(self, label):
        self.data.set_visible(label)
        idx = (self.data.scaled_groups == self.data.current_miller_index).iselection()[
            0
        ]
        self.update_via_textbox(idx)
        self.fig.canvas.draw_idle()

    def clear(self):
        self.ax.cla()

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

    if not params.input.experiments or not params.input.reflections:
        if args:
            for arg in args:
                assert os.path.isfile(arg), arg
                reader = any_reflection_file(arg)
                if reader.file_type() == "ccp4_mtz":
                    data = MTZModel(arg)
                    break
                else:
                    try:
                        data = XDSModel(arg)
                    except Exception:
                        pass
                    else:
                        break
        else:
            parser.print_help()
            sys.exit()
    else:
        reflections, experiments = reflections_and_experiments_from_files(
            params.input.reflections, params.input.experiments
        )

        data = Model(reflections[0], experiments)
    _ = Viewer(data)


if __name__ == "__main__":
    run_viewer()
