import sys
import numpy as np
import collections
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons, CheckButtons, TextBox
from dials.array_family import flex
from dials.util.filter_reflections import filter_reflection_table
from dxtbx.serialize import load
import matplotlib
from orderedset import OrderedSet

class Model(object):

    def __init__(self, refls, expts):

        dose = flex.ceil(refls["xyzobs.px.value"].parts()[2])
        refls["dose"] = dose.iround()

        scaled_refls = filter_reflection_table(refls, intensity_choice=["scale"], partiality_threshold=0.4)
        intensity_s = scaled_refls["intensity.scale.value"]
        sigma_s = scaled_refls["intensity.scale.variance"] ** 0.5
        dose_s = scaled_refls["dose"]

        from cctbx import miller, crystal, uctbx

        def map_indices_to_asu(miller_indices, space_group, uc):
            """Map the indices to the asymmetric unit."""
            crystal_symmetry = crystal.symmetry(
                space_group=space_group,
                unit_cell=uc,
            )
            miller_set = miller.set(
                crystal_symmetry=crystal_symmetry, indices=miller_indices, anomalous_flag=False
            )
            miller_set_in_asu = miller_set.map_to_asu()
            return miller_set_in_asu.indices(), miller_set_in_asu.d_spacings().data()

        sg = expts[0].crystal.get_space_group()
        uc = expts[0].crystal.get_unit_cell()

        asu_index_s, d_s = map_indices_to_asu(scaled_refls["miller_index"], sg, uc)

        integrated_refls = filter_reflection_table(refls, intensity_choice=["profile"], partiality_threshold=0.4)
        intensity_u = integrated_refls["intensity.prf.value"]
        sigma_u = integrated_refls["intensity.prf.variance"] ** 0.5
        dose_u = integrated_refls["dose"]

        asu_index_u, d_u = map_indices_to_asu(integrated_refls["miller_index"], sg, uc)

        isel_s = flex.sort_permutation(d_s, reverse=True)
        sorted_asu_s = asu_index_s.select(isel_s)
        scaled_groups = list(OrderedSet(sorted_asu_s))

        # use sorted indices
        crystal_symmetry = crystal.symmetry(space_group=sg, unit_cell=uc, )
        miller_set = miller.set(
            crystal_symmetry=crystal_symmetry, indices=flex.miller_index(scaled_groups), anomalous_flag=False
        )
        d_spacings = miller_set.d_spacings().data()

        isel_u = flex.sort_permutation(d_u, reverse=True)
        sorted_asu_u = asu_index_u.select(isel_u)


        sel = flex.bool(len(sorted_asu_u), True)
        current_idx = sorted_asu_u[0]
        for i, v in enumerate(sorted_asu_u):
            if v == current_idx:
                continue
            elif v not in scaled_groups:
                sel[i] = False
            else:
                current_idx = v

        unscaled_groups = list(OrderedSet(sorted_asu_u.select(sel)))

        assert len(unscaled_groups) == len(scaled_groups)

        Data = collections.namedtuple("Data",
            ["intensity", "sigma", "dose", "asu_index"])

        self.scaled_data = Data(
            intensity=intensity_s.select(isel_s),
            sigma=sigma_s.select(isel_s),
            dose=dose_s.select(isel_s),
            asu_index=sorted_asu_s,
        )
        self.unscaled_data = Data(
            intensity=intensity_u.select(isel_u).select(sel),
            sigma=sigma_u.select(isel_u).select(sel),
            dose=dose_u.select(isel_u).select(sel),
            asu_index=sorted_asu_u.select(sel),
        )

        self.scaled_groups = scaled_groups
        self.d_spacings = d_spacings
        self.visibility = [True, False]
        self.current_miller_index = self.scaled_groups[0]

    def get_scaled_data(self, idx=0):
        miller_idx = self.scaled_groups[idx]
        print(miller_idx)
        sel = self.scaled_data.asu_index == miller_idx
        I = self.scaled_data.intensity.select(sel)
        s = self.scaled_data.sigma.select(sel)
        d = self.scaled_data.dose.select(sel)
        return d, I, s, self.visibility[0]

    def get_unscaled_data(self, idx=0):
        miller_idx = self.scaled_groups[idx]
        sel = self.unscaled_data.asu_index == miller_idx
        I = self.unscaled_data.intensity.select(sel)
        s = self.unscaled_data.sigma.select(sel)
        d = self.unscaled_data.dose.select(sel)
        return d, I, s, self.visibility[1]

    def set_visible(self, label):
        if label == "unscaled":
            self.visibility[1] = not self.visibility[1]
        elif label == "scaled":
            self.visibility[0] = not self.visibility[0]

class Viewer(object):

    def __init__(self, data):
        fig, self.ax = plt.subplots()
        plt.subplots_adjust(left=0.3, bottom=0.40)

        self.fig = fig
        self.data = data
        self.plot_main_chart()
        self.messages = []

        axcolor = 'blanchedalmond'
        axgroupno = plt.axes([0.05, 0.1, 0.70, 0.06], facecolor=axcolor, title='Resolution slider')

        msg = axgroupno.text((len(self.data.scaled_groups) // 2) * 0.85, 0.25, "(click me)")
        self.messages.append(msg)
        self.gno = Slider(
            axgroupno, '', 0, len(data.scaled_groups)-1, valinit=0, valstep=1)

        self.gno.valtext.set_text('')

        self.gno.on_changed(self.update)

        rax = plt.axes([0.05, 0.75, 0.15, 0.10], facecolor=axcolor)
        self.check = CheckButtons(rax, ('scaled', 'unscaled'), self.data.visibility)
        self.check.on_clicked(self.change_type)

        axbox = plt.axes([0.8, 0.1, 0.15, 0.06], title='Go to Miller index', facecolor=axcolor)
        self.text_box = TextBox(axbox, '',
            initial=str(self.data.scaled_groups[0]), color=axcolor)
        self.text_box.on_submit(self.submit)

        plt.show()

    def plot_main_chart(self, n=0):
        self.ax.set_title("Scaled intensity explorer")
        x, y, err, vis = self.data.get_scaled_data(n)
        if vis:
            self.ax.errorbar(x, y, yerr=err, fmt="o", visible=vis, color='k', label="scaled")
        x1, y1, err1, vis = self.data.get_unscaled_data(n)
        if vis:
            self.ax.errorbar(x1, y1, yerr=err1, fmt="o", visible=vis, color='b', label="unscaled")
        self.ax.margins(x=0.1)
        self.ax.set_xlabel("Image number/dose")
        self.ax.set_ylabel("Intensity")
        self.ax.legend(loc="best")

    def change_type(self, label):
        self.data.set_visible(label)
        idx = self.data.scaled_groups.index(self.data.current_miller_index)
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
            str(self.data.scaled_groups[n]) +", d = %.2f" % self.data.d_spacings[n],
            color='g')
        self.gno.valtext.set_text('')
        self.text_box.text_disp.set_color('k')
        # remove any help messages
        for _ in range(0, len(self.messages)):
            m = self.messages.pop()
            m.remove()

    def update_via_textbox(self, val):
        n = int(val)
        self.clear()
        self.plot_main_chart(n)
        self.data.current_miller_index = self.data.scaled_groups[n]
        self.gno.ax.xaxis.label.set_color('k')
        self.fig.canvas.draw_idle()

    def submit(self, text):
        miller_index = eval(text)
        try:
            idx = self.data.scaled_groups.index(miller_index)
        except ValueError:
            self.text_box.text_disp.set_color('r')
        else:
            self.text_box.text_disp.set_color('g')
            self.update_via_textbox(idx)



def run_viewer(args):

    refls = flex.reflection_table.from_file(args[0])
    expts = load.experiment_list(args[1], check_format=False)
    data = Model(refls, expts)
    _ = Viewer(data)


if __name__ == "__main__":
    run_viewer(sys.argv[1:])