"""Definitions of screw axes with methods for scoring against data."""
from scitbx.array_family import flex
from dials.util.observer import Observer, Subject, singleton
from jinja2 import Environment, ChoiceLoader, PackageLoader
from dials_scratch.jbe.sys_abs.plots import plot_screw_axes

@singleton
class ScrewAxisObserver(Observer):

    """Observer to record data used in screw axis analysis."""

    def update(self, screw_axis):
        self.data[screw_axis.name] = {
            'miller_axis_vals': screw_axis.miller_axis_vals,
            'i_over_sigma': screw_axis.i_over_sigma,
        }

    def generate_html_report(self, filename):
        """Generate a html report using the data."""
        screw_axes_graphs = plot_screw_axes(self.data)
        self.data['screw_axes'] = screw_axes_graphs
        loader = ChoiceLoader(
            [
                PackageLoader("dials", "templates"),
                PackageLoader("dials", "static", encoding="utf-8"),
            ]
        )
        env = Environment(loader=loader)
        template = env.get_template("systematic_absences_report.html")
        html = template.render(
            page_title="DIALS systematic absences report",
            screw_axes_graphs=self.data['screw_axes'],
        )
        with open(filename, "wb") as f:
            f.write(html.encode("ascii", "xmlcharrefreplace"))



class ScrewAxis(Subject):

    """Definition of a generic screw axis."""

    axis_idx = None #x=0, y=1, z=2
    axis_repeat = None # repeat of present reflections e.g =4 for 41, =2 for 42
    name = None
    orthogonal_vectors = (None, None) # two vectors orthogonal to the screw axis,
    # to allow determination of screw axis (these two not necessariy orthogonal
    # to each other).
    exclude_in_sum = [] # should any be excluded (e.g. don't test l=4 for 42 screw)

    def __init__(self):
        super(ScrewAxis, self).__init__(events=["selected data for scoring"])
        self.equivalent_axes = []
        self.n_refl_used = ()
        self.miller_axis_vals = []
        self.i_over_sigma = []

    def add_equivalent_axis(self, equivalent):
        """Add a symmetry equivalent axis."""
        self.equivalent_axes.append(equivalent)

    def select_axial_reflections(self, miller_indices):
        """Select reflections along the screw axis."""
        indices = miller_indices.as_vec3_double()

        v1 = flex.vec3_double(indices.size(), self.orthogonal_vectors[0])
        v2 = flex.vec3_double(indices.size(), self.orthogonal_vectors[1])

        a1 = v1.dot(indices)
        a2 = v2.dot(indices)
        selection = (a1 == 0.0) & (a2 == 0.0)
        return selection

    @Subject.notify_event(event="selected data for scoring")
    def get_all_suitable_reflections(self, reflection_table):
        """Select suitable reflections for testing the screw axis."""
        refl = reflection_table
        sel = self.select_axial_reflections(refl['miller_index'])
        miller_idx = refl['miller_index'].select(sel)
        self.i_over_sigma = refl['intensity'].select(sel) / (
            refl['variance'].select(sel) ** 0.5)
        self.miller_axis_vals = miller_idx.as_vec3_double().parts()[self.axis_idx]
        if self.equivalent_axes:
            for a in self.equivalent_axes:
                sel = a.select_axial_reflections(refl['miller_index'])
                miller_idx = refl['miller_index'].select(sel)
                extra_i_over_sigma = refl['intensity'].select(sel) / (
                    refl['variance'].select(sel) ** 0.5)
                extra_miller_axis_vals = miller_idx.as_vec3_double().parts()[a.axis_idx]
                self.miller_axis_vals.extend(extra_miller_axis_vals)
                self.i_over_sigma.extend(extra_i_over_sigma)

        if self.exclude_in_sum:
            for repeat_to_exclude in self.exclude_in_sum:
                #need to ignore those which might be due to screw with lower translation
                include = (self.miller_axis_vals.iround() % repeat_to_exclude != 0)
                self.miller_axis_vals = self.miller_axis_vals.select(include)
                self.i_over_sigma = self.i_over_sigma.select(include)


    def score_axis(self, reflection_table):
        """Score the axis give a reflection table of data."""
        self.get_all_suitable_reflections(reflection_table)

        expected_sel = (self.miller_axis_vals.iround() % self.axis_repeat == 0)

        expected = self.i_over_sigma.select(expected_sel)
        expected_abs = self.i_over_sigma.select(~expected_sel)

        # Limit to best #n reflections to avoid weak at high res?
        self.n_refl_used = (expected.size(), expected_abs.size())
        expected_sum = flex.sum(expected)
        expected_abs_sum = flex.sum(expected_abs)
        #need some uncertainty estimates to help when all values low!
        if expected_sum < 0.0:
            score = 0
        else:
            F0 = expected_abs_sum + expected_sum
            F_frac = expected_sum - expected_abs_sum #fractional F e.g F1/2 for 2_1, F1/3 for 3_1 etc
            score = max(min(F_frac/F0, 1.0), 0.0)
        #print("score for %s: %s" % (self.name, score))
        #print(expected_abs_sum, expected_sum)

        return score

class ScrewAxis21c(ScrewAxis):

    """Definition of a 21c screw axis"""

    axis_idx = 2
    axis_repeat = 2
    name = "21c"
    orthogonal_vectors = ((1, 0, 0), (0, 1, 0))


class ScrewAxis21b(ScrewAxis):

    """Definition of a 21b screw axis"""

    axis_idx = 1
    axis_repeat = 2
    name = "21b"
    orthogonal_vectors = ((1, 0, 0), (0, 0, 1))

class ScrewAxis21a(ScrewAxis):

    """Definition of a 21a screw axis"""

    axis_idx = 0
    axis_repeat = 2
    name = "21a"
    orthogonal_vectors = ((0, 1, 0), (0, 0, 1))

class ScrewAxis41c(ScrewAxis):

    """Definition of a 41c screw axis"""

    axis_idx = 2
    axis_repeat = 4
    name = "41c"
    orthogonal_vectors = ((0, 1, 0), (1, 0, 0))


class ScrewAxis42c(ScrewAxis):

    """Definition of a 42c screw axis"""

    axis_idx = 2
    axis_repeat = 2
    exclude_in_sum = [4]
    name = "42c"
    orthogonal_vectors = ((0, 1, 0), (1, 0, 0))


class ScrewAxis41b(ScrewAxis):

    """Definition of a 41b screw axis"""

    axis_idx = 1
    axis_repeat = 4
    name = "41c"
    orthogonal_vectors = ((0, 0, 1), (1, 0, 0))


class ScrewAxis41a(ScrewAxis):

    """Definition of a 41a screw axis"""

    axis_idx = 0
    axis_repeat = 4
    name = "41c"
    orthogonal_vectors = ((0, 1, 0), (0, 0, 1))


class ScrewAxis31c(ScrewAxis):

    """Definition of a 31c screw axis"""

    axis_idx = 2
    axis_repeat = 3
    name = "31c"
    orthogonal_vectors = ((0, 1, 0), (1, 0, 0))


class ScrewAxis61c(ScrewAxis):

    """Definition of a 61c screw axis"""

    axis_idx = 2
    axis_repeat = 6
    name = "61c"
    orthogonal_vectors = ((0, 1, 0), (1, 0, 0))


class ScrewAxis62c(ScrewAxis):

    """Definition of a 62c screw axis"""

    axis_idx = 2
    axis_repeat = 3
    name = "62c"
    orthogonal_vectors = ((0, 1, 0), (1, 0, 0))
    exclude_in_sum = [6]


class ScrewAxis63c(ScrewAxis):

    """Definition of a 63c screw axis"""

    axis_idx = 2
    axis_repeat = 2
    name = "63c"
    orthogonal_vectors = ((0, 1, 0), (1, 0, 0))
    exclude_in_sum = [6, 3]
