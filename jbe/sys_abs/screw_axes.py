"""Definitions of screw axes with methods for scoring against data."""
import math
import logging
from scitbx.array_family import flex
from dials.util.observer import Observer, Subject, singleton
from jinja2 import Environment, ChoiceLoader, PackageLoader
from dials_scratch.jbe.sys_abs.plots import plot_screw_axes

logger = logging.getLogger('dials.absences')

@singleton
class ScrewAxisObserver(Observer):

    """Observer to record data used in screw axis analysis."""

    def update(self, screw_axis):
        self.data[screw_axis.name] = {
            'miller_axis_vals': screw_axis.miller_axis_vals,
            'i_over_sigma': screw_axis.i_over_sigma,
            'intensities': screw_axis.intensities,
            'sigmas': screw_axis.sigmas,
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


def fourfold_axis_test():
    # first test 4_1, if not 4_1 then try 4_2
    pass

class ScrewAxis(Subject):

    """Definition of a generic screw axis."""

    axis_idx = None #x=0, y=1, z=2
    axis_repeat = None # repeat of present reflections e.g =4 for 41, =2 for 42
    name = None
    orthogonal_vectors = (None, None) # two vectors orthogonal to the screw axis,
    # to allow determination of screw axis (these two not necessariy orthogonal
    # to each other).

    def __init__(self):
        super(ScrewAxis, self).__init__(events=["selected data for scoring"])
        self.equivalent_axes = []
        self.n_refl_used = (0.0, 0.0)
        self.miller_axis_vals = []
        self.i_over_sigma = []
        self.intensities = []
        self.sigmas = []
        self.mean_I_sigma_abs = 0.0
        self.mean_I_sigma = 0.0
        self.mean_I_abs = 0.0
        self.mean_I = 0.0

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
        self.intensities = refl['intensity'].select(sel)
        self.sigmas = refl['variance'].select(sel) ** 0.5

        if self.equivalent_axes:
            for a in self.equivalent_axes:
                sel = a.select_axial_reflections(refl['miller_index'])
                miller_idx = refl['miller_index'].select(sel)
                extra_i_over_sigma = refl['intensity'].select(sel) / (
                    refl['variance'].select(sel) ** 0.5)
                extra_i = refl['intensity'].select(sel)
                extra_sig = refl['variance'].select(sel) ** 0.5
                extra_miller_axis_vals = miller_idx.as_vec3_double().parts()[a.axis_idx]
                self.miller_axis_vals.extend(extra_miller_axis_vals)
                self.i_over_sigma.extend(extra_i_over_sigma)
                self.intensities.extend(extra_i)
                self.sigmas.extend(extra_sig)

    def score_axis(self, reflection_table, significance_level=0.95):
        """Score the axis give a reflection table of data."""
        assert significance_level in [0.95, 0.975, 0.99]
        self.get_all_suitable_reflections(reflection_table)

        expected_sel = (self.miller_axis_vals.iround() % self.axis_repeat == 0)

        expected = self.i_over_sigma.select(expected_sel)
        expected_abs = self.i_over_sigma.select(~expected_sel)
        self.n_refl_used = (expected.size(), expected_abs.size())

        if not expected or not expected_abs:
            return 0.0

        # Limit to best #n reflections to avoid weak at high res - use wilson B?
        self.n_refl_used = (expected.size(), expected_abs.size())

        # z = (sample mean - population mean) / standard error
        S_E_abs = 1.0 # errors probably correlated so say standard error = 1
        S_E_pres = 1.0 # / expected.size() ** 0.5

        self.mean_I_sigma_abs = flex.mean(expected_abs)
        self.mean_I_sigma = flex.mean(expected)

        self.mean_I = flex.mean(self.intensities.select(expected_sel))
        self.mean_I_abs = flex.mean(self.intensities.select(~expected_sel))

        z_score_absent = self.mean_I_sigma_abs / S_E_abs
        z_score_present = self.mean_I_sigma / S_E_pres

        # get a p-value for z > z_score
        P_absent = 0.5 * (1.0 + math.erf(z_score_absent / (2**0.5)))
        P_present = 0.5 * (1.0 + math.erf(z_score_present / (2**0.5)))

        # sanity check - is most of intensity in 'expected' channel?
        intensity_test = self.mean_I_sigma > (20.0 * self.mean_I_sigma_abs)

        cutoffs = {0.95 : 1.645, 0.975: 1.960, 0.99 : 2.326}
        cutoff = cutoffs[significance_level]

        if z_score_absent > cutoff and not intensity_test: # z > 1.65 in only 5% of cases for normal dist
            # significant nonzero intensity where expected absent.
            return (1.0 - P_absent) * P_present
        elif z_score_absent > cutoff:
            # results appear inconsistent - significant i_over_sigma_abs, but this
            # is still low compared to i_over_sigma_expected
            # try removing the highest absent reflection in case its an outlier
            logger.info(
"""Tests for %s appear inconsistent (significant nonzero intensity for 'absent'
reflections, but majority of intensity in reflection condition). Performing
outlier check.""", self.name)
            sel = flex.sort_permutation(expected_abs)
            sorted_exp_abs = expected_abs.select(sel)
            mean_i_sigma_abs = flex.mean(sorted_exp_abs[:-1])
            if (mean_i_sigma_abs / S_E_abs) > cutoff:
                logger.info("""Removed a reflection but still significant nonzero intensity
for 'absent' reflections.""")
                # Still high intensity of absent, so return as before
                return (1.0 - P_absent) * P_present
            logger.info("""Removed a reflection, intensities of 'absent' reflections no longer
significantly different to zero.""")
            self.mean_I_sigma_abs = mean_i_sigma_abs
            self.mean_I_abs = flex.mean(self.intensities.select(~expected_sel).select(sel)[:-1])
            # else was maybe misled by outlier, can continue to next test
        if z_score_present > cutoff: # evidence with confidence
            return P_present
        else:
            logger.info("""No evidence to suggest a screw axis for %s, but insufficient
evidence to rule out completely, possibly due to limited data.""", self.name)
            return 0.0


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
    name = "42c"
    orthogonal_vectors = ((0, 1, 0), (1, 0, 0))


class ScrewAxis41b(ScrewAxis):

    """Definition of a 41b screw axis"""

    axis_idx = 1
    axis_repeat = 4
    name = "41b"
    orthogonal_vectors = ((0, 0, 1), (1, 0, 0))


class ScrewAxis41a(ScrewAxis):

    """Definition of a 41a screw axis"""

    axis_idx = 0
    axis_repeat = 4
    name = "41a"
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


class ScrewAxis63c(ScrewAxis):

    """Definition of a 63c screw axis"""

    axis_idx = 2
    axis_repeat = 2
    name = "63c"
    orthogonal_vectors = ((0, 1, 0), (1, 0, 0))
