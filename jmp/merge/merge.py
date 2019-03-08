from __future__ import division
from libtbx import phil
from dxtbx.model.experiment_list import ExperimentList
from scipy.optimize import least_squares
from collections import Counter, defaultdict
from dials.array_family import flex
from math import floor, sqrt, ceil, log, pi
import numpy
import cPickle as pickle
import logging

logger = logging.getLogger("dials")

phil_scope = phil.parse(
    """

  tolerance = 1
    .type = float
    .help = "Tolerance to test for convergence (difference in log(L))"
  
  max_iter = 10000
    .type = int
    .help = "Maximum number of iterations"

  scales = None
    .type = str
    .help = "The scales pickle file"

  scaling {
    d_min = 4
      .type = float
      .help = "Minimum resolution for scaling"
  }

"""
)


class Scaler(object):
    def __init__(self, tolerance=1e-3, max_iter=10000):

        # Set the parameters
        self.tolerance = tolerance
        self.max_iter = max_iter

        # Initialise the reflections
        self.reflections = None

        # Initialise the merged intensities
        self.h_unique = None
        self.I_mean = None
        self.V_mean = None

        # Initialise the scales
        self.image_scale = None
        self.image_bfactor = None
        self.image_mosaicity = None
        self.image_pfactor = None

    def scale(self, experiments, reflections, reference=None):

        # Match reflections with reference
        if reference is not None:
            reflections, reference = self._match_with_reference(
                reflections, experiments, reference
            )

        # Preprocess the data
        experiments, reflections = self._preprocess(experiments, reflections, reference)

        # Save the reflections
        self.reflections = reflections

        # Initialise the starting values
        # self._initialise(self.reflections, reference)

        # Solve equations
        self._solve(
            self.reflections["intensity.sum.value"],
            self.reflections["intensity.sum.variance"],
            self.reflections["theta"],
            self.reflections["epsilon"],
            self.reflections["correction"],
            self.reflections["wavelength"],
            self.reflections["image"],
            self.reflections["reflection_id"],
        )

    def _match_with_reference(self, reflections, experiments, reference):

        # Compute the miller indices in the asu
        reflections.compute_miller_indices_in_asu(experiments)

        # Item to hold matches
        class Item(object):
            def __init__(self):
                self.a = None
                self.b = flex.size_t()

        # Search for matches
        lookup = defaultdict(Item)
        for i, h in enumerate(reference["miller_index"]):
            lookup[h].a = i
        for i, h in enumerate(reflections["miller_index_asu"]):
            lookup[h].b.append(i)

        # Compute the selections
        selection1 = flex.size_t()
        selection2 = flex.size_t()
        for key, value in lookup.iteritems():
            if value.a is not None and len(value.b) > 0:
                selection1.append(value.a)
                selection2.extend(value.b)

        # Select the reflections
        reference = reference.select(selection1)
        reflections = reflections.select(selection2)

        # Return
        return reflections, reference

    def _preprocess(self, experiments, reflections, reference):

        # Filter reflections by:
        # - integrated
        # - non-zero variance
        # - valid id
        selection1 = reflections.get_flags(reflections.flags.integrated_sum)
        selection2 = reflections["intensity.sum.variance"] > 0
        selection3 = reflections["id"] >= 0
        selection = selection1 & selection2 & selection3
        reflections = reflections.select(selection)

        # Ensure we aren't missing experiments
        num_images = len(set(reflections["id"]))
        if num_images != len(experiments):
            new_experiments = ExperimentList()
            new_reflections = flex.reflection_table()
            counter = 0
            for experiment, indices in reflections.iterate_experiments_and_indices(
                experiments
            ):
                if len(indices) > 0:
                    subset = reflections.select(indices)
                    subset["id"] = flex.size_t(len(subset), counter)
                    new_experiments.append(experiment)
                    new_reflections.extend(subset)
                    counter += 1
            experiments = new_experiments
            reflections = new_reflections

        # Convert id to image
        reflections["image"] = flex.size_t(list(reflections["id"]))

        def compute_epsilon(reflections):
            if "s2" not in reflections:
                A = flex.mat3_double([e.crystal.get_A() for e in experiments])
                s0 = flex.vec3_double([e.beam.get_s0() for e in experiments])
                A = A.select(reflections["image"])
                s0 = s0.select(reflections["image"])
                h = flex.vec3_double(list(reflections["miller_index"]))
                r = A * h
                s2 = s0 + r
                s1 = s0.norms() * s2 / s2.norms()
                reflections["s1"] = s1
                reflections["s2"] = s2
            else:
                s2 = reflections["s2"]
                s1 = reflections["s1"]
            epsilon = s2.norms() - s1.norms()
            reflections["epsilon"] = epsilon

        def compute_theta(reflections):
            if "theta" not in reflections:
                panels = [
                    experiments[i].detector[j]
                    for i, j in zip(reflections["image"], reflections["panel"])
                ]
                s0 = flex.vec3_double([e.beam.get_s0() for e in experiments])
                s0 = s0.select(reflections["image"])
                x, y, _ = reflections["xyzcal.px"].parts()
                xy = flex.vec2_double(zip(x, y))
                two_theta = flex.double(
                    panels[i].get_two_theta_at_pixel(s0[i], xy[i])
                    for i in range(len(s0))
                )
                reflections["theta"] = two_theta / 2.0
                reflections["wavelength"] = 1.0 / s0.norms()

        # Compute the resolution
        reflections.compute_d(experiments)

        # Compute distance to Ewald sphere
        compute_epsilon(reflections)

        # Compute the theta angle
        compute_theta(reflections)

        # Compute the miller indices in the asu
        reflections.compute_miller_indices_in_asu(experiments)

        # Compute the unique miller indices
        h_unique = flex.miller_index(sorted(list(set(reflections["miller_index_asu"]))))
        h_lookup = dict((h, i) for i, h in enumerate(h_unique))
        h_index = flex.size_t()
        for h in reflections["miller_index_asu"]:
            h_index.append(h_lookup[h])

        # Save the reflection data
        reflections["reflection_id"] = h_index

        # Create a lookup for reflection index
        lookup_index = [flex.size_t() for i in range(len(h_unique))]
        for i, h in enumerate(reflections["reflection_id"]):
            lookup_index[h].append(i)

        # Create a lookup for image number
        num_images = len(set(reflections["id"]))
        assert num_images == len(experiments)
        assert num_images == flex.max(reflections["id"]) + 1
        lookup_image = [flex.size_t() for i in range(num_images)]
        for i, j in enumerate(reflections["id"]):
            lookup_image[j].append(i)

        # Create reference lookup
        reference_lookup = flex.size_t()
        for h0 in h_unique:
            for i, h in enumerate(reference["miller_index"]):
                if h == h0:
                    reference_lookup.append(i)
                    break

        # Get reference intensities
        I_ref = flex.double()
        for h in range(len(lookup_index)):
            I_ref.append(reference[reference_lookup[h]]["intensity.ref.value"])

        # Save the reference intensity
        self.I_mean = I_ref

        # Save the arrays
        self.h_unique = h_unique
        self.lookup_index = lookup_index
        self.lookup_image = lookup_image

        # Compute multiplicity and # obs on image
        multiplicity = flex.int(len(x) for x in lookup_index)
        obs_on_image = flex.int(len(x) for x in lookup_image)
        assert flex.sum(multiplicity) == len(reflections)
        assert flex.sum(obs_on_image) == len(reflections)

        # Print some info
        logger.info("# Total reflections:  %d" % len(reflections))
        logger.info("# Unique reflections: %d" % len(self.h_unique))
        logger.info("# Images: %d" % len(self.lookup_image))
        logger.info("Min multiplicity: %d" % flex.min(multiplicity))
        logger.info("Max multiplicity: %d" % flex.max(multiplicity))
        logger.info("Min observations on image: %d" % flex.min(obs_on_image))
        logger.info("Max observations on image: %d" % flex.max(obs_on_image))

        # Return reflections
        return experiments, reflections

    def _initialise(self, reflections, reference=None):

        # Get the intensity data
        I = reflections["intensity.sum.value"]

        # Get number of classes
        num_images = len(self.lookup_image)

        # Initialise the scale factors
        g_old = flex.double(num_images, 1.0)

        # Initialise the intensities
        c_old = flex.double()
        for selection in self.lookup_index:
            I_i = I.select(selection)
            c_old.append(flex.sum(I_i) / len(I_i))

        # Initialise the partiality parameters
        v_old = 1.0
        m_old = sqrt(1e-7)

        # Set the initial values
        self.image_scale = g_old
        self.partiality_v = v_old
        self.partiality_m = m_old
        # self.I_mean = c_old

    def _solve_for_image(
        self,
        I_obs,
        V_obs,
        correction,
        epsilon,
        theta,
        wavelength,
        J_ref,
        I_ref,
        V0,
        S0,
        B0,
        G0,
    ):
        from scipy.optimize import least_squares
        import numpy as np

        def function(p):
            """
      The runction of residuals

      """
            # Get the parameters
            V, S, B, G = flex.double(p)

            # Get the reference intensity
            I_i = I_ref.select(J_ref)

            # Compute the B factor correction
            b = flex.exp(-2 * B * (flex.sin(theta) / wavelength) ** 2)

            # Compute the expected partiality for the observation
            p = (1 + (1 / V) * (epsilon / S) ** 2) ** (-(V + 1) / 2)

            # from matplotlib import pylab
            # pylab.scatter(epsilon, I_obs/I_i)
            # pylab.show()

            # Compute the expected value of the observation
            mu = correction * G * b * p * I_i

            # Compute the residuals
            residuals = (I_obs - mu) / flex.sqrt(V_obs)

            # Print out some information
            cost = flex.sum(residuals ** 2)
            logger.info(
                "%d: cost=%.5e V=%.5f S=%.5f B=%.3f G=%.3f"
                % (self.it, cost, V, S, B, G)
            )

            # Update the iteration counter
            self.it += 1

            # Return the residuals
            return residuals

        def jacobian(p):
            """
      Compute the Jacobian matrix of residuals

      """

            # Get the parameters
            V, S, B, G = flex.double(p)

            # Get the reference intensity
            I_i = I_ref.select(J_ref)

            # Compute the B factor correction
            b = flex.exp(-2 * B * (flex.sin(theta) / wavelength) ** 2)

            # Compute the expected partiality for the observation
            p = (1 + (1 / V) * (epsilon / S) ** 2) ** (-(V + 1) / 2)

            # Compute the expected value of the observation
            mu = G * b * p * I_i

            # Compute the derivative of p(V,S) wrt V
            num = (
                V
                * (S ** 2 * V + epsilon ** 2)
                * flex.log(1 + (1 / V) * (epsilon / S) ** 2)
                - (V + 1) * epsilon ** 2
            )
            den = 2 * V * (S ** 2 * V + epsilon ** 2)
            dp_dV = -p * num / den

            # Compute the derivative of p(V,S) wrt S
            num = (
                (1 + (1 / V) * (epsilon / S) ** 2) ** (-(V + 1) / 2 - 1)
                * (V + 1)
                * epsilon ** 2
            )
            den = S ** 3 * V
            dp_dS = num / den

            # The derivatives of each expected value wrt to parameters that effect them
            dmu_dV = correction * G * b * I_i * dp_dV / flex.sqrt(V_obs)
            dmu_dS = correction * G * b * I_i * dp_dS / flex.sqrt(V_obs)
            dmu_dB = (
                correction
                * mu
                * (-2 * (flex.sin(theta) / wavelength) ** 2)
                / flex.sqrt(V_obs)
            )
            dmu_dG = correction * b * p * I_i / flex.sqrt(V_obs)

            # Construct the Jacobian
            n_par = 4
            J = flex.double(flex.grid(len(I_obs), n_par), 0)
            for j in range(len(I_obs)):
                J[j, 0] = 0  # -dmu_dV[j]
                J[j, 1] = 0  # -dmu_dS[j]
                J[j, 2] = 0  # -dmu_dB[j]
                J[j, 3] = -dmu_dG[j]

            # Return the Jacobian
            return J.as_numpy_array()

        # Get number of images and reflections etc
        n_ref = len(I_ref)

        # Set the iteration counter
        self.it = 0

        # The initial parameters
        p0 = flex.double([V0, S0, B0, G0])

        # The parameter bounds
        bounds = ([1e-7, 1e-20, -np.inf, -np.inf], [np.inf, np.inf, 0, np.inf])

        # Do the least squares minimization
        result = least_squares(function, p0, jac=jacobian, bounds=bounds)

        # Get the parameter estimates
        V, S, B, G = result.x
        return V, S, B, G

    def _solve_for_all(
        self, I_obs, V_obs, correction, epsilon, theta, J_img, J_ref, I_ref, wavelength
    ):

        n_ref = len(I_ref)
        n_img = flex.max(J_img) + 1

        # The initial parameters
        V = flex.double(1.0 for i in range(n_img))
        S = flex.double(1.0 for i in range(n_img))
        B = flex.double(0.0 for i in range(n_img))
        G = flex.double(1.0 for i in range(n_img))

        # Vt = 1.0
        # St = 1e-7
        # Gt = 1.0
        for i in range(n_img):

            selection = J_img == i
            logger.info("Image: %d, # Reflections: %d" % (i, selection.count(True)))

            # V[i] = Vt
            # S[i] = St
            # G[i] = Gt

            V_i, S_i, B_i, G_i = self._solve_for_image(
                I_obs.select(selection),
                V_obs.select(selection),
                correction.select(selection),
                epsilon.select(selection),
                theta.select(selection),
                wavelength.select(selection),
                J_ref.select(selection),
                I_ref,
                V[i],
                S[i],
                B[i],
                G[i],
            )

            V[i] = V_i
            S[i] = S_i
            B[i] = B_i
            G[i] = G_i

            Vt = flex.sum(V[: i + 1]) / (i + 1)
            St = flex.sum(S[: i + 1]) / (i + 1)
            Gt = flex.sum(G[: i + 1]) / (i + 1)

        # Return the result
        return V, S, B, G

    def _compute_partiality_and_intensity(
        self,
        I_obs,
        V_obs,
        correction,
        epsilon,
        theta,
        wavelength,
        J_img,
        J_ref,
        V,
        S,
        B,
        G,
    ):
        """
    Compute the partiality and intensity from the parameters

    """

        # Get the sigma, b factor and scale for each reflection
        V_i = V.select(J_img)
        S_i = S.select(J_img)
        B_i = B.select(J_img)
        G_i = G.select(J_img)

        # Compute the B factor correction
        b_i = flex.exp(-2 * B_i * (flex.sin(theta) / wavelength) ** 2)

        # Compute the expected partiality for the observation
        p_i = flex.double(
            (1 + (1 / V_i[i]) * (epsilon[i] / S_i[i]) ** 2) ** (-(V_i[i] + 1) / 2)
            for i in range(len(V_i))
        )

        # Compute the intensities#
        n_ref = flex.max(J_ref) + 1
        I_est = flex.double(n_ref)
        for i in range(n_ref):
            selection = J_ref == i
            g = G_i.select(selection)
            b = b_i.select(selection)
            p = p_i.select(selection)
            I = I_obs.select(selection)
            var_I = V_obs.select(selection)
            c = correction.select(selection)

            num = flex.sum(I * c * g * b * p / var_I)
            den = flex.sum((c * g * b * p) ** 2 / var_I)

            I_est[i] = num / den

        # Return the intensity
        return p_i, I_est

    def _solve(
        self, I_obs, V_obs, theta, epsilon, correction, wavelength, img_index, ref_index
    ):

        # Solve for all images
        V, S, B, G = self._solve_for_all(
            I_obs,
            V_obs,
            correction,
            epsilon,
            theta,
            img_index,
            ref_index,
            self.I_mean,
            wavelength,
        )
        # with open("scales.pickle") as infile:
        #   result = pickle.load(infile)
        #   V = result['image_pfactor']
        #   S = result['image_mosaicity']
        #   B = result['image_bfactor']
        #   G = result['image_scale']

        # Compute the intensity
        P, I_est = self._compute_partiality_and_intensity(
            I_obs,
            V_obs,
            correction,
            epsilon,
            theta,
            wavelength,
            img_index,
            ref_index,
            V,
            S,
            B,
            G,
        )

        from matplotlib import pylab

        xmin = min(self.I_mean)
        xmax = max(self.I_mean)
        pylab.scatter(self.I_mean, I_est, s=1, alpha=0.5)
        pylab.plot([xmin, xmax], [xmin, xmax], color="black")
        pylab.show()

        # Set the results
        self.I_mean = I_est
        self.image_scale = G
        self.image_bfactor = B
        self.image_mosaicity = S
        self.image_pfactor = V

        # Set some reflection data
        self.reflections["partiality"] = P
        self.reflections["image_scale"] = G.select(img_index)
        self.reflections["image_bfactor"] = B.select(img_index)
        self.reflections["image_mosaicity"] = S.select(img_index)
        self.reflections["image_pfactor"] = V.select(img_index)
        self.reflections["Imean"] = I_est.select(ref_index)


def scale_and_merge(experiments, reflections, params, reference):

    # reflections = reflections[0:1000]

    # Compute corrections (here LP is just 1/P)
    # I_true * correction = I
    reflections.compute_corrections(experiments)
    reflections["correction"] = reflections["qe"] / reflections["lp"]

    reflections = reflections.select(reflections["d"] > params.scaling.d_min)

    # Configure the scaler
    scaler = Scaler(tolerance=params.tolerance, max_iter=params.max_iter)

    # Run the scaler
    scaler.scale(experiments, reflections, reference)

    # Get the reflections
    reflections = scaler.reflections

    # Get the results
    h_unique = scaler.h_unique
    I_mean = scaler.I_mean
    V_mean = scaler.V_mean

    result = {
        "h_unique": h_unique,
        "I_mean": I_mean,
        "image_scale": scaler.image_scale,
        "image_bfactor": scaler.image_bfactor,
        "image_mosaicity": scaler.image_mosaicity,
        "image_pfactor": scaler.image_pfactor,
    }

    reflections.as_pickle("output_reflections.pickle")

    pickle.dump(result, open("scales.pickle", "wb"))
