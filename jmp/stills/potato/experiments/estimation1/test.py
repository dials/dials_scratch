from __future__ import print_function
from numpy.random import normal, poisson, seed, uniform, multivariate_normal
from math import exp, pi, sqrt, log, sin, cos
from dials.array_family import flex
from scitbx import matrix
from dials_scratch.jmp.stills.potato.derivatives import (
    compute_change_of_basis_operation,
)
from dials_scratch.jmp.stills.potato.derivatives import MosaicityParameterisation
from dials_scratch.jmp.stills.potato.derivatives import ProfileModel
from dials_scratch.jmp.stills.potato.derivatives import ReflectionProfileModel
from dials_scratch.jmp.stills.potato.derivatives import estimate_parameters
from dials_scratch.jmp.stills.potato.experiments.generate_simple import generate_simple


def test_non_axis_aligned_random_orientation_likelihood():

    seed(0)

    M = matrix.sqr((1, 0, 0, 0.2, sqrt(2), 0, 0.3, 0.4, sqrt(3)))
    # M = matrix.sqr((
    #   1,   0,       0,
    #   0, 1, 0,
    #   0, 0,  1))

    sigma = M * M.transpose()

    params = M

    print(sigma)

    scale = 1e-4

    s0 = 10000 * matrix.col((0, 0, 1))

    s0 = s0 * scale
    sigma *= scale

    print("Generating")
    s2_list, ctot_list, Sobs_list = generate_simple(s0, sigma)

    def compute_simple(s0, s2_list):
        sigma = matrix.sqr((0, 0, 0, 0, 0, 0, 0, 0, 0))
        for s2 in s2_list:
            sigma += (s2 - s2.normalize() * s0.length()) * (
                s2 - s2.normalize() * s0.length()
            ).transpose()
        sigma /= len(s2_list)

        print("SIMPLE")
        print(sigma)
        return sigma

    sigma_sim = compute_simple(s0, s2_list)
    params = estimate_parameters(s0, s2_list, ctot_list, Sobs_list)

    M = matrix.sqr(
        (params[0], 0, 0, params[1], params[2], 0, params[3], params[4], params[5])
    )
    sigma_cal = M * M.transpose()

    # plot_parameters(s0, s2_list, ctot_list, Sobs_list, params)

    print(sigma)

    def kl_divergence(A, B):
        return 0.5 * (
            (B.inverse() * A).trace() - 3 + log(B.determinant() / A.determinant())
        )

    A = sigma
    B = sigma_cal
    C = sigma_sim
    print(kl_divergence(A, B))
    print(kl_divergence(A, C))
    # A = (sigma_cal - sigma)
    # print sqrt((A.transpose()*A).trace())


test_non_axis_aligned_random_orientation_likelihood()
