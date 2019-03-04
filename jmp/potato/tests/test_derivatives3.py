from __future__ import division
from __future__ import print_function
from scitbx import matrix
from math import log, pi
from random import uniform, randint
from dials_scratch.jmp.potato.model import compute_change_of_basis_operation
from dials_scratch.jmp.potato.parameterisation import SimpleMosaicityParameterisation
from dxtbx.model.experiment_list import ExperimentListFactory
from dials.algorithms.refinement.parameterisation.crystal_parameters import (
    CrystalUnitCellParameterisation,
)
from dials.algorithms.refinement.parameterisation.crystal_parameters import (
    CrystalOrientationParameterisation,
)
from dials_scratch.jmp.potato.profile_refiner import ReflectionData
from dials.array_family import flex
from dials_scratch.jmp.potato.parameterisation import ModelState
from dials_scratch.jmp.potato.parameterisation import ReflectionModelState


def first_derivative(func, x, h):
    return (-func(x + 2 * h) + 8 * func(x + h) - 8 * func(x - h) + func(x - 2 * h)) / (
        12 * h
    )


def second_derivative(func, x, y=None, h=None):
    if y is None:
        A = func(x + 2 * h)
        B = func(x + h)
        C = func(x)
        D = func(x - h)
        E = func(x - 2 * h)
        return (-(1 / 12) * (A + E) + (4 / 3) * (B + D) - (5 / 2) * C) / h ** 2
    else:
        A = func(x - h, y - h)
        B = func(x - h, y)
        C = func(x, y - h)
        D = func(x, y)
        E = func(x, y + h)
        F = func(x + h, y)
        G = func(x + h, y + h)
        return (A - B - C + 2 * D - E - F + G) / (2 * h ** 2)


def generate_data(experiments, reflections):

    from random import seed

    seed(0)

    index = randint(0, len(reflections))

    h = reflections[index]["miller_index"]

    s0 = matrix.col(experiments[0].beam.get_s0())

    U_param = CrystalOrientationParameterisation(experiments[0].crystal)
    B_param = CrystalUnitCellParameterisation(experiments[0].crystal)

    U = matrix.sqr(experiments[0].crystal.get_U())
    B = matrix.sqr(experiments[0].crystal.get_B())
    r = U * B * matrix.col(h)
    s2 = s0 + r
    mobs = (
        s2 + matrix.col((uniform(0, 1e-3), uniform(0, 1e-3), uniform(0, 1e-3)))
    ).normalize() * s0.length()

    b1, b2, b3, b4, b5, b6 = (
        uniform(1e-3, 3e-3),
        uniform(0.0, 1e-3),
        uniform(1e-3, 3e-3),
        uniform(0.0, 1e-3),
        uniform(0.0, 1e-3),
        uniform(1e-3, 3e-3),
    )

    params = (b1, b2, b3, b4, b5, b6)

    S_param = SimpleMosaicityParameterisation(params)
    L_param = (uniform(1e-3, 2e-3),)
    W_param = (uniform(1e-3, 2e-3), uniform(0, 1e-3), uniform(1e-3, 2e-3))
    ctot = randint(100, 1000)

    T = matrix.sqr((uniform(1e-3, 2e-3), 0, uniform(1e-6, 2e-6), uniform(1e-3, 2e-3)))
    Sobs = T * T.transpose()

    params = [S_param, U_param, B_param, L_param, W_param]

    return params, s0, h, ctot, mobs, Sobs


def test_first_derivatives(experiment, models, s0, h, ctot, mobs, Sobs):

    from dials_scratch.jmp.potato.refiner import ReflectionLikelihood

    U_params = models[1].get_param_vals()
    B_params = models[2].get_param_vals()
    M_params = flex.double(models[0].parameters())
    L_params = flex.double(models[3])
    W_params = flex.double(models[4])

    state = ModelState(experiment.crystal)
    state.set_U_params(U_params)
    state.set_B_params(B_params)
    state.set_M_params(M_params)
    state.set_L_params(L_params)
    state.set_W_params(W_params)

    def compute_L(parameters):
        state.set_active_parameters(parameters)
        model = ReflectionModelState(state, s0, h)

        rd = ReflectionLikelihood(model, s0, mobs, h, ctot, matrix.col((0, 0)), Sobs)
        return rd.log_likelihood()

    step = 1e-5

    parameters = state.get_active_parameters()

    dL_num = []
    for i in range(len(parameters)):

        def f(x):
            p = [pp for pp in parameters]
            p[i] = x
            return compute_L(p)

        dL_num.append(first_derivative(f, parameters[i], step))

    state.set_active_parameters(parameters)
    model = ReflectionModelState(state, s0, h)

    rd = ReflectionLikelihood(model, s0, mobs, h, ctot, matrix.col((0, 0)), Sobs)

    dL_cal = list(rd.first_derivatives())

    assert all(abs(n - c) < 1e-3 for n, c in zip(dL_num, dL_cal))

    print("OK")


def test_reflection_model(experiment, models, s0, h):

    U_params = models[1].get_param_vals()
    B_params = models[2].get_param_vals()
    M_params = flex.double(models[0].parameters())
    L_params = flex.double(models[3])
    W_params = flex.double(models[4])

    state = ModelState(experiment.crystal)
    state.set_U_params(U_params)
    state.set_B_params(B_params)
    state.set_M_params(M_params)
    state.set_L_params(L_params)
    state.set_W_params(W_params)

    model = ReflectionModelState(state, s0, h)

    r = model.get_r()
    sigma = model.get_sigma()
    dr_dp = model.get_dr_dp()
    dS_dp = model.get_dS_dp()

    def compute_sigma(parameters):
        state.set_active_parameters(parameters)
        model = ReflectionModelState(state, s0, h)
        return model.get_sigma()

    def compute_r(parameters):
        state.set_active_parameters(parameters)
        model = ReflectionModelState(state, s0, h)
        return model.get_r()

    step = 1e-6

    parameters = state.get_active_parameters()

    dr_num = []
    for i in range(len(parameters)):

        def f(x):
            p = [pp for pp in parameters]
            p[i] = x
            return compute_r(p)

        dr_num.append(first_derivative(f, parameters[i], step))

    for n, c in zip(dr_num, dr_dp):
        assert all(abs(nn - cc) < 1e-7 for nn, cc in zip(n, c))

    ds_num = []
    for i in range(len(parameters)):

        def f(x):
            p = [pp for pp in parameters]
            p[i] = x
            return compute_sigma(p)

        ds_num.append(first_derivative(f, parameters[i], step))

    for n, c in zip(ds_num, dS_dp):
        assert all(abs(nn - cc) < 1e-7 for nn, cc in zip(n, c))

    print("OK")


def read_experiments():
    experiments = ExperimentListFactory.from_json_file("experiments.json")
    return experiments


def test():

    experiments = read_experiments()
    experiments[0].scan.set_oscillation((0, 0.01), deg=True)
    reflections = flex.reflection_table.from_predictions_multi(experiments)

    for i in range(1):

        models, s0, h, ctot, mobs, Sobs = generate_data(experiments, reflections)
        test_reflection_model(experiments[0], models, s0, h)
        test_first_derivatives(experiments[0], models, s0, h, ctot, mobs, Sobs)


if __name__ == "__main__":
    test()
