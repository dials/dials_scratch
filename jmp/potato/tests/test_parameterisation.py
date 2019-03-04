from __future__ import division
from __future__ import print_function
from dials_scratch.jmp.potato.parameterisation import ModelState
from dxtbx.model.experiment_list import ExperimentListFactory
from dials.array_family import flex


def test_model_state_with_fixed(
    crystal,
    fix_rlp_mosaicity=False,
    fix_wavelength_spread=False,
    fix_angular_mosaicity=False,
    fix_unit_cell=False,
    fix_orientation=False,
):

    state = ModelState(
        crystal,
        fix_rlp_mosaicity=fix_rlp_mosaicity,
        fix_wavelength_spread=fix_wavelength_spread,
        fix_angular_mosaicity=fix_angular_mosaicity,
        fix_unit_cell=fix_unit_cell,
        fix_orientation=fix_orientation,
    )

    assert state.is_orientation_fixed() == fix_orientation
    assert state.is_unit_cell_fixed() == fix_unit_cell
    assert state.is_rlp_mosaicity_fixed() == fix_rlp_mosaicity
    assert state.is_wavelength_spread_fixed() == fix_wavelength_spread
    assert state.is_angular_mosaicity_fixed() == fix_angular_mosaicity

    U = state.get_U()
    B = state.get_B()
    A = state.get_A()
    M = state.get_M()
    L = state.get_L()
    W = state.get_W()

    assert len(U) == 9
    assert len(B) == 9
    assert len(A) == 9
    assert len(M) == 9
    assert len(L) == 9
    assert len(W) == 9

    U_params = state.get_U_params()
    B_params = state.get_B_params()
    M_params = state.get_M_params()
    L_params = state.get_L_params()
    W_params = state.get_W_params()

    assert len(U_params) == 3
    assert len(B_params) == 2
    assert len(M_params) == 6
    assert len(L_params) == 1
    assert len(W_params) == 3

    dU = state.get_dU_dp()
    dB = state.get_dB_dp()
    dM = state.get_dM_dp()
    dL = state.get_dL_dp()
    dW = state.get_dW_dp()

    assert len(dU) == 3
    assert len(dB) == 2
    assert len(dM) == 6
    assert len(dL) == 1
    assert len(dW) == 3

    params = state.get_active_parameters()

    expected_len = 0
    if not fix_rlp_mosaicity:
        expected_len += 6
    if not fix_wavelength_spread:
        expected_len += 1
    if not fix_angular_mosaicity:
        expected_len += 3
    if not fix_unit_cell:
        expected_len += 2
    if not fix_orientation:
        expected_len += 3

    assert len(params) == expected_len
    new_params = params
    state.set_active_parameters(new_params)

    print("OK")


def test_model_state(experiments):

    crystal = experiments[0].crystal

    test_model_state_with_fixed(crystal, fix_rlp_mosaicity=True)
    test_model_state_with_fixed(crystal, fix_wavelength_spread=True)
    test_model_state_with_fixed(crystal, fix_angular_mosaicity=True)
    test_model_state_with_fixed(crystal, fix_unit_cell=True)
    test_model_state_with_fixed(crystal, fix_orientation=True)


def read_experiments():
    experiments = ExperimentListFactory.from_json_file("experiments.json")
    return experiments


def test():

    experiments = read_experiments()
    experiments[0].scan.set_oscillation((0, 0.01), deg=True)

    test_model_state(experiments)


if __name__ == "__main__":
    test()
