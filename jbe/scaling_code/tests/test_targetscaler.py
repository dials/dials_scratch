"""
Test for targetscaler and fixed Ih target function.
"""
import pytest
from dials.array_family import flex
from dials.util.options import OptionParser
from libtbx import phil
from libtbx.test_utils import approx_equal
from dxtbx.model.experiment_list import ExperimentList
from dxtbx.model import Crystal, Scan, Beam, Goniometer, Detector, Experiment
from dials_scratch.jbe.scaling_code.model.scaling_model_factory import \
  create_scaling_model
from dials_scratch.jbe.scaling_code.scaler_factory import TargetScaler
from dials_scratch.jbe.scaling_code.scaler import SingleScalerBase
from dials_scratch.jbe.scaling_code.target_function import ScalingTargetFixedIH
from dials_scratch.jbe.scaling_code.parameter_handler import create_apm
from test_basis_and_target_function import calculate_gradient_fd

def generated_refl():
  '''function to generate input for datamanagers'''
  #these miller_idx/d_values don't make physical sense, but I didn't want to
  #have to write the tests for lots of reflections.
  reflections = flex.reflection_table()
  reflections['intensity.prf.value'] = flex.double([1.0, 10.0, 1.0])
  reflections['intensity.prf.variance'] = flex.double([1.0, 10.0, 1.0])
  reflections['miller_index'] = flex.miller_index([(1, 0, 0), (0, 0, 1),
    (1, 0, 0)]) #don't change
  reflections['d'] = flex.double([0.8, 2.0, 2.0]) #don't change
  reflections['lp'] = flex.double([1.0, 1.0, 1.0])
  reflections['dqe'] = flex.double([1.0, 1.0, 1.0])
  reflections['partiality'] = flex.double([1.0, 1.0, 1.0])
  reflections['xyzobs.px.value'] = flex.vec3_double([(0.0, 0.0, 0.0),
    (0.0, 0.0, 5.0), (0.0, 0.0, 10.0)])
  reflections['s1'] = flex.vec3_double([(0.0, 0.1, 1.0), (0.0, 0.1, 1.0),
    (0.0, 0.1, 1.0)])
  reflections.set_flags(flex.bool([True, True, True]),
    reflections.flags.integrated)
  return [reflections]

def generated_refl_2():
  """Generate a second reflection table."""
  #these miller_idx/d_values don't make physical sense, but I didn't want to
  #have to write the tests for lots of reflections.
  reflections = flex.reflection_table()
  reflections['intensity.prf.value'] = flex.double([2.0, 10.0, 100.0])
  reflections['intensity.prf.variance'] = flex.double([1.0, 10.0, 100.0])
  reflections['miller_index'] = flex.miller_index([(1, 0, 0), (0, 0, 5),
    (5, 0, 0)]) #don't change
  reflections['d'] = flex.double([0.8, 2.0, 2.0]) #don't change
  reflections['lp'] = flex.double([1.0, 1.0, 1.0])
  reflections['dqe'] = flex.double([1.0, 1.0, 1.0])
  reflections['partiality'] = flex.double([1.0, 1.0, 1.0])
  reflections['xyzobs.px.value'] = flex.vec3_double([(0.0, 0.0, 0.0),
    (0.0, 0.0, 5.0), (0.0, 0.0, 10.0)])
  reflections['s1'] = flex.vec3_double([(0.0, 0.1, 1.0), (0.0, 0.1, 1.0),
    (0.0, 0.1, 1.0)])
  reflections['inverse_scale_factor'] = flex.double([1.0, 1.0, 1.0])
  reflections.set_flags(flex.bool([True, True, True]),
    reflections.flags.integrated)
  return [reflections]

#@pytest.fixture(scope='module')
def generated_single_exp():
  """Generate and experiments object."""
  experiments = ExperimentList()
  exp_dict = {"__id__" : "crystal", "real_space_a": [1.0, 0.0, 0.0],
              "real_space_b": [0.0, 1.0, 0.0], "real_space_c": [0.0, 0.0, 2.0],
              "space_group_hall_symbol": " C 2y"}
  crystal = Crystal.from_dict(exp_dict)
  scan = Scan(image_range=[0, 90], oscillation=[0.0, 1.0])
  beam = Beam(s0=(0.0, 0.0, 1.01))
  goniometer = Goniometer((1.0, 0.0, 0.0))
  detector = Detector()
  experiments.append(Experiment(beam=beam, scan=scan, goniometer=goniometer,
    detector=detector, crystal=crystal))
  return experiments

#@pytest.fixture(scope='module')
def generated_two_exp():
  """Generate a two experiments Experimentlist."""
  experiments = ExperimentList()
  exp_dict = {"__id__" : "crystal", "real_space_a": [1.0, 0.0, 0.0],
              "real_space_b": [0.0, 1.0, 0.0], "real_space_c": [0.0, 0.0, 2.0],
              "space_group_hall_symbol": " C 2y"}
  crystal = Crystal.from_dict(exp_dict)
  scan = Scan(image_range=[0, 90], oscillation=[0.0, 1.0])
  beam = Beam(s0=(0.0, 0.0, 1.01))
  goniometer = Goniometer((1.0, 0.0, 0.0))
  goniometer_2 = Goniometer((1.0, 1.0, 0.0))
  detector = Detector()
  experiments.append(Experiment(beam=beam, scan=scan, goniometer=goniometer,
    detector=detector, crystal=crystal))
  experiments.append(Experiment(beam=beam, scan=scan, goniometer=goniometer_2,
    detector=detector, crystal=crystal))
  return experiments

@pytest.fixture(scope='module')
def generated_param():
  """Generate a params phil scope."""
  phil_scope = phil.parse('''
      include scope dials_scratch.jbe.scaling_code.scaling_options.phil_scope
  ''', process_includes=True)

  optionparser = OptionParser(phil=phil_scope, check_format=False)
  parameters, _ = optionparser.parse_args(args=None, quick_parse=True,
    show_diff_phil=False)
  parameters.__inject__('model', 'KB')
  return parameters

def generated_target_input(generated_refl, generated_refl_2,
    generated_single_exp, generated_param):
  """Generate suitable input for a targetscaler."""
  refl = generated_refl
  refl_2 = generated_refl_2
  refl.append(refl_2[0])
  exp = generated_single_exp
  param = generated_param
  #refl.append(refl_2)
  return (refl, exp, param)

def test_TargetScalerFactory():
  """Test for successful initialisation of TargetedScaler."""
  (test_reflections, test_experiments, params) = generated_target_input(
    generated_refl(), generated_refl_2(), generated_two_exp(), generated_param())

  assert len(test_reflections) == 2
  experiments = create_scaling_model(params, test_experiments, test_reflections)
  ss = SingleScalerBase(params, experiments[0], test_reflections[0])
  ss2 = SingleScalerBase(params, experiments[0], test_reflections[1])
  targetscaler = TargetScaler(params, experiments, [ss],
    experiments, [ss2])
  assert isinstance(targetscaler, TargetScaler)
  assert list(targetscaler.Ih_table.Ih_values) == list(flex.double([10.0, 1.0, 1.0]))
  assert list(targetscaler.unscaled_scalers[0].Ih_table.asu_miller_index) == (
    list(flex.miller_index([(1, 0, 0)])))
  assert list(targetscaler.unscaled_scalers[0].Ih_table.Ih_values) == list(flex.double([1.0]))

  # Test fixed Ih scaling target
  apm_factory = create_apm(targetscaler)
  apm = apm_factory.make_next_apm()

  # Create a scaling target and check gradients and residuals.
  target_function = ScalingTargetFixedIH(targetscaler, apm)
  target_function.predict()
  resid = target_function.calculate_residuals()
  assert list(resid) == [1.0]#   (weight x (2.0 - (1.0 * 1.0))**2)
  gradients = target_function.calculate_gradients()
  fd_grad = calculate_gradient_fd(target_function)
  assert approx_equal(list(fd_grad), list(gradients))
