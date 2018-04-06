'''
This code tests for Ih_table and joint_Ih_table data structures.
This also provides a test for the scaler, which must be successfully
initialised in order to provide input for the Ih_table.
'''
import pytest
from Ih_table import SingleIhTable, JointIhTable
from dials.array_family import flex
from dials.util.options import OptionParser
from libtbx import phil
from libtbx.test_utils import approx_equal
from dxtbx.model import Crystal, Scan, Beam, Goniometer, Detector, Experiment,\
  ExperimentList
from dials_scratch.jbe.scaling_code.model.scaling_model_factory import \
  create_scaling_model
from dials_scratch.jbe.scaling_code.scaler_factory import create_scaler

@pytest.fixture(scope='module')
def generate_refl_1():
  """Generate a test reflection table."""
  reflections = flex.reflection_table()
  reflections['intensity.prf.value'] = flex.double(
    [100.0, 100.0, 80.0, 60.0, 30.0, 40.0, 60.0])
  reflections['intensity.prf.variance'] = flex.double(
    [90.0, 100.0, 90.0, 60.0, 30.0, 50.0, 50.0])
  #note the variance values is what should come out as Ih_values if unity weights
  reflections['miller_index'] = flex.miller_index([(1, 0, 0), (0, 0, 1),
    (-1, 0, 0), (0, 2, 0), (0, 4, 0), (0, 0, -2), (0, 0, 2)])
  reflections['d'] = flex.double([5.0, 5.0, 5.0, 5.0, 2.5, 2.5, 2.5])
  reflections['lp'] = flex.double([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
  reflections['dqe'] = flex.double([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
  reflections['partiality'] = flex.double([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
  reflections.set_flags(flex.bool([True, True, True, True, True, True, True]),
    reflections.flags.integrated)
  return [reflections]

@pytest.fixture(scope='module')
def generate_refl_2():
  """Generate another test reflection table."""
  reflections = flex.reflection_table()
  reflections['intensity.prf.value'] = flex.double([60.0, 30.0])
  reflections['intensity.prf.variance'] = flex.double([60.0, 30.0])
  reflections['miller_index'] = flex.miller_index([(1, 0, 0), (0, 4, 0)])
  reflections['d'] = flex.double([5.0, 2.5])
  reflections['lp'] = flex.double([1.0, 1.0])
  reflections['dqe'] = flex.double([1.0, 1.0])
  reflections['partiality'] = flex.double([1.0, 1.0])
  reflections.set_flags(flex.bool([True, True]), reflections.flags.integrated)
  return [reflections]

@pytest.fixture(scope='module')
def generate_experiments():
  """Generate a test experiments object."""
  experiments = ExperimentList()
  exp_dict = {"__id__" : "crystal", "real_space_a": [5.0, 0.0, 0.0],
              "real_space_b": [0.0, 10.0, 0.0], "real_space_c": [0.0, 0.0, 5.0],
              "space_group_hall_symbol": " C 2y"}
  crystal = Crystal.from_dict(exp_dict)
  scan = Scan(image_range=[0, 90], oscillation=[0.0, 1.0])
  beam = Beam(s0=(0.0, 0.0, 1.01))
  goniometer = Goniometer((1.0, 0.0, 0.0))
  detector = Detector()
  experiments.append(Experiment(beam=beam, scan=scan, goniometer=goniometer,
    detector=detector, crystal=crystal))
  return experiments

@pytest.fixture(scope='module')
def generate_params():
  """Generate a test params object."""
  phil_scope = phil.parse('''
      include scope dials_scratch.jbe.scaling_code.scaling_options.phil_scope
  ''', process_includes=True)
  optionparser = OptionParser(phil=phil_scope, check_format=False)
  parameters, _ = optionparser.parse_args(args=None, quick_parse=True,
    show_diff_phil=False)
  parameters.__inject__('model', 'KB')
  return parameters

@pytest.fixture(scope='module')
def generate_test_scaler(generate_refl_1, generate_experiments, generate_params):
  """Generate a test scaler."""
  (params, test_experiments, test_reflections) = (
    generate_params, generate_experiments, generate_refl_1)
  experiments = create_scaling_model(params, test_experiments, test_reflections)
  scaler = create_scaler(params, experiments, test_reflections)
  return scaler

@pytest.fixture(scope='module')
def generate_second_test_scaler(generate_refl_2, generate_experiments,
    generate_params):
  """Generate a second test scaler."""
  (params, test_experiments, test_reflections) = (
    generate_params, generate_experiments, generate_refl_2)
  experiments = create_scaling_model(params, test_experiments, test_reflections)
  scaler = create_scaler(params, experiments, test_reflections)
  return scaler

@pytest.fixture
def single_test_input(generate_test_scaler):
  """Generate input for testing Ih_table."""
  scaler = generate_test_scaler
  weights = flex.double([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
  return scaler.reflection_table, weights

@pytest.fixture
def joint_test_input(generate_test_scaler, generate_second_test_scaler):
  """Generate input for testing joint_Ih_table."""
  scaler = generate_test_scaler
  weights = flex.double([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
  scaler_2 = generate_second_test_scaler
  weights_2 = flex.double([1.0, 1.0])
  scaler.Ih_table = SingleIhTable(scaler.reflection_table, weights)
  scaler_2.Ih_table = SingleIhTable(scaler_2.reflection_table, weights_2)
  return scaler, scaler_2

def test_Ih_table(single_test_input):
  """Test for Ih_table datastructure. Upon initialisation, Ih_table should set
  unity scale factors and calculate Ih_values. It should also create the
  a h_index_matrix."""
  (reflection_table, weights) = single_test_input
  Ih_table = SingleIhTable(reflection_table, weights)

  assert Ih_table.id_ == "IhTableBase"

  # Tests calc_Ih, assign_h_matrices, interface
  assert Ih_table.size == 7
  assert (Ih_table.inverse_scale_factors == 1.0).count(True) == Ih_table.size
  assert (Ih_table.weights == 1.0).count(True) == Ih_table.size
  assert list(Ih_table.asu_miller_index) == list(flex.miller_index([(0, 0, 1),
    (0, 0, 2), (0, 0, 2), (0, 2, 0), (0, 4, 0), (1, 0, 0), (1, 0, 0)]))
  assert list(Ih_table.Ih_values) == list(flex.double(
    [100.0, 50.0, 50.0, 60.0, 30.0, 90.0, 90.0]))
  assert list(Ih_table.variances) == list(flex.double(
    [100.0, 50.0, 50.0, 60.0, 30.0, 90.0, 90.0]))
  assert list(Ih_table.intensities) == list(flex.double(
    [100.0, 40.0, 60.0, 60.0, 30.0, 100.0, 80.0]))

  assert Ih_table.h_index_matrix[0, 0] == 1
  assert Ih_table.h_index_matrix[1, 1] == 1
  assert Ih_table.h_index_matrix[2, 1] == 1
  assert Ih_table.h_index_matrix[3, 2] == 1
  assert Ih_table.h_index_matrix[4, 3] == 1
  assert Ih_table.h_index_matrix[5, 4] == 1
  assert Ih_table.h_index_matrix[6, 4] == 1
  assert Ih_table.h_index_matrix.non_zeroes == 7
  assert Ih_table.h_index_matrix.n_cols == 5
  assert Ih_table.h_index_matrix.n_rows == 7

  # Test calc_nh function.
  Ih_table.calc_nh()
  assert list(Ih_table.n_h) == [1.0, 2.0, 2.0, 1.0, 1.0, 2.0, 2.0]

  # Test selection function
  sel = flex.bool([True, True, False, False, False, False, False])
  Ih_table = Ih_table.select(sel)
  assert Ih_table.size == 2
  assert (Ih_table.inverse_scale_factors == 1.0).count(True) == Ih_table.size
  assert (Ih_table.weights == 1.0).count(True) == Ih_table.size
  assert list(Ih_table.asu_miller_index) == list(flex.miller_index(
    [(0, 0, 1), (0, 0, 2)]))
  assert list(Ih_table.Ih_values) == list(flex.double([100.0, 50.0]))
  assert list(Ih_table.variances) == list(flex.double([100.0, 50.0]))
  assert list(Ih_table.intensities) == list(flex.double([100.0, 40.0]))
  assert Ih_table.h_index_matrix[0, 0] == 1
  assert Ih_table.h_index_matrix[1, 1] == 1
  assert Ih_table.h_index_matrix.non_zeroes == 2
  assert Ih_table.h_index_matrix.n_cols == 2
  assert Ih_table.h_index_matrix.n_rows == 2

  # Test for second method to initialise without specifying weights - weights
  # should bee set to inverse variances if no weights are given.
  Ih_table = SingleIhTable(reflection_table)
  expected_weights = 1.0/flex.double([100.0, 50.0, 50.0, 60.0, 30.0, 90.0, 90.0])
  assert list(Ih_table.weights) == list(expected_weights)
  # Test that one can set the weights to unity
  Ih_table.weights = weights
  assert (Ih_table.weights == 1.0).count(True) == Ih_table.size

  # Test that one can apply an error model, with params that reset w to 1/var
  Ih_table.update_error_model([1.0, 0.0])
  assert approx_equal(list(Ih_table.weights), list(expected_weights))

  # Test for functionality of having preset Ih_values, set to a tenth of what
  # they would be. Test that these are set after Ih table is initialised.
  reflection_table['Ih_values'] = flex.double(
    [10.0, 5.0, 5.0, 6.0, 3.0, 9.0, 9.0])
  Ih_table = SingleIhTable(reflection_table, weights)
  assert list(Ih_table.Ih_values) == list(flex.double(
    [10.0, 5.0, 5.0, 6.0, 3.0, 9.0, 9.0]))

def test_Ih_table_nonzero_weights(single_test_input):
  """Test for 'nonzero_Weights' attribute and how this changes during selection.
  The purpose of this is to indicate the relationship of the Ih_table data to
  the original input reflection table."""
  (reflection_table, weights) = single_test_input
  weights[0] = 0.0
  Ih_table = SingleIhTable(reflection_table, weights)
  assert list(Ih_table.nonzero_weights) == list(flex.bool([False, True, True,
    True, True, True, True]))
  assert Ih_table.size == 6
  Ih_table = Ih_table.select(flex.bool([True, True, True, False, False, False]))
  assert Ih_table.size == 3
  assert list(Ih_table.nonzero_weights) == list(flex.bool([False, True, True,
    True, False, False, False]))


def test_joint_Ih_table(joint_test_input):
  """Test that the joint_Ih_table datastructure correctly combined the data
  from two reflection tables."""
  (dm1, dm2) = joint_test_input
  Ih_table = JointIhTable([dm1, dm2])

  # Test for correct setup and calc_Ih method.
  assert list(Ih_table.asu_miller_index) == list(flex.miller_index(
    [(0, 0, 1), (0, 0, 2), (0, 0, 2), (0, 2, 0), (0, 4, 0), (1, 0, 0),
    (1, 0, 0), (0, 4, 0), (1, 0, 0)]))
  assert list(Ih_table.Ih_values) == list(flex.double(
    [100.0, 50.0, 50.0, 60.0, 30.0, 80.0, 80.0, 30.0, 80.0]))
  assert Ih_table.size == 9

  # Test for correct setup of joint h_index_matrix.
  assert Ih_table.h_index_matrix[0, 0] == 1
  assert Ih_table.h_index_matrix[1, 1] == 1
  assert Ih_table.h_index_matrix[2, 1] == 1
  assert Ih_table.h_index_matrix[3, 2] == 1
  assert Ih_table.h_index_matrix[4, 3] == 1
  assert Ih_table.h_index_matrix[5, 4] == 1
  assert Ih_table.h_index_matrix[6, 4] == 1
  assert Ih_table.h_index_matrix[7, 3] == 1
  assert Ih_table.h_index_matrix[8, 4] == 1
  assert Ih_table.h_index_matrix.non_zeroes == 9
  assert Ih_table.h_index_matrix.n_cols == 5
  assert Ih_table.h_index_matrix.n_rows == 9

  # Test setting of error models and updating weights.
  dm1.Ih_table.update_error_model([0.5, 0.0])
  dm2.Ih_table.update_error_model([0.5, 0.0])
  Ih_table.update_weights_from_error_models()
  expected_weights = 4.0/flex.double([100.0, 50.0, 50.0, 60.0, 30.0, 90.0,
    90.0, 30.0, 60.0])
  assert approx_equal(list(Ih_table.weights), list(expected_weights))
