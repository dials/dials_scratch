'''
This code tests for Ih_table and joint_Ih_table data structures.
This also provides a test for the scaler, which must
be successfully initialised in order to provide a feed in for the
Ih_table.
'''
import pytest
from Ih_table import SingleIhTable, JointIhTable
from dials.array_family import flex
from dials.util.options import OptionParser
from libtbx import phil
from dxtbx.model import Crystal, Scan, Beam, Goniometer, Detector, Experiment, ExperimentList
import ScalingModelFactory as ScalingModelFactory
import ScalerFactory as ScalerFactory

@pytest.fixture(scope='module')
def generate_refl_1():
  '''generate a test reflection table'''
  reflections = flex.reflection_table()
  reflections['intensity.prf.value'] = flex.double(
    [100.0, 100.0, 80.0, 60.0, 30.0, 40.0, 60.0])
  reflections['intensity.prf.variance'] = flex.double(
    [90.0, 100.0, 90.0, 60.0, 30.0, 50.0, 50.0])
  #note the variance values is what should come out as Ih_values if unity weights
  reflections['miller_index'] = flex.miller_index(
    [(1, 0, 0), (0, 0, 1), (-1, 0, 0), (0, 2, 0), (0, 4, 0), (0, 0, -2), (0, 0, 2)])
  reflections['d'] = flex.double([5.0, 5.0, 5.0, 5.0, 2.5, 2.5, 2.5])
  reflections['lp'] = flex.double([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
  reflections['dqe'] = flex.double([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
  reflections['partiality'] = flex.double([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
  reflections.set_flags(flex.bool([True, True, True, True, True, True, True]),
    reflections.flags.integrated)
  return [reflections]

@pytest.fixture(scope='module')
def generate_refl_2():
  '''generate another test reflection table'''
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
  '''generate a test experiments object'''
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
  '''generate a test params object'''
  phil_scope = phil.parse('''
      include scope dials_scratch.jbe.scaling_code.scaling_options.phil_scope
  ''', process_includes=True)
  optionparser = OptionParser(phil=phil_scope, check_format=False)
  parameters, _ = optionparser.parse_args(args=None, quick_parse=True,
    show_diff_phil=False)
  parameters.__inject__('scaling_model', 'KB')
  return parameters

@pytest.fixture(scope='module')
def generate_test_scaler(generate_refl_1, generate_experiments, generate_params):
  '''generate a test scalerr'''
  params, test_experiments, test_reflections = generate_params, generate_experiments, generate_refl_1
  experiments = ScalingModelFactory.Factory.create(params, test_experiments, test_reflections)
  scaler = ScalerFactory.Factory.create(params, experiments, test_reflections)
  return scaler

@pytest.fixture(scope='module')
def generate_second_test_scaler(generate_refl_2, generate_experiments, generate_params):
  '''generate a second test scaler'''
  params, test_experiments, test_reflections = generate_params, generate_experiments, generate_refl_2
  experiments = ScalingModelFactory.Factory.create(params, test_experiments, test_reflections)
  scaler = ScalerFactory.Factory.create(params, experiments, test_reflections)
  return scaler

@pytest.fixture
def single_test_input(generate_test_scaler):
  '''generate input for testing Ih_table'''
  scaler = generate_test_scaler
  weights = flex.double([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
  return scaler.reflection_table, weights

@pytest.fixture
def joint_test_input(generate_test_scaler, generate_second_test_scaler):
  '''generate input for testing joint_Ih_table'''
  scaler = generate_test_scaler
  weights = flex.double([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
  scaler_2 = generate_second_test_scaler
  weights_2 = flex.double([1.0, 1.0])
  scaler._Ih_table = SingleIhTable(scaler.reflection_table, weights)
  scaler_2._Ih_table = SingleIhTable(scaler_2.reflection_table, weights_2)
  return scaler, scaler_2


def test_Ih_table(single_test_input):
  '''test for Ih_table. Upon initialisation, Ih_table should set unity
  scale factors and calculate Ih_values. It should also create the
  h_index arrays and h_index_matrix'''
  (reflection_table, weights) = single_test_input
  Ih_table = SingleIhTable(reflection_table, weights)

  assert Ih_table.size == 7
  assert (Ih_table.inverse_scale_factors == 1.0).count(True) == Ih_table.size
  assert (Ih_table.weights == 1.0).count(True) == Ih_table.size
  assert list(Ih_table.asu_miller_index) == list(flex.miller_index(
    [(0, 0, 1), (0, 0, 2), (0, 0, 2), (0, 2, 0), (0, 4, 0), (1, 0, 0), (1, 0, 0)]))
  assert list(Ih_table.h_index_counter_array) == list(flex.int([1, 2, 1, 1, 2]))
  assert list(Ih_table.h_index_cumulative_array) == list(flex.int([0, 1, 3, 4, 5, 7]))
  assert list(Ih_table.n_h) == list(flex.double([1, 2, 2, 1, 1, 2, 2]))
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

  #test selection function
  sel = flex.bool([True, True, False, False, False, False, False])
  Ih_table = Ih_table.select(sel)
  assert Ih_table.size == 2
  assert (Ih_table.inverse_scale_factors == 1.0).count(True) == Ih_table.size
  assert (Ih_table.weights == 1.0).count(True) == Ih_table.size
  assert list(Ih_table.asu_miller_index) == list(flex.miller_index(
    [(0, 0, 1), (0, 0, 2)]))
  assert list(Ih_table.h_index_counter_array) == list(flex.int([1, 1]))
  assert list(Ih_table.h_index_cumulative_array) == list(flex.int([0, 1, 2]))
  assert list(Ih_table.n_h) == list(flex.double([1, 1]))
  assert list(Ih_table.Ih_values) == list(flex.double([100.0, 50.0]))
  assert list(Ih_table.variances) == list(flex.double([100.0, 50.0]))
  assert list(Ih_table.intensities) == list(flex.double([100.0, 40.0]))
  assert Ih_table.h_index_matrix[0, 0] == 1
  assert Ih_table.h_index_matrix[1, 1] == 1
  assert Ih_table.h_index_matrix.non_zeroes == 2
  assert Ih_table.h_index_matrix.n_cols == 2
  assert Ih_table.h_index_matrix.n_rows == 2

  '''test for functionality of setting weights or not'''
  Ih_table = SingleIhTable(reflection_table)
  #test that weights are set to inverse variances if no weights are given.
  expected_weights = 1.0/flex.double([100.0, 50.0, 50.0, 60.0, 30.0, 90.0, 90.0])
  assert list(Ih_table.weights) == list(expected_weights)
  #now test that one can set the weights to unity
  Ih_table.weights = weights
  assert (Ih_table.weights == 1.0).count(True) == Ih_table.size
  #now test that one can apply an error model, with params that reset w to 1/var
  Ih_table.update_aimless_error_model([1.0, 0.0])
  assert (list(abs(Ih_table.weights - expected_weights)) <
    list(flex.double([1e-6]*Ih_table.size)))

  '''test for functionality of setting Ih_values - set to a tenth of what they
  would otherwise be. Test that these are set after Ih table is initialised'''
  reflection_table['Ih_values'] = flex.double([10.0, 5.0, 5.0, 6.0, 3.0, 9.0, 9.0])
  Ih_table = SingleIhTable(reflection_table, weights)
  assert list(Ih_table.Ih_values) == list(flex.double(
    [10.0, 5.0, 5.0, 6.0, 3.0, 9.0, 9.0]))


def test_joint_Ih_table(joint_test_input):
  '''test that the two reflection tables have been sorted/combined correctly'''
  (dm1, dm2) = joint_test_input
  Ih_table = JointIhTable([dm1, dm2])

  assert list(Ih_table.asu_miller_index) == list(flex.miller_index([(0, 0, 1),
    (0, 0, 2), (0, 0, 2), (0, 2, 0), (0, 4, 0), (0, 4, 0), (1, 0, 0), (1, 0, 0),
    (1, 0, 0)]))
  assert list(Ih_table.Ih_values) == list(flex.double(
    [100.0, 50.0, 50.0, 60.0, 30.0, 30.0, 80.0, 80.0, 80.0]))
  assert Ih_table.size == 9
  assert list(Ih_table.h_index_counter_array) == list(flex.int([1, 2, 1, 2, 3]))
  assert list(Ih_table.h_index_cumulative_array) == list(flex.int([0, 1, 3, 4, 6, 9]))
  assert list(Ih_table._h_idx_count_list[0]) == list(flex.int([1, 2, 1, 1, 2]))
  assert list(Ih_table._h_idx_count_list[1]) == list(flex.int([0, 0, 0, 1, 1]))
  assert list(Ih_table._h_idx_cumulative_list[0]) == list(flex.int([0, 1, 3, 4, 5, 7]))
  assert list(Ih_table._h_idx_cumulative_list[1]) == list(flex.int([0, 0, 0, 0, 1, 2]))

  #now test h_expand matrices
  assert Ih_table.h_index_expand_list[0][0, 0] == 1
  assert Ih_table.h_index_expand_list[0][1, 1] == 1
  assert Ih_table.h_index_expand_list[0][2, 2] == 1
  assert Ih_table.h_index_expand_list[0][3, 3] == 1
  assert Ih_table.h_index_expand_list[0][4, 4] == 1
  assert Ih_table.h_index_expand_list[0][5, 6] == 1
  assert Ih_table.h_index_expand_list[0][6, 7] == 1
  assert Ih_table.h_index_expand_list[0].non_zeroes == 7
  assert Ih_table.h_index_expand_list[0].n_cols == 9
  assert Ih_table.h_index_expand_list[0].n_rows == 7

  assert Ih_table.h_index_expand_list[1][0, 5] == 1
  assert Ih_table.h_index_expand_list[1][1, 8] == 1
  assert Ih_table.h_index_expand_list[1].non_zeroes == 2
  assert Ih_table.h_index_expand_list[1].n_cols == 9
  assert Ih_table.h_index_expand_list[1].n_rows == 2
