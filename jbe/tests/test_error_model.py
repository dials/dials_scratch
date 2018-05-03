"""
Tests for the error model.
"""
from math import sqrt
import pytest
from libtbx.test_utils import approx_equal
from dials_scratch.jbe.scaling_code.error_model.error_model import BasicErrorModel
from dials_scratch.jbe.scaling_code.error_model.error_model_target import ErrorModelTarget
from dials_scratch.jbe.scaling_code.Ih_table import SingleIhTable
from dials.array_family import flex
from cctbx.sgtbx import space_group

@pytest.fixture()
def large_reflection_table():
  """Create a reflection table."""
  return generate_refl_1()

@pytest.fixture(scope='module')
def test_sg():
  """Create a space group object."""
  return space_group("P 1")

def generate_refl_1():
  """Generate a test reflection table. Note tha the variance values are chosen
  as the 'True' Ih_values, which would be found if unity weights were chosen
  in this example."""
  reflections = flex.reflection_table()
  reflections['intensity'] = flex.double([1.0, 2.0, 3.0, 4.0, 5.0,
    6.0, 7.0, 8.0, 9.0, 10.0])
  reflections['inverse_scale_factor'] = flex.double(10, 1.0)
  reflections['variance'] = flex.double([1.0, 2.0, 3.0, 4.0, 5.0,
    6.0, 7.0, 8.0, 9.0, 10.0])
  reflections['miller_index'] = flex.miller_index([(1, 0, 0), (0, 0, 1),
    (-3, 0, 0), (0, 2, 0), (1, 0, 0), (0, 0, -20), (0, 0, 2), (10, 0, 0),
    (0, 0, 10), (1, 0, 0)])
  reflections.set_flags(flex.bool(10, True), reflections.flags.integrated)
  return reflections

def test_errormodel(large_reflection_table, test_sg):
  """Test the initialisation and methods of the error model."""
  Ih_table = SingleIhTable(large_reflection_table, test_sg)

  error_model = BasicErrorModel(Ih_table, n_bins=10)
  for i in range(error_model.summation_matrix.n_cols):
    assert error_model.summation_matrix[i, i] == 1
  assert error_model.summation_matrix.non_zeroes == large_reflection_table.size()
  assert error_model.bin_counts == flex.double(large_reflection_table.size(), 1)
  assert error_model.n_h == Ih_table.n_h

  # Test calc sigmaprime
  x0 = 1.0
  x1 = 0.1
  error_model.sigmaprime = error_model.calc_sigmaprime([x0, x1])
  assert error_model.sigmaprime == x0 * ((Ih_table.variances +
    ((x1*Ih_table.intensities)**2))**0.5)

  # Test calc delta_hl
  error_model.sigmaprime = error_model.calc_sigmaprime([1.0, 0.0]) #Reset
  # Calculate example for three elements, with intensities 1, 5 and 10 and
  # variances 1, 5 and 10 using he formula
  # delta_hl = sqrt(n_h - 1 / n_h) * (Ihl/ghl - Ih) / sigmaprime
  error_model.delta_hl = error_model.calc_deltahl()
  expected_deltas = [(-17.0/13.0) * sqrt(2.0/3.0), 0.0, 0.0, 0.0,
    (7.0/13.0) * sqrt(10.0/3.0), 0.0, 0.0, 0.0, 0.0,
    (10.0/13.0) * sqrt(20.0/3.0)]
  assert approx_equal(list(error_model.delta_hl), expected_deltas)

  # Test bin variance calculation on example with fewer bins.
  error_model = BasicErrorModel(Ih_table, n_bins=2)
  for i in range(5):
    assert error_model.summation_matrix[i, 0] == 1
    assert error_model.summation_matrix[i+5, 1] == 1
  error_model.sigmaprime = error_model.calc_sigmaprime([1.0, 0.0])
  error_model.delta_hl = error_model.calc_deltahl()
  bin_vars = error_model.calculate_bin_variances()
  mu_0 = sqrt(2.0/3.0) * ((7.0 * sqrt(5.0)) - 17.0)/ 65.0
  a = ((3.0 * mu_0**2) + (((7.0/13.0) * sqrt(10.0/3.0)) - mu_0)**2 +
    (((-17.0/13.0) * sqrt(2.0/3.0)) - mu_0)**2)/5.0
  expected_bin_vars = [a, 320.0/507.0]
  assert approx_equal(list(bin_vars), expected_bin_vars)

  # Test updating call
  error_model = BasicErrorModel(Ih_table, n_bins=2)
  error_model.update_for_minimisation([1.0, 0.0])
  assert approx_equal(list(error_model.bin_variances), expected_bin_vars)

def test_error_model_target(large_reflection_table, test_sg):
  """Test the error model target."""
  Ih_table = SingleIhTable(large_reflection_table, test_sg)

  error_model = BasicErrorModel(Ih_table, n_bins=2)
  error_model.update_for_minimisation([1.0, 0.05])
  target = ErrorModelTarget(error_model)
  # Test residual calculation
  residuals = target.calculate_gradients()
  assert residuals == (flex.double(2, 1.0) - error_model.bin_variances)**2

  # Test gradient calculation against finite differences.
  gradients = target.calculate_gradients()
  gradient_fd = calculate_gradient_fd(target)
  assert approx_equal(gradients, gradient_fd)

def calculate_gradient_fd(target):
  """Calculate gradient array with finite difference approach."""
  delta = 1.0e-6
  gradients = flex.double([0.0] * len(target.x))
  #iterate over parameters, varying one at a time and calculating the gradient
  for i in range(len(target.x)):
    target.x[i] -= 0.5 * delta
    target.predict()
    R_low = target.calculate_residuals()
    target.x[i] += delta
    target.predict()
    R_upper = target.calculate_residuals()
    target.x[i] -= 0.5 * delta
    target.predict()
    gradients[i] = (flex.sum(R_upper) - flex.sum(R_low)) / delta
  return gradients