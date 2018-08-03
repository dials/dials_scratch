#!/usr/bin/env python

"""
An example of doing linear regression (y = mx + c curve fitting) with LSTBX
"""
from scitbx.lstbx import normal_eqns
from scitbx.array_family import flex

def linear_regression(x, y, fit_intercept=True, fit_gradient=True):
  """Perform linear regression by least squares to model how a response
  variable (y) is linearly-related to an explanatory variable (x).

  Args:
      x: Sequence of values of the explanatory variable
      y: Sequence of values of the response variable (observations)
      fit_gradient (bool): Fit the gradient as a parameter of the model
          (otherwise assume unit gradient).
      fit_intercept (bool): Fit the intercept as a parameter of the model
          (otherwise assume an intercept of zero).

  Returns:
      list: Values of the model parameters
  """
  from scitbx import sparse
  n_obs = len(y)
  assert len(x) == n_obs
  n_param = [fit_gradient, fit_intercept].count(True)
  assert n_param > 0
  eqns = normal_eqns.linear_ls(n_param)

  # The method add_equations requires a sparse matrix, so must convert to that
  # even though we build the design matrix as dense.
  a = sparse.matrix(n_obs, n_param)
  icol = 0

  if fit_intercept:
    col = flex.double(n_obs, 1)
    col.reshape(flex.grid(n_obs,1))
    a.assign_block(col, 0, icol)
    icol += 1

  if fit_gradient:
    col = flex.double(x)
    col.reshape(flex.grid(n_obs,1))
    a.assign_block(col, 0, icol)

  # Unweighted fit only, for now
  w = flex.double(n_obs, 1)

  eqns.add_equations(right_hand_side=y, design_matrix=a, weights=w)
  eqns.solve()

  return list(eqns.solution())

def test_linear_regression():
  """Tests linear fit using R's 'cars' dataset"""

  speed = flex.double([4, 4, 7, 7, 8, 9, 10, 10, 10, 11])
  dist = flex.double([2,10,4,22,16,10,18,26,34,17])

  import pytest
  assert linear_regression(dist, speed) == pytest.approx([5.40711598, 0.16307447])

def test_hyperbola_linear_fit():
  """Test fitting a hyperbola with formula y = sqrt(x^2 + c) by a change of
  variables followed by a linear fit. Warning: this is not recommended as the
  change of variables implicitly assigns weights. Better to fit the function
  by NLLS methods"""

  # Create hyperbola with intercept of 10
  c = 10
  x = flex.double([1,2,3,4,5,6,8,9,10])
  y = flex.sqrt(flex.pow2(x) + c)

  # Change of variables to do a linear fit
  u = flex.pow2(x)
  v = flex.pow2(y)

  # Perform fit of the intercept only (could equivalently just use the mean
  # of v - u)
  intercept = linear_regression(u, (v - u), fit_gradient=False)[0]
  import pytest
  assert intercept == pytest.approx(10)
