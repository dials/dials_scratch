#!/usr/bin/env python

# LIBTBX_SET_DISPATCHER_NAME dev.dials.plot_Fo_vs_Fc

"""
Create a plot of Fo vs Fc similar to that shown by Figure 6 in
https://doi.org/10.1107/S2059798317010348

Usage: dev.dials.plot_Fo_vs_Fc hklin=refined.mtz
"""

from __future__ import division, print_function, absolute_import
import sys
from libtbx.utils import Sorry
from dials.util.options import OptionParser
#from libtbx.table_utils import simple_table
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from iotbx import mtz
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

      """

  from scitbx.lstbx import normal_eqns
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

class Script(object):
  '''A class for running the script.'''

  def __init__(self):
    '''Initialise the script.'''
    from libtbx.phil import parse
    import libtbx.load_env

    # The phil scope
    phil_scope = parse('''
      hklin = None
        .type = path
        .help = "MTZ file (containing observed and calculated structure factors)"

      Fo = F
        .type = str
        .help = "MTZ column name for Fobs"

      Fc = FC
        .type = str
        .help = "MTZ column name for Fcalc"

      max_Fc = 300
        .type = float
        .help = "Plot and perform fit with data up this value of Fc"

      plot_filename = Fo_vs_Fc.pdf
        .type = str
        .help = "Filename for plot"

      fit_hyperbola = True
        .type = bool
        .help = "Calculate and show the fit of a hyperbolic function given by"
                "|Fo|^2 = |Fc|^2 + |Fe|^2, where |Fe| describes the error term"
                "containing information about dynamic scattering and other"
                "effects"

      show_y_eq_x = True
        .type = bool
        .help = "Plot y=x as a dashed line"
    ''', process_includes=True)

    # The script usage
    #import __main__
    usage = ("usage: dev.dials.plot_Fo_vs_Fc hklin=refined.mtz")

    # Create the parser
    self.parser = OptionParser(
      usage=usage,
      phil=phil_scope,
      epilog=__doc__)

    self.model_fit=None

    return

  def _extract_data_from_mtz(self):
    try:
      m = mtz.object(self.params.hklin)
    except RuntimeError:
      raise Sorry('Could not read {0}'.format(self.params.hklin))

    try:
      fobs = m.extract_reals(self.params.Fo)
      fc = m.extract_reals(self.params.Fc)
    except RuntimeError:
      raise Sorry('Columns {0} not found in available labels: {1}'.format(
        ", ".join([self.params.Fo, self.params.Fc]),
        ", ".join(m.column_labels())))

    self.fobs = fobs.data
    self.fc = fc.data

    if self.params.max_Fc:
      sel = self.fc <= self.params.max_Fc
      self.fobs = self.fobs.select(sel)
      self.fc = self.fc.select(sel)

    return

  def _plot(self):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    minor_loc = MultipleLocator(10)
    ax.yaxis.set_minor_locator(minor_loc)
    ax.xaxis.set_minor_locator(minor_loc)
    ax.grid(True, which='minor')
    ax.set_axisbelow(True)
    ax.set_aspect('equal')
    ax.set_xlabel(r'$F_c$')
    ax.set_ylabel(r'$F_o$')
    ax.scatter(self.fc, self.fobs, s=1, c="indianred")

    if self.params.max_Fc:
      ax.set_xlim((0, self.params.max_Fc))
      ax.set_ylim((0, self.params.max_Fc))

    if self.params.show_y_eq_x:
      ax.plot(ax.get_xlim(), ax.get_ylim(), ls="--", c="0.0", linewidth=0.1)

    if self.model_fit:
      x = flex.double_range(0, int(ax.get_xlim()[1]))
      y = self.model_fit(x)
      ax.plot(x, y, c="0.0", linewidth=0.1)

    print("Saving plot to {0}".format(self.params.plot_filename))
    plt.savefig(self.params.plot_filename)

  def run(self):
    '''Execute the script.'''

    # Parse the command line
    self.params, options = self.parser.parse_args(show_diff_phil=True)

    if self.params.hklin is None:
      self.parser.print_help()
      sys.exit()

    self._extract_data_from_mtz()

    if self.params.fit_hyperbola:
      # The hyperbola formula has a single parameter, |Fe|, which can be
      # determined as the intercept of a linear regression line after a
      # change of variables
      x = flex.pow2(self.fc)
      y = flex.pow2(self.fobs)
      intercept = linear_regression(x, (y - x), fit_gradient=False)[0]

      # Set the model_fit function using the determined intercept
      def hyperbola(x, c):
        return flex.sqrt(flex.pow2(x) + c)
      from functools import partial
      self.model_fit = partial(hyperbola, c=intercept)

    if self.params.plot_filename:
      self._plot()

    return

def test_linear_fit():

  speed = flex.double([4, 4, 7, 7, 8, 9, 10, 10, 10, 11])
  dist = flex.double([2,10,4,22,16,10,18,26,34,17])

  import pytest
  assert linear_regression(dist, speed) == pytest.approx([5.40711598, 0.16307447])

def test_hyperbola_fit():

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

if __name__ == '__main__':
  from dials.util import halraiser

  test_hyperbola_fit()
  test_linear_fit()

  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
