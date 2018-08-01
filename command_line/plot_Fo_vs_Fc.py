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

def linear_fit(y, x):
  from scitbx.lstbx import normal_eqns
  from scitbx import sparse
  eqns = normal_eqns.linear_ls(2)
  rhs = y
  nobs = len(rhs)

  # design matrix has to be sparse
  a = sparse.matrix(nobs, 2)
  beta0 = flex.double(nobs, 1)
  beta0.reshape(flex.grid(nobs,1))
  a.assign_block(beta0, 0, 0)
  beta1 = flex.double(x)
  beta1.reshape(flex.grid(nobs,1))
  a.assign_block(beta1, 0, 1)

  w = flex.double(nobs, 1)

  eqns.add_equations(right_hand_side=rhs, design_matrix=a, weights=w)
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

      max_F = 500
        .type = float
        .help = "Set plot limits to include up to this value"

      fit_line = linear hyperbolic
        .type = choice
        .help = "Type of model for a least-squares fit"

      plot_filename = Fo_vs_Fc.pdf
        .type = str
        .help = "Filename for plot"
    ''', process_includes=True)

    # The script usage
    import __main__
    usage = ("usage: dev.dials.plot_Fo_vs_Fc hklin=refined.mtz")

    # Create the parser
    self.parser = OptionParser(
      usage=usage,
      phil=phil_scope,
      epilog=__doc__)

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
    if self.params.max_F:
      ax.set_xlim((0, self.params.max_F))
      ax.set_ylim((0, self.params.max_F))
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

    print(linear_fit(self.fobs, self.fc))

    if self.params.plot_filename:
      self._plot()

    return

if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
