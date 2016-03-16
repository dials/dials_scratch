from __future__ import division

import iotbx.phil
from dials.util.options import OptionParser

help_message = '''
'''

phil_scope = iotbx.phil.parse("""
d_min = None
  .type = float(value_min=0)
""", process_includes=True)


def run(args):
  import libtbx.load_env
  usage = "%s [options]" %libtbx.env.dispatcher_name

  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    check_format=False,
    epilog=help_message)

  params, options, args = parser.parse_args(show_diff_phil=True,
                                            return_unhandled=True)

  assert len(args) == 1
  from iotbx.reflection_file_reader import any_reflection_file

  intensities = None

  f = args[0]

  arrays = any_reflection_file(f).as_miller_arrays(merge_equivalents=False)
  for ma in arrays:
    print ma.info().labels
    if ma.info().labels == ['I', 'SIGI']:
      intensities = ma
    elif ma.info().labels == ['IMEAN', 'SIGIMEAN']:
      intensities = ma
    elif ma.info().labels == ['I(+)', 'SIGI(+)', 'I(-)', 'SIGI(-)']:
      intensities = ma

  assert intensities is not None

  if params.d_min is not None:
    intensities = intensities.resolution_filter(d_min=params.d_min)

  from cctbx.array_family import flex

  # see also:
  #   cctbx/miller/merge_equivalents.h
  #   cctbx/miller/equivalent_reflection_merging.tex

  # this should calculate the external variance, i.e. V(y) = sum(v_i)
  merging_external = intensities.merge_equivalents(use_internal_variance=False)
  multiplicities = merging_external.redundancies().data()
  external_sigmas = merging_external.array().sigmas()
  # sigmas should be bigger not smaller
  external_sigmas *= flex.sqrt(multiplicities.as_double())

  # set the sigmas to 1, and calculate the mean intensities and internal variances
  intensities_copy = intensities.customized_copy(
    sigmas=flex.double(intensities.size(), 1))
  merging_internal = intensities_copy.merge_equivalents()
  merged_intensities = merging_internal.array()
  internal_sigmas = merging_internal.array().sigmas()
  # sigmas should be bigger not smaller
  internal_sigmas *= flex.sqrt(multiplicities.as_double())

  # select only those reflections with sufficient repeat observations
  sel = (multiplicities > 3)
  external_sigmas = external_sigmas.select(sel)
  internal_sigmas = internal_sigmas.select(sel)
  merged_intensities = merged_intensities.select(sel)

  # what we want to plot/do linear regression with
  y = flex.pow2(internal_sigmas/merged_intensities.data())
  x = flex.pow2(external_sigmas/merged_intensities.data())

  sel = (x < 1) & (y < 1)
  x = x.select(sel)
  y = y.select(sel)

  # set backend before importing pyplot
  import matplotlib
  #matplotlib.use('Agg')

  linreg = flex.linear_regression(x, y)
  linreg.show_summary()
  import math
  print 1/math.sqrt(linreg.slope() * linreg.y_intercept())

  #x = -flex.log10(x)
  #y = -flex.log10(y)

  x = 1/x
  y = 1/y

  from matplotlib import pyplot
  pyplot.scatter(x, y, marker='+', s=20, alpha=1, c='black')
  pyplot.show()
  pyplot.clf()

  # chi^2 plot vs resolution
  # i.e. <var(int)>/<var(ext)>
  # where var(ext) and var(int) are as defined in equations 4 & 5 respectively
  # in Blessing (1997)

  internal_var = merged_intensities.customized_copy(
    data=flex.pow2(internal_sigmas))
  external_var = merged_intensities.customized_copy(
    data=flex.pow2(external_sigmas))

  n_bins = 10
  internal_var.setup_binner(n_bins=n_bins)
  external_var.use_binning_of(internal_var)

  mean_internal = internal_var.mean(use_binning=True)
  mean_external = external_var.mean(use_binning=True)

  y = [mean_internal.data[i+1]/mean_external.data[i+1] for i in range(n_bins)]
  x = [mean_internal.binner.bin_centers(2)]

  pyplot.scatter(x, y)
  pyplot.xlabel('1/d^2')
  pyplot.ylabel('<var(int)>/<var(ext)>')
  pyplot.show()
  pyplot.clf()

  return


if __name__ == '__main__':
  import sys
  from libtbx.utils import show_times_at_exit
  show_times_at_exit()
  run(sys.argv[1:])
