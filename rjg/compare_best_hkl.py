from __future__ import division

from libtbx.phil import command_line
from libtbx.utils import Sorry
import iotbx.phil
from cctbx import sgtbx, uctbx
# from dials.util.command_line import Importer
from dials.util.options import OptionParser
from dials.array_family import flex

help_message = '''
'''

phil_scope = iotbx.phil.parse('''
space_group = None
  .type = space_group
unit_cell = None
  .type = unit_cell
''', process_includes=True)


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

  space_group = params.space_group
  if space_group is None:
    space_group = sgtbx.space_group()
  else:
    space_group = space_group.group ()

  unit_cell = params.unit_cell
  if unit_cell is None:
    unit_cell = space_group.info().any_compatible_unit_cell(volume=100000)
    print unit_cell

  assert len(args) == 2
  from cctbx import crystal, miller
  cs = crystal.symmetry(
    space_group=space_group, unit_cell=unit_cell)
  intensities = []
  for filename in args:
    hkl, i, sigi = parse_best_hkl(filename)
    ms = miller.set(cs, hkl)
    ma = miller.array(ms, data=i, sigmas=sigi)
    intensities.append(ma)
    #ma.show_summary()

  # Two subplots, the axes array is 1-d
  from matplotlib import pyplot

  ma1, ma2 = intensities
  hist1 = flex.histogram(ma1.data(), n_slots=100)
  hist2 = flex.histogram(ma2.data(), n_slots=100)
  f, axarr = pyplot.subplots(2, sharex=True)
  axarr[0].bar(hist1.slot_centers() - 0.5 * hist1.slot_width(), hist1.slots(), align="center",
               width=hist1.slot_width(), color='black', edgecolor=None)
  axarr[1].bar(hist2.slot_centers() - 0.5 * hist2.slot_width(), hist2.slots(), align="center",
               width=hist2.slot_width(), color='black', edgecolor=None)
  pyplot.savefig('hist_intensities.png')
  pyplot.clf()

  hist1 = flex.histogram(ma1.data()/ma1.sigmas(), n_slots=100)
  hist2 = flex.histogram(ma2.data()/ma2.sigmas(), n_slots=100)
  f, axarr = pyplot.subplots(2, sharex=True)
  axarr[0].bar(hist1.slot_centers() - 0.5 * hist1.slot_width(), hist1.slots(), align="center",
               width=hist1.slot_width(), color='black', edgecolor=None)
  axarr[1].bar(hist2.slot_centers() - 0.5 * hist2.slot_width(), hist2.slots(), align="center",
               width=hist2.slot_width(), color='black', edgecolor=None)
  pyplot.savefig('hist_isigi.png')
  pyplot.clf()

  print ma1.d_max_min()
  print ma2.d_max_min()
  ma1.setup_binner(n_bins=20)
  ma2.setup_binner(n_bins=20)

  imean1 = ma1.mean(use_binning=True)
  imean2 = ma2.mean(use_binning=True)
  f, axarr = pyplot.subplots(2, sharex=True)
  axarr[0].scatter(imean1.binner.bin_centers(2), imean1.data[1:-1])
  axarr[1].scatter(imean2.binner.bin_centers(2), imean2.data[1:-1])
  ax = pyplot.gca()
  xticks = ax.get_xticks()
  xticks_d = [
    '%.2f' %uctbx.d_star_sq_as_d(ds2) if ds2 > 0 else 0 for ds2 in xticks]
  ax.set_xticklabels(xticks_d)
  pyplot.xlabel('d spacing (A)')
  pyplot.savefig('imean_vs_resolution.png')
  pyplot.clf()

  isigi1 = ma1.i_over_sig_i(use_binning=True)
  isigi2 = ma2.i_over_sig_i(use_binning=True)
  f, axarr = pyplot.subplots(2, sharex=True)
  axarr[0].scatter(isigi1.binner.bin_centers(2), isigi1.data[1:-1])
  axarr[1].scatter(isigi2.binner.bin_centers(2), isigi2.data[1:-1])
  ax = pyplot.gca()
  xticks = ax.get_xticks()
  xticks_d = [
    '%.2f' %uctbx.d_star_sq_as_d(ds2) if ds2 > 0 else 0 for ds2 in xticks]
  ax.set_xticklabels(xticks_d)
  pyplot.xlabel('d spacing (A)')
  pyplot.savefig('isigi_vs_resolution.png')
  pyplot.clf()

  best_cb_op = None
  best_count = 0
  for i_op, op in enumerate(space_group.build_derived_reflection_intensity_group(False).all_ops()):
    if not op.t().is_zero():
      continue
    cb_op = sgtbx.change_of_basis_op(op)#.inverse())

    ma1, ma2 = intensities
    ma1, ma2 = ma1.common_sets(ma2.change_basis(cb_op))
    #print cb_op
    #print ma1.size(), ma2.size()
    if ma1.size() > best_count:
      best_cb_op = cb_op
      best_count = ma1.size()

  print "Best cb_op: %s (%i matches)" %(best_cb_op, best_count)
  ma1, ma2 = intensities
  ma1, ma2 = ma1.common_sets(ma2.change_basis(best_cb_op))

  from matplotlib import pyplot
  pyplot.scatter(ma1.data(), ma2.data(), marker='+', alpha=0.5)
  m = max(pyplot.xlim()[1], pyplot.ylim()[1])
  pyplot.plot((0,m), (0, m), c='black')
  pyplot.xlabel(args[0])
  pyplot.ylabel(args[1])
  pyplot.savefig('scatter_intensities.png')
  pyplot.clf()

  pyplot.scatter(ma1.sigmas(), ma2.sigmas(), marker='+', alpha=0.5)
  m = max(pyplot.xlim()[1], pyplot.ylim()[1])
  pyplot.plot((0,m), (0, m), c='black')
  pyplot.savefig('scatter_sigmas.png')
  pyplot.clf()

  pyplot.scatter(flex.pow2(ma1.sigmas()), flex.pow2(ma2.sigmas()),
    marker='+', alpha=0.5)
  m = max(pyplot.xlim()[1], pyplot.ylim()[1])
  pyplot.plot((0,m), (0, m), c='black')
  pyplot.savefig('scatter_variances.png')
  pyplot.clf()

  isigi1 = ma1.data()/ma1.sigmas()
  isigi2 = ma2.data()/ma2.sigmas()
  pyplot.scatter(isigi1, isigi2, marker='+', alpha=0.5)
  m = max(pyplot.xlim()[1], pyplot.ylim()[1])
  pyplot.plot((0,m), (0, m), c='black')
  pyplot.savefig('scatter_i_sig_i.png')
  pyplot.clf()

  return

def parse_best_hkl(filename):
  hkl = flex.miller_index()
  intensities = flex.double()
  sigmas = flex.double()
  with open(filename, 'rb') as f:
    for line in f.readlines():
      tokens = line.split()
      assert len(tokens) == 5, tokens
      hkl.append(tuple(int(t) for t in tokens[:3]))
      intensities.append(float(tokens[3]))
      sigmas.append(float(tokens[4]))
  return hkl, intensities, sigmas



if __name__ == '__main__':
  import sys
  from libtbx.utils import show_times_at_exit
  show_times_at_exit()
  run(sys.argv[1:])
