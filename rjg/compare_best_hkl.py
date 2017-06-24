from __future__ import division

from libtbx.phil import command_line
from libtbx.utils import Sorry
import iotbx.phil
from cctbx import sgtbx
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

  best_cb_op = None
  best_count = 0  for i_op, op in enumerate(space_group.build_derived_reflection_intensity_group(False).all_ops()):
    if not op.t().is_zero():      continue    cb_op = sgtbx.change_of_basis_op(op)#.inverse())

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

