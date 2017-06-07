from __future__ import absolute_import, division
from cctbx.array_family import flex

import iotbx.phil

help_message = '''
'''


phil_scope = iotbx.phil.parse(
'''
space_group = None
  .type = space_group
max_deviation = 5
  .type = float(value_min=0)
  .help = "Maximum angular deviation from previous experiments orientation matrix (in degrees)"
output {
  prefix = reindexed_experiments_
    .type = str
}
''')


def run(args):

  from dials.util.options import OptionParser
  from dials.util.options import flatten_experiments
  import libtbx.load_env

  usage = "%s [options] experiments.json" %libtbx.env.dispatcher_name

  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    read_experiments=True,
    check_format=False,
    epilog=help_message)

  params, options = parser.parse_args(show_diff_phil=True)
  experiments = flatten_experiments(params.input.experiments)
  if len(experiments) <= 1:
    parser.print_help()
    return

  from dials.algorithms.indexing.compare_orientation_matrices import \
       difference_rotation_matrix_axis_angle
  crystals = []
  for experiment in experiments:
    crystal = experiment.crystal
    if params.space_group is not None:
      crystal.set_space_group(params.space_group.group())
    crystals.append(crystal)

  angles = flex.double()

  import math
  padding = int(math.ceil(math.log10(len(experiments))))
  output_template = '%s%%0%ii.json' %(params.output.prefix, padding)

  prev_expt = experiments[0]
  for i in range(1, len(experiments)):
    R_ij, axis, angle, cb_op = difference_rotation_matrix_axis_angle(
      prev_expt.crystal, experiments[i].crystal)
    angles.append(angle)
    #print i, angle
    if abs(angle) > params.max_deviation:
      continue
    experiments[i].crystal = experiments[i].crystal.change_basis(cb_op)
    prev_expt = experiments[i]

    from dxtbx.serialize import dump
    dump.experiment_list(experiments[i:i+1], output_template %i)

  from matplotlib import pyplot
  n, bins, patches = pyplot.hist(angles.as_numpy_array(), 100)
  pyplot.show()

if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
