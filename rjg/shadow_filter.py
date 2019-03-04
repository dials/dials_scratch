from __future__ import division
from __future__ import print_function
import copy
from dials.array_family import flex
import libtbx.phil

help_message = '''

'''

phil_scope= libtbx.phil.parse('''
negate = False
  .type = bool

output {
  reflections = shadowed.pickle
    .type = path
  filter_hkl = None
    .type = path
}
''')


def run(args):

  from dials.util.options import OptionParser
  from dials.util.options import flatten_experiments
  from dials.util.options import flatten_reflections
  import libtbx.load_env

  usage = "%s [options] datablock.json" %(
    libtbx.env.dispatcher_name)

  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    read_experiments=True,
    read_reflections=True,
    check_format=True,
    epilog=help_message)

  params, options = parser.parse_args(show_diff_phil=True)
  experiments = flatten_experiments(params.input.experiments)
  reflections = flatten_reflections(params.input.reflections)

  if (len(reflections) == 0 or len(experiments) == 0):
    parser.print_help()
    exit(0)

  reflections = reflections[0]
  assert len(experiments) == 1

  experiment = experiments[0]

  from dials.command_line.check_strategy import filter_shadowed_reflections
  sel = filter_shadowed_reflections(experiments, reflections)
  print("%i/%i (%.2f%%) shadowed reflections" %(
    sel.count(True), sel.size(), 100*sel.count(True)/sel.size()))

  if params.negate:
    sel = ~sel
  shadowed = reflections.select(sel)
  shadowed.as_pickle(params.output.reflections)

  if params.output.filter_hkl is not None:

    from rstbx.cftbx.coordinate_frame_helpers import align_reference_frame
    from scitbx import matrix
    detector = experiment.detector

    if len(detector) > 1:
      fast = detector[0].get_parent_fast_axis()
      slow = detector[0].get_parent_slow_axis()
      Rd = align_reference_frame(fast, (1,0,0), slow, (0,1,0))
      origin = Rd * matrix.col(detector[0].get_parent_origin())
    else:
      fast = detector[0].get_fast_axis()
      slow = detector[0].get_slow_axis()
      Rd = align_reference_frame(fast, (1,0,0), slow, (0,1,0))
      origin = Rd * matrix.col(detector[0].get_origin())

    with open(params.output.filter_hkl, 'wb') as f:

      for ref in shadowed:
        p = detector[ref['panel']]
        ox, oy = p.get_raw_image_offset()
        h, k, l = ref['miller_index']
        x, y, z = ref['xyzcal.px']
        dx, dy, dz = (2, 2, 2)
        print("%i %i %i %.1f %.1f %.1f %.1f %.1f %.1f" %(
          h, k, l, x+ox, y+oy, z, dx, dy, dz), file=f)


if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
