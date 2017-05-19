from __future__ import absolute_import, division
from cctbx.array_family import flex

import iotbx.phil

import os

help_message = '''
'''


phil_scope = iotbx.phil.parse(
'''
filename = distances.txt
  .type = path
plot = False
  .type = bool
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
  if len(experiments) == 0:
    parser.print_help()
    return

  image_numbers = flex.size_t()
  distances = flex.double()

  for expt in experiments:
    detector = expt.detector
    assert len(detector) == 1
    panel = detector[0]
    from scitbx import matrix
    fast = matrix.col(panel.get_fast_axis())
    slow = matrix.col(panel.get_slow_axis())
    normal = fast.cross(slow)
    origin = matrix.col(panel.get_origin())
    distance = origin.dot(normal)
    fast_origin = - (origin - distance * normal).dot(fast)
    slow_origin = - (origin - distance * normal).dot(slow)
    distances.append(distance)
    imageset = expt.imageset
    path = imageset.get_path(0)
    image_number = int(os.path.splitext(path)[0].split('_')[-1])
    image_numbers.append(image_number)

  print 'Writing image number and distances to %s' %params.filename
  with open(params.filename, 'wb') as f:
    for image_number, distance in zip(image_numbers, distances):
      print >> f, image_number, distance

  if params.plot:
    from matplotlib import pyplot
    pyplot.scatter(image_numbers, distances, s=5)
    pyplot.xlabel('image number')
    pyplot.ylabel('distance')
    pyplot.show()


if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
