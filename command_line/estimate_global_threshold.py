# LIBTBX_SET_DISPATCHER_NAME dev.dials.estimate_global_threshold
from __future__ import absolute_import, division
from __future__ import print_function

import iotbx.phil
from scitbx.array_family import flex
from dials.util.options import OptionParser
from dials.util.options import flatten_datablocks
from libtbx.utils import Sorry

help_message = '''
'''

phil_scope = iotbx.phil.parse("""\
plot = False
  .type = bool
""", process_includes=True)


def run(args):
  import libtbx.load_env
  from libtbx.utils import Sorry
  usage = "%s [options] datablock.json" %libtbx.env.dispatcher_name

  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    read_datablocks=True,
    check_format=True,
    read_datablocks_from_images=True,
    epilog=help_message)

  params, options = parser.parse_args(show_diff_phil=False)

  # Log the diff phil
  diff_phil = parser.diff_phil.as_str()
  if diff_phil is not '':
    print('The following parameters have been modified:\n')
    print(diff_phil)

  datablocks = flatten_datablocks(params.input.datablock)

  if len(datablocks) == 0:
    parser.print_help()
    return
  elif len(datablocks) > 1:
    raise Sorry("Only one DataBlock can be processed at a time")
  else:
    imagesets = []
    for datablock in datablocks:
      imagesets.extend(datablock.extract_imagesets())

  assert len(imagesets) == 1
  imageset = imagesets[0]
  image = imageset[0]
  assert len(image) == 1
  from dials.extensions.kabsch_spotfinder_threshold_ext import estimate_global_threshold
  threshold = estimate_global_threshold(image[0], plot=params.plot)
  print("Estimated global_threshold: %i" %(threshold))
  return


if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
