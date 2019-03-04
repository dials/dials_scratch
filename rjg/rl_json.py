# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

from __future__ import division
from __future__ import print_function

import libtbx.load_env
import logging
logger = logging.getLogger(libtbx.env.dispatcher_name)

import libtbx.phil
from scitbx.array_family import flex

phil_scope= libtbx.phil.parse("""
include scope dials.command_line.reciprocal_lattice_viewer.phil_scope
output {
  json = rlp.json
    .type = path
  compact = True
    .type = bool
}
""", process_includes=True)

help_message = """\
"""

from dials.command_line.reciprocal_lattice_viewer import render_3d

class ReciprocalLatticeJson(render_3d):

  def __init__(self, settings=None):
    render_3d.__init__(self)
    if settings is not None:
      self.settings = settings
    else:
      self.settings = settings()

  def load_models(self, imagesets, reflections):
    self.imagesets = imagesets
    self.reflections_input = reflections
    self.map_points_to_reciprocal_space()

  def as_dict(self):
    rlp = list(self.reflections['rlp'])
    if 'imageset_id' in self.reflections:
      imageset_id = list(self.reflections['imageset_id'])
      expt_id = list(self.reflections['id'])
    else:
      imageset_id = list(self.reflections['id'])
      expt_id = None

    indexed = self.reflections.get_flags(self.reflections.flags.indexed)
    d = {
      'rlp': rlp,
      'imageset_id': imageset_id,
      'experiment_id': expt_id,
    }
    return d

  def as_json(self, filename=None, compact=False):
    import json
    d = self.as_dict()
    if compact:
      text = json.dumps(d, separators=(',',':'), ensure_ascii=True)
    else:
      text = json.dumps(d, indent=2, ensure_ascii=True)

    if filename is not None:
      from libtbx import smart_open
      with smart_open.for_writing(filename) as f:
        f.write(text)
    else:
      return text

def run(args):

  from dials.util.options import OptionParser
  from dials.util.options import flatten_datablocks
  from dials.util.options import flatten_experiments
  from dials.util.options import flatten_reflections
  from dials.util import log

  usage = "%s [options] datablock.json reflections.pickle" %(
    libtbx.env.dispatcher_name)

  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    read_datablocks=True,
    read_experiments=True,
    read_reflections=True,
    check_format=False,
    epilog=help_message)

  params, options = parser.parse_args()
  datablocks = flatten_datablocks(params.input.datablock)
  experiments = flatten_experiments(params.input.experiments)
  reflections = flatten_reflections(params.input.reflections)

  if (len(datablocks) == 0 and len(experiments) == 0) or len(reflections) == 0:
    parser.print_help()
    exit(0)

  ## Configure the logging
  #log.config(info='dials.rl_png.log')

  # Log the diff phil
  diff_phil = parser.diff_phil.as_str()
  if diff_phil is not '':
    logger.info('The following parameters have been modified:\n')
    logger.info(diff_phil)

  reflections = reflections[0]

  if len(datablocks) == 0 and len(experiments) > 0:
    imagesets = experiments.imagesets()
  else:
    imagesets = []
    for datablock in datablocks:
      imagesets.extend(datablock.extract_imagesets())

  f = ReciprocalLatticeJson(settings=params)
  f.load_models(imagesets, reflections)
  f.as_json(filename=params.output.json, compact=params.output.compact)
  print()


if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
