from __future__ import division

import libtbx.load_env
import logging
logger = logging.getLogger(libtbx.env.dispatcher_name)

try:
  # try importing scipy.linalg before any cctbx modules, otherwise we
  # sometimes get a segmentation fault/core dump if it is imported after
  # scipy.linalg is a dependency of sklearn.cluster.DBSCAN
  import scipy.linalg # import dependency
except ImportError, e:
  pass

import copy

from libtbx.phil import command_line
import iotbx.phil
from dials.util.options import OptionParser
from dials.util.options import flatten_reflections
from dials.util.options import flatten_datablocks
from dials.array_family import flex

import matplotlib
matplotlib.use('Agg')

help_message = '''
'''

phil_scope = iotbx.phil.parse('''\
figsize = 12, 8
  .type = floats(size=2, value_min=0)
steps = 25
  .type = int(value_min=1)
''')


def run(args):
  import libtbx.load_env
  from libtbx.utils import Sorry
  from dials.util import log
  usage = "%s [options] datablock.json strong.pickle" %libtbx.env.dispatcher_name

  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    read_reflections=True,
    read_datablocks=True,
    read_datablocks_from_images=True,
    check_format=False,
    epilog=help_message)

  params, options = parser.parse_args(show_diff_phil=False)

  from dials.util.version import dials_version
  logger.info(dials_version())

  # Log the diff phil
  diff_phil = parser.diff_phil.as_str()
  if diff_phil is not '':
    logger.info('The following parameters have been modified:\n')
    logger.info(diff_phil)

  datablocks = flatten_datablocks(params.input.datablock)
  reflections = flatten_reflections(params.input.reflections)

  if len(datablocks) == 0:
    parser.print_help()
    return
  imagesets = []
  for datablock in datablocks:
    imagesets.extend(datablock.extract_imagesets())

  if len(reflections) == 0:
    raise Sorry("No reflection lists found in input")
  if len(reflections) > 1:
    assert len(reflections) == len(imagesets)
    for i in range(len(reflections)):
      reflections[i]['imageset_id'] = flex.int(len(reflections[i]), i)
      if i > 0:
        reflections[0].extend(reflections[i])

  reflections_input = reflections[0]

  for imageset in imagesets:
    if (imageset.get_goniometer() is not None and
        imageset.get_scan() is not None and
        imageset.get_scan().get_oscillation()[1] == 0):
      imageset.set_goniometer(None)
      imageset.set_scan(None)

  from dials.algorithms.indexing.indexer import indexer_base

  reflections = flex.reflection_table()

  for i, imageset in enumerate(imagesets):
    if 'imageset_id' not in reflections_input:
      reflections_input['imageset_id'] = reflections_input['id']
    sel = (reflections_input['imageset_id'] == i)
    refl = indexer_base.map_spots_pixel_to_mm_rad(
      reflections_input.select(sel),
      imageset.get_detector(), imageset.get_scan())
    indexer_base.map_centroids_to_reciprocal_space(
      refl, imageset.get_detector(), imageset.get_beam(),
      imageset.get_goniometer())
    refl['entering'] = indexer_base.calculate_entering_flags(
      refl, beam=imageset.get_beam(),
      goniometer=imageset.get_goniometer())
    reflections.extend(refl)

  filter_ice(reflections, steps=params.steps)
  return

def filter_ice(reflections, steps):

  from cctbx import miller, sgtbx, uctbx
  from matplotlib import pyplot as plt

  d_spacings = 1/reflections['rlp'].norms()
  d_star_sq = uctbx.d_as_d_star_sq(d_spacings)

  from dials.algorithms.spot_finding.per_image_analysis import \
    ice_rings_selection

  from dials.algorithms.integration import filtering

  ice_uc = uctbx.unit_cell((4.498, 4.498, 7.338,90,90,120))
  ice_sg = sgtbx.space_group_info(number=194).group()
  ice_generator = miller.index_generator(ice_uc, ice_sg.type(), False, flex.min(d_spacings))
  ice_indices = ice_generator.to_array()
  ice_d_spacings = flex.sorted(ice_uc.d(ice_indices))
  ice_d_star_sq = uctbx.d_as_d_star_sq(ice_d_spacings)

  cubic_ice_uc = uctbx.unit_cell((6.358, 6.358, 6.358,90,90,90))
  cubic_ice_sg = sgtbx.space_group_info(number=227).group()
  cubic_ice_generator = miller.index_generator(cubic_ice_uc, cubic_ice_sg.type(), False, flex.min(d_spacings))
  cubic_ice_indices = cubic_ice_generator.to_array()
  cubic_ice_d_spacings = flex.sorted(cubic_ice_uc.d(cubic_ice_indices))
  cubic_ice_d_star_sq = uctbx.d_as_d_star_sq(cubic_ice_d_spacings)

  import numpy
  widths = flex.double(numpy.geomspace(start=0.0001, stop=0.01, num=steps))
  n_spots = flex.double()
  total_intensity = flex.double()
  for width in widths:
    d_min = flex.min(d_spacings)

    ice_filter = filtering.PowderRingFilter(ice_uc, ice_sg, d_min, width)
    ice_sel = ice_filter(d_spacings)

    n_spots.append(ice_sel.count(False))
    if 'intensity.sum.value' in reflections:
      total_intensity.append(flex.sum(reflections['intensity.sum.value'].select(~ice_sel)))

  fig, axes = plt.subplots(nrows=2, figsize=(12,8), sharex=True)
  axes[0].plot(widths, n_spots, label='#spots', marker='+')
  if total_intensity.size():
    axes[1].plot(widths, total_intensity, label='total intensity', marker='+')
  axes[0].set_ylabel('# spots remaining')
  axes[1].set_xlabel('Ice ring width (1/d^2)')
  axes[1].set_ylabel('Total intensity')
  for ax in axes: ax.set_xlim(0, flex.max(widths))
  plt.savefig('ice_ring_filtering.png')
  plt.clf()
  return


if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
