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

from dials.algorithms.indexing.indexer import max_cell_phil_str

phil_scope = iotbx.phil.parse('''\
figsize = 12, 8
  .type = floats(size=2, value_min=0)
%s
''' %max_cell_phil_str)


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

  filter_ice(reflections)
  return

def filter_ice(reflections, n_bins=200):

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
  widths = flex.double(numpy.geomspace(start=0.0001, stop=0.01, num=50))
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

  hist = flex.histogram(d_star_sq, n_slots=n_bins)
  potential_ice_sel = flex.bool(n_bins, False)
  for i, (count, dss) in enumerate(zip(hist.slots(), hist.slot_centers())):
    for ice_dss in ice_d_star_sq:
      if ((dss - hist.slot_width()) < ice_dss) and ((dss + hist.slot_width()) > ice_dss):
        potential_ice_sel[i] = True

  print flex.mean(hist.slots().select(potential_ice_sel).as_double())
  non_ice_slots = hist.slots().select(~potential_ice_sel)
  print flex.mean(non_ice_slots.select(non_ice_slots > 0).as_double())

  fig = plt.figure(figsize=(12,8))
  plt.bar(hist.slot_centers(), hist.slots(), align="center",
          width=hist.slot_width(), zorder=10, color='black', edgecolor=None)
  plt.bar(hist.slot_centers().select(potential_ice_sel),
          hist.slots().select(potential_ice_sel), align="center",
          width=hist.slot_width(), zorder=10, color='red', edgecolor=None)
  ylim = plt.ylim()
  for dss in ice_d_star_sq:
    plt.plot([dss, dss], (0, ylim[1]), c='r', zorder=0, linestyle=':', alpha=0.3)
  for dss in cubic_ice_d_star_sq:
    plt.plot([dss, dss], (0, ylim[1]), c='g', zorder=1, linestyle=':', alpha=0.3)
  plt.ylim(ylim)
  ax = plt.gca()
  xticks = ax.get_xticks()
  xticks_d = [uctbx.d_star_sq_as_d(x) for x in xticks]
  ax.set_xticklabels(['%.2f' %x for x in xticks_d])
  plt.xlabel('d spacing (A^-1)')
  plt.ylabel('Frequency')
  plt.tight_layout()
  plt.savefig('d_spacings_hist.png')
  plt.clf()

if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
