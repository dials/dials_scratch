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

help_message = '''
'''

phil_scope = iotbx.phil.parse('''\
max_cell_estimation
    .expert_level = 1
{
  filter_ice = True
    .type = bool
    .help = "Filter out reflections at typical ice ring resolutions"
            "before max_cell estimation."
  filter_overlaps = True
    .type = bool
    .help = "Filter out reflections with overlapping bounding boxes before"
            "max_cell estimation."
  overlaps_border = 0
    .type = int(value_min=0)
    .help = "Optionally add a border around the bounding boxes before finding"
            "overlaps."
  multiplier = 1.3
    .type = float(value_min=0)
    .help = "Multiply the estimated maximum basis vector length by this value."
    .expert_level = 2
  step_size = 45
    .type = float(value_min=0)
    .help = "Step size, in degrees, of the blocks used to peform the max_cell "
            "estimation."
    .expert_level = 2
  nearest_neighbor_percentile = 0.05
    .type = float(value_min=0)
    .help = "Percentile of NN histogram to use for max cell determination."
    .expert_level = 2
}
''', process_includes=True)


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

  # Configure the logging
  #log.config(
    #params.verbosity,
    #info=params.output.log,
    #debug=params.output.debug_log)

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
    #raise Sorry("Multiple reflections lists provided in input")
    assert len(reflections) == len(imagesets)
    for i in range(len(reflections)):
      reflections[i]['imageset_id'] = flex.int(len(reflections[i]), i)
      if i > 0:
        reflections[0].extend(reflections[i])

  #assert(len(reflections) == 1)
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

  params = params.max_cell_estimation
  from dials.algorithms.indexing.indexer import find_max_cell
  result = find_max_cell(
    reflections, max_cell_multiplier=params.multiplier,
    step_size=params.step_size,
    nearest_neighbor_percentile=params.nearest_neighbor_percentile,
    filter_ice=params.filter_ice,
    filter_overlaps=params.filter_overlaps,
    overlaps_border=params.overlaps_border)
  max_cell = result.max_cell
  print "Found max_cell: %.1f Angstrom" %max_cell
  result.plot_histogram()
  plot_d_spacings(reflections)
  plot_direct_space_distances(result.direct, result.d_spacings)
  #logger.info("Found max_cell: %.1f Angstrom" %(max_cell))

  return


def plot_d_spacings(reflections):
  from cctbx import miller, sgtbx, uctbx
  from matplotlib import pyplot as plt

  d_spacings = 1/reflections['rlp'].norms()
  d_star_sq = uctbx.d_as_d_star_sq(d_spacings)

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

  perm = flex.sort_permutation(d_star_sq)
  plt.plot(d_star_sq.select(perm), flex.int_range(perm.size()), zorder=10)
  for dss in ice_d_star_sq:
    plt.plot([dss, dss], plt.ylim(), c='r', zorder=0)
  for dss in cubic_ice_d_star_sq:
    plt.plot([dss, dss], plt.ylim(), c='g', zorder=1)
  ax = plt.gca()
  xticks = ax.get_xticks()
  xticks_d = [uctbx.d_star_sq_as_d(x) for x in xticks]
  ax.set_xticklabels(['%.2f' %x for x in xticks_d])
  plt.xlabel('d spacing (A^-1)')
  plt.ylabel('Cumulative frequency')
  plt.savefig('d_spacings.png')
  plt.clf()

  intensities = reflections['intensity.sum.value']
  plt.scatter(uctbx.d_as_d_star_sq(d_spacings), intensities, marker='.',
              c='black', s=1)
  ax = plt.gca()
  xticks = ax.get_xticks()
  xticks_d = [uctbx.d_star_sq_as_d(x) for x in xticks]
  ax.set_xticklabels(['%.2f' %x for x in xticks_d])
  plt.xlabel('d spacing (A^-1)')
  plt.ylabel('Intensity')
  plt.savefig('d_vs_intensity.png')
  plt.clf()

  hist = flex.histogram(d_star_sq, n_slots=200)
  plt.bar(hist.slot_centers(), hist.slots(), align="center",
          width=hist.slot_width(), zorder=10, color='black', edgecolor=None)
  for dss in ice_d_star_sq:
    plt.plot([dss, dss], plt.ylim(), c='r', zorder=0)
  for dss in cubic_ice_d_star_sq:
    plt.plot([dss, dss], plt.ylim(), c='g', zorder=0)
  ax = plt.gca()
  xticks = ax.get_xticks()
  xticks_d = [uctbx.d_star_sq_as_d(x) for x in xticks]
  ax.set_xticklabels(['%.2f' %x for x in xticks_d])
  plt.xlabel('d spacing (A^-1)')
  plt.ylabel('Frequency')
  plt.savefig('d_spacings_hist.png')
  plt.clf()

def plot_direct_space_distances(direct, d_spacings):
  from cctbx import uctbx
  from matplotlib import pyplot as plt
  plt.plot(direct, flex.double_range(direct.size())/direct.size())
  ax = plt.gca()
  ax.set_xscale('log')
  plt.xlabel('Direct space distance (A)')
  plt.ylabel('Percentile')
  plt.savefig('ordered.png')
  plt.clf()

  fig = plt.figure(figsize=(12,8))
  plt.scatter(uctbx.d_as_d_star_sq(d_spacings), direct, marker='.', c='black', s=1)
  ax = plt.gca()
  ax.set_yscale('log')
  xticks = ax.get_xticks()
  xticks_d = [uctbx.d_star_sq_as_d(x) for x in xticks]
  ax.set_xticklabels(['%.2f' %x for x in xticks_d])
  plt.ylabel('Direct space distance (A)')
  plt.xlabel('d spacing(A^-1)')
  plt.savefig('d_vs_direct.png')
  plt.clf()


if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
