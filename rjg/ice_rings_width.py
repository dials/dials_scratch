from __future__ import division
from __future__ import print_function

import libtbx.load_env
import logging
logger = logging.getLogger(libtbx.env.dispatcher_name)

try:
  # try importing scipy.linalg before any cctbx modules, otherwise we
  # sometimes get a segmentation fault/core dump if it is imported after
  # scipy.linalg is a dependency of sklearn.cluster.DBSCAN
  import scipy.linalg # import dependency
except ImportError as e:
  pass

import copy

from libtbx.phil import command_line
import iotbx.phil
from dials.util.options import OptionParser
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
    read_datablocks=True,
    read_datablocks_from_images=True,
    check_format=True,
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

  if len(datablocks) == 0:
    parser.print_help()
    return
  imagesets = []
  for datablock in datablocks:
    imagesets.extend(datablock.extract_imagesets())

  for imageset in imagesets:
    if (imageset.get_goniometer() is not None and
        imageset.get_scan() is not None and
        imageset.get_scan().get_oscillation()[1] == 0):
      imageset.set_goniometer(None)
      imageset.set_scan(None)

  from dials.algorithms.indexing.indexer import indexer_base

  estimate_ice_rings_width(imagesets, params.steps)
  return

def estimate_ice_rings_width(imagesets, steps):

  from cctbx import miller, sgtbx, uctbx
  from matplotlib import pyplot as plt

  imageset = imagesets[0]
  detector = imageset.get_detector()
  beam = imageset.get_beam()

  from dials.util import masking
  params = masking.phil_scope.extract()
  params.resolution_range = []
  params.ice_rings.filter = True

  import numpy
  widths = flex.double(numpy.geomspace(start=0.0001, stop=0.01, num=steps))
  total_intensity = flex.double()
  n_pixels = flex.double()
  for width in widths:
    params.ice_rings.width = width
    generator = masking.MaskGenerator(params)
    mask = generator.generate(imageset)
    image = imageset.get_corrected_data(0)
    tot_intensity = 0
    n_pix = 0
    for im, m in zip(image, mask):
      im = im.as_1d()
      m = m.as_1d()
      print(m.count(True), m.count(False))
      print(flex.sum(im), flex.sum(im.select(m)), flex.sum(im.select(~m)))
      tot_intensity += flex.sum(im.select(m))
      n_pix += m.count(True)
    total_intensity.append(tot_intensity)
    n_pixels.append(n_pix)
  average_intensity = total_intensity / n_pixels

  fig, axes = plt.subplots(nrows=2, figsize=(12,8), sharex=True)
  axes[0].plot(widths, average_intensity, label='average intensity', marker='+')
  axes[1].plot(widths, total_intensity, label='total intensity', marker='+')
  axes[0].set_ylabel('Average intensity per pixel')
  axes[1].set_xlabel('Ice ring width (1/d^2)')
  axes[1].set_ylabel('Total intensity')
  for ax in axes: ax.set_xlim(0, flex.max(widths))
  plt.savefig('ice_rings_width.png')
  plt.clf()
  return


if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
