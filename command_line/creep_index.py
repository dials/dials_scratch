#!/usr/bin/env cctbx.python

#
#  Copyright (C) (2017) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

# LIBTBX_SET_DISPATCHER_NAME dev.dials.creep_index


"""
Script to attempt to index a full sweep in cases where large changes in the
model during data collection (such as beam centre drift) makes it difficult
to index the sweep in one pass.

We assume an experiments.json is available that can be used to index the first
image. One this is indexed, a new experiments.json is produced with a refined
model. This experiments.json can then be used to try to index the next image,
and so on.

"""

from __future__ import division

from dxtbx.model.experiment.experiment_list import Experiment
from dxtbx.model.experiment.experiment_list import ExperimentListFactory
from dxtbx.model.experiment.experiment_list import ExperimentListDumper
from dials.util.options import flatten_reflections
from dials.array_family import flex
from dials.command_line.show import beam_centre_mm
#from scitbx import matrix
from libtbx import phil
from libtbx import easy_run
import os
from math import sqrt

#from libtbx.utils import Sorry
#import matplotlib
# Offline backend
#matplotlib.use("Agg")
#matplotlib.rc('font', family='serif')
#matplotlib.rc('font', serif='Times New Roman')
#from matplotlib import pyplot as plt

help_message = '''

Index images in a dataset sequentially.

Example::

  dev.dials.creep_index experiments.json strong.pickle

'''

phil_scope = phil.parse('''
  images_per_block=5
    .type=int
''')


class Script(object):

  def __init__(self):
    '''Check script input'''

    import libtbx.load_env
    from dials.util.options import OptionParser

    # The script usage
    usage = ("usage: {0} [options] [param.phil] experiments.json "
             "strong.pickle").format(libtbx.env.dispatcher_name)

    parser = OptionParser(
      usage=usage,
      phil=phil_scope,
      read_reflections=True,
      read_experiments=True,
      check_format=False,
      epilog=help_message)

    params, options = parser.parse_args(show_diff_phil=True)

    if len(params.input.experiments) != 1:
      raise Sorry("Please provide a single experiment list as input")

    self._current_exp_path = params.input.experiments[0].filename

    el = params.input.experiments[0].data
    scan = el.scans()
    if len(scan) != 1:
      raise Sorry("Currently only a single scan is supported")
    self._scan = scan[0]

    reflections = flatten_reflections(params.input.reflections)
    if len(reflections) > 1:
      raise Sorry("Please provide a single reflection table as input")

    self._strong_path = params.input.reflections[0].filename

    self._images_per_block = params.images_per_block
    return

  def _index_current_block(self, start, stop):

    fmt_dic = {'experiments':self._current_exp_path,
               'indexed':self._strong_path,
               'start':start,
               'stop':stop}

    cmd = ('dials.index {experiments} {indexed} '
           'scan_range={start},{stop} '
           'output.experiments=experiments_{start}_{stop}.json '
           'output.reflections=indexed_{start}_{stop}.pickle '
           'output.log=None output.debug_log=None').format(**fmt_dic)
    result = easy_run.fully_buffered(command=cmd)

    # check if indexing worked
    new_exp_path = 'experiments_{start}_{stop}.json'.format(**fmt_dic)
    indexed_path = 'indexed_{start}_{stop}.pickle'.format(**fmt_dic)
    tst = [os.path.exists(e) for e in (new_exp_path, indexed_path)]
    if tst.count(True) != 2:
      # print the errors
      #for l in result.stderr_lines: print l
      return (None, None)
    return new_exp_path, indexed_path

  def __call__(self):

    first, last = self._scan.get_array_range()
    start = first

    while True:

      # keep trying to index, extending the block size looking for success
      nblocks=1
      while True:
        stop = min(start + nblocks * self._images_per_block, last)
        new_exp_path, indexed_path = self._index_current_block(start, stop)
        # exit if successful indexing
        if [new_exp_path, indexed_path].count(None) == 0: break
        # exit if no more data to include
        if stop == last:break
        nblocks += 1

      if stop == last:break

      # it worked, so update the pointer to the current model
      self._current_exp_path = new_exp_path

      # load the indexed results
      el = ExperimentListFactory.from_json_file(self._current_exp_path,
        check_format=False)
      rt = flex.reflection_table.from_pickle(indexed_path)
      indexed = rt.select(rt.get_flags(rt.flags.indexed))
      bc = [beam_centre_mm(e.detector, e.beam) for e in el]
      bc = [e.detector[panel_id].millimeter_to_pixel((x, y)) for e, (panel_id, (x, y)) in zip(el, bc)]
      indexed_rmsds = self._rmsds(indexed)
      cell = el[0].crystal.get_unit_cell().parameters()

      msg = self._current_exp_path + " " + "%.2f %.2f %.2f %.2f %.2f %.2f " % cell

      for (x, y) in bc:
        msg += "Beam centre (px): (%.2f,%.2f) " % (x, y)

      msg += "RMSD X: {:.4f} RMSD Y: {:.4f} RMSD Z: {:.4f}".format(*indexed_rmsds)
      print msg

      # update the scan range in the experiments and save
      for exp in el:
        exp.scan.swap(exp.scan[start:stop])

      dump = ExperimentListDumper(el)
      dump.as_json("sliced_" + self._current_exp_path)

      # update start image for the next block
      start = stop


  def _rmsds(self, reflections):
    """calculate unweighted RMSDs for the specified reflections"""

    x_calc, y_calc, z_calc = reflections['xyzcal.px'].parts()
    x_obs, y_obs, z_obs = reflections['xyzobs.px.value'].parts()

    x_resid2 = (x_calc - x_obs)**2
    y_resid2 = (y_calc - y_obs)**2
    z_resid2 = (z_calc - z_obs)**2

    resid_x = flex.sum(x_resid2)
    resid_y = flex.sum(y_resid2)
    resid_z = flex.sum(z_resid2)
    n = len(reflections)

    rmsds = (sqrt(resid_x / n),
              sqrt(resid_y / n),
              sqrt(resid_z / n))
    return rmsds

if __name__ == "__main__":

  run = Script()
  run()
