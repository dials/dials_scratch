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

from dxtbx.model.experiment.experiment_list import ExperimentListFactory
from dxtbx.model.experiment.experiment_list import ExperimentListDumper
from dials.util.options import flatten_reflections
from dials.array_family import flex
from dials.command_line.show import beam_centre_mm
from libtbx import phil
from libtbx import easy_run
from libtbx.table_utils import simple_table
import os
from math import sqrt

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
    self._num_strong = len(params.input.reflections[0].data)
    print "{0} strong spots read from {1}".format(self._num_strong, self._strong_path)

    self._images_per_block = params.images_per_block

    self._all_indexed = None
    return

  def _index_current_block(self, job_id, start, stop):

    fmt_dic = {'experiments':self._current_exp_path,
               'indexed':self._strong_path,
               'start':start,
               'stop':stop,
               'job_id':job_id}

    cmd = ('dials.index {experiments} {indexed} '
           'scan_range={start},{stop} '
           'output.experiments=latest_experiments.json '
           'output.reflections=indexed_{job_id:03d}.pickle '
           'output.log=None output.debug_log=None').format(**fmt_dic)
    result = easy_run.fully_buffered(command=cmd)

    # check if indexing worked
    new_exp_path = 'latest_experiments.json'
    indexed_path = 'indexed_{job_id:03d}.pickle'.format(**fmt_dic)
    tst = [os.path.exists(e) for e in (new_exp_path, indexed_path)]
    if tst.count(True) != 2:
      return (None, None)
    return new_exp_path, indexed_path

  def __call__(self):

    # set up variables we need to determine blocks
    first, last = self._scan.get_array_range()
    start = first
    stop = first

    # set up table
    header = ["Job", "Scan\nrange", "#Idx", "Cell", "Beam centre\n(fast,slow)",
      "Distance\n(mm)", "RMSD_X\n(px)", "RMSD_Y\n(px)", "RMSD_Z\n(px)"]
    rows = []
    num_indexed = []
    filelist_lines = []
    job_id = 1

    while True:

      # finish if we already processed all blocks
      if stop == last:break

      # keep trying to index, extending the block size looking for success
      nblocks=1
      while True:
        stop = start + nblocks * self._images_per_block
        # if within one block size from the end, include all images up to the end
        if stop >= last - self._images_per_block: stop = last
        new_exp_path, indexed_path = self._index_current_block(job_id, start, stop)
        # exit if successful indexing
        if [new_exp_path, indexed_path].count(None) == 0: break
        # exit if no more data to include
        if stop == last:break
        nblocks += 1

      # exit if the last job failed
      if new_exp_path is None: break

      # it worked, so update the pointer to the current model
      self._current_exp_path = new_exp_path
      print "Job {0} completed: scan range ({1},{2})".format(job_id, start, stop)

      # load the indexed results
      el = ExperimentListFactory.from_json_file(self._current_exp_path,
        check_format=False)
      assert len(el) == 1
      exp = el[0]

      rt = flex.reflection_table.from_pickle(indexed_path)
      if self._all_indexed is None:
        self._all_indexed = rt
      else:
        self._all_indexed.extend(rt)

      indexed = rt.select(rt.get_flags(rt.flags.indexed))
      num_indexed.append(len(indexed))
      panel_id, (x,y) = beam_centre_mm(exp.detector, exp.beam)
      bc = exp.detector[panel_id].millimeter_to_pixel((x, y))
      indexed_rmsds = self._rmsds(indexed)
      cell = exp.crystal.get_unit_cell().parameters()
      pnl_dists = [p.get_distance() for p in exp.detector]
      rows.append(["{0}".format(job_id),
                   "({0},{1})".format(start, stop),
                   "%i" % num_indexed[-1],
                   "%.2f %.2f %.2f %.2f %.2f %.2f" % cell,
                   "%.2f %.2f" % bc,
                   (" ".join(["%.2f" % d for d in pnl_dists])),
                   "%.4f" % indexed_rmsds[0],
                   "%.4f" % indexed_rmsds[1],
                   "%.4f" % indexed_rmsds[2]])

      # update the scan range in the experiments and save
      for exp in el:
        exp.scan.swap(exp.scan[start:stop])

      dump = ExperimentListDumper(el)
      dump.as_json("experiments_{0:03d}.json".format(job_id))

      filelist_lines.append("block{0:03d} experiments_{0:03d}.json ".format(job_id) + indexed_path)

      # set up for the next block
      start = stop
      job_id += 1

    # print results
    st = simple_table(rows, header)
    print st.format()
    total_num_indexed = sum(num_indexed)
    print "{0} indexed reflections out of {1} strong spots ({2:.1f}%)".format(
      total_num_indexed, self._num_strong, 100.*total_num_indexed/self._num_strong)

    # write out the list of files
    print "writing out list of experiments and indexed reflections to filelist.txt"
    with open("filelist.txt", "w") as f:
      f.write("\n".join(filelist_lines))
      f.write("\n")

    # write out the full reflection table of indexed reflections
    print "writing out all indexed reflections to all_indexed.pickle"
    self._all_indexed.as_pickle('all_indexed.pickle')

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
