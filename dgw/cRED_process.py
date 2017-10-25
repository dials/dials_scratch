#!/usr/bin/env python
#
#  Copyright (C) 2017 Diamond Light Source and STFC Rutherford Appleton
#  Laboratory, UK.
#
#  Author: David Waterman
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

"""
Processing protocols for cRED lysozyme data. Pass in a directory (or
directories) of datasets, in which images consist of files with the extension
.img
"""

from __future__ import division
import os
import glob
from libtbx.utils import Sorry
from dials.test import cd
from libtbx import easy_run, Auto
from dials.array_family import flex
from scitbx import matrix
from dxtbx.model.experiment_list import ExperimentListFactory

class Script(object):
  '''A class for running the script.'''

  def __init__(self):
    '''Initialise the script.'''
    from dials.util.options import OptionParser
    from libtbx.phil import parse
    import libtbx.load_env

    # The phil scope
    phil_scope = parse('''
      data_directory = None
        .type = str
        .help = "Parent directory containing subdirectories of datasets"
        .multiple = True

      output_directory = "dials-proc"
        .type = str
        .help = "Name for a directory of output (must not exist already)"

      nproc = Auto
        .type = int
        .help = "Set nproc for spot finding and integration"

      find_spots
      {
        d_min = 1.9
          .type = float
          .help = "resolution limit for spot-finding"
      }

      integrate
      {
        d_min = 1.9
          .type = float
          .help = "resolution limit for integration"
      }''', process_includes=True)

    # The script usage
    import __main__
    usage  = ("usage: dials.python {0} [options] [param.phil] "
              "data_directory=/path/to/datadir".format(__main__.__file__))

    # Create the parser
    self.parser = OptionParser(
      usage=usage,
      phil=phil_scope,
      epilog=__doc__)

    return

  @staticmethod
  def is_dataset_dir(data_dir):
    if not os.path.isdir(data_dir): return False

    # check there are at least 5 images with the extension .img
    has_images = len(glob.glob(os.path.join(data_dir, '*.img'))) > 5

    return has_images

  def run(self):
    '''Execute the script.'''

    # Parse the command line
    self.params, options = self.parser.parse_args(show_diff_phil=True)

    if len(self.params.data_directory) == 0:
      raise Sorry('Please provide input in the form:\n  {0}'.format(self.parser.usage))

    # Get full paths to datasets
    self.params.data_directory = [os.path.realpath(e) for e in self.params.data_directory]

    # Check that these do contain datasets
    have_data = [self.is_dataset_dir(e) for e in self.params.data_directory]
    if have_data.count(False) > 0:
      no_data = [e for e, h in zip(
        self.params.data_directory, have_data) if not h]
      raise Sorry('No data found in found in:\n  ' + "\n  ".join(no_data))

    # Set nproc
    if self.params.nproc is Auto:
      from libtbx.introspection import number_of_processors
      self.params.nproc = number_of_processors(return_value_if_unknown=-1)

    # Look for CCP4
    if 'CCP4' not in os.environ:
      raise Sorry('Please set up the CCP4 environment')

    # Create processing directory and cd to it
    if os.path.isdir(self.params.output_directory):
      raise Sorry("Output directory already exists")
    try:
      os.mkdir(self.params.output_directory)
    except OSError:
      raise Sorry("Failed to create output directory")
    os.chdir(self.params.output_directory)

    for job in self.params.data_directory:

      job_dir = os.path.basename(job)

      # For each job, do processing in its own directory
      try:
        with cd(job_dir):
          self.process(job)
      except OSError:
        raise Sorry("Failed to create job directory {0}".format(job_dir))

    return

  def _run_one_job(self, command, keywords=None):
    print '> ' + command
    try:
      result = easy_run.fully_buffered(command=command,
            stdin_lines=keywords).raise_if_errors()
    except RuntimeError:
      print "Failed job"
      return None
    return result

  @staticmethod
  def _get_aimless_summary(logfile):

    summary = []
    with open(logfile, 'r') as f:
      for line in f:
        if line.startswith('<!--SUMMARY_BEGIN-->'): break
      for line in f:
        if line.startswith('    $$ <!--SUMMARY_END-->'): break
        summary.append(line)
    return summary


  def process(self, job):
    '''Perform processing tasks for one dataset'''

    print "Processing job {0}".format(job)

    # Import
    path = os.path.join(job, "*.img")
    cmd = ('dials.import {0} ').format(path)
    if self._run_one_job(cmd) is None: return

    # Spot finding
    cmd = 'dials.find_spots datablock.json filter.d_min={0} nproc={1}'.format(
      self.params.find_spots.d_min, self.params.nproc)
    if self._run_one_job(cmd) is None: return
    with open('dials.find_spots.log', 'r') as log:
      while True:
        line = log.readline()
        if line.startswith('Histogram'): break
      for i in range(14):
        print log.readline(),

    # Attempt indexing in P1
    n_indexed = 0
    best_indexing = None
    cmd = ('dials.index datablock.json strong.pickle '
           'output.experiments=P1_experiments.json '
           'output.reflections=P1_indexed.pickle output.log=P1_dials.index.log')
    if self._run_one_job(cmd) is not None:
      indexed = flex.reflection_table.from_pickle('P1_indexed.pickle')
      n_indexed = indexed.get_flags(indexed.flags.indexed).count(True)
      best_indexing = ('P1_experiments.json','P1_indexed.pickle')
      el = ExperimentListFactory.from_json_file('P1_experiments.json')
      cell = el[0].crystal.get_unit_cell().parameters()
      cell_s = ', '.join(['{:.3f}'.format(e) for e in cell])
      print "{0} reflections indexed in P1 with unit cell ({1})".format(
        n_indexed, cell_s)

    # Attempt indexing in P222
    cmd = ('dials.index datablock.json strong.pickle '
           'output.experiments=P222_experiments.json '
           'output.reflections=P222_indexed.pickle '
           'output.log=P222_dials.index.log space_group=P222')
    if self._run_one_job(cmd) is not None:
      indexed = flex.reflection_table.from_pickle('P222_indexed.pickle')
      n_indexed2 = indexed.get_flags(indexed.flags.indexed).count(True)
      if n_indexed2 > n_indexed:
        best_indexing = ('P222_experiments.json','P222_indexed.pickle')
        n_indexed = n_indexed2
      el = ExperimentListFactory.from_json_file('P222_experiments.json')
      cell = el[0].crystal.get_unit_cell().parameters()
      cell_s = ', '.join(['{:.3f}'.format(e) for e in cell])
      print "{0} reflections indexed in P222 with unit cell ({1})".format(
        n_indexed2, cell_s)

    if best_indexing is None: return
    print "Selecting {0} and {1} for further processing".format(
      *best_indexing)

    # Convert best indexing solution to P1 and P222 versions
    cmd = ('dials.reindex {0} {1} space_group={2} '
           'output.experiments=best_experiments_conv_to_{2}.json '
           'output.reflections=best_reflections_conv_to_{2}.pickle ')
    if self._run_one_job(cmd.format(best_indexing[0],
        best_indexing[1], 'P1')) is None: return
    if self._run_one_job(cmd.format(best_indexing[0],
        best_indexing[1], 'P222')) is None: return

    # Further processing in sub directories
    for sub_job in ('P1', 'P222'):
      try:
        with cd(sub_job):
          experiments = 'best_experiments_conv_to_{0}.json'.format(sub_job)
          reflections = 'best_reflections_conv_to_{0}.pickle'.format(sub_job)
          print "Processing best indexing solution converted to {0}".format(
              sub_job)
          self._refine_onwards(experiments, reflections)
      except OSError:
        raise Sorry("Failed to create sub job directory {0}".format(sub_job))

    # Compare summary tables from aimless.log files from each job
    import difflib
    logpaths = [os.path.join('P1', 'aimless.log'),
                os.path.join('P222', 'aimless.log')]
    if [os.path.exists(e) for e in logpaths].count(True) == 2:
      aimlessP1 = self._get_aimless_summary(logpaths[0])
      aimlessP222 = self._get_aimless_summary(logpaths[1])
      diff = difflib.HtmlDiff(wrapcolumn=80).make_file(
          aimlessP1, aimlessP222, logpaths[0], logpaths[1])
      with open('aimless-diff.html', 'w') as f:
        f.writelines(diff)

    print
    return

  @staticmethod
  def _get_detector_geometry(detector):
    # use panel 0
    panel0 = detector[0]
    fast = matrix.col(panel0.get_fast_axis())
    slow = matrix.col(panel0.get_slow_axis())
    normal = fast.cross(slow).normalize()
    origin = matrix.col(panel0.get_origin())
    distance = origin.dot(normal)
    return {'fast': fast, 'slow': slow, 'origin': origin, 'distance': distance}

  def _refine_onwards(self, experiments, reflections):
    '''Perform processing tasks in a sub-directory for files matching the
    base_name in the parent directory'''

    exp_path = os.path.join('..', experiments)
    ref_path = os.path.join('..', reflections)

    # Static refinement, two rounds to improve outlier rejection. Increase
    # max_iterations from the default
    cmd = ('dials.refine {0} {1} max_iterations=100 '
           'output.experiments=static_01_refined_experiments.json '
           'output.reflections=static_01_refined_reflections.pickle '
           'output.log=static_01_dials.refine.log').format(exp_path, ref_path)
    if self._run_one_job(cmd) is None: return
    cmd = ('dials.refine static_01_refined_experiments.json '
           'static_01_refined_reflections.pickle '
           'max_iterations=100 '
           'output.experiments=static_02_refined_experiments.json '
           'output.reflections=static_02_refined_reflections.pickle '
           'output.log=static_02_dials.refine.log')
    if self._run_one_job(cmd) is None: return

    # Print cell and detector geometry
    el = ExperimentListFactory.from_json_file('static_02_refined_experiments.json')
    cell = el[0].crystal.get_unit_cell().parameters()
    cell_s = ', '.join(['{:.3f}'.format(e) for e in cell])
    print "Refined unit cell: ({0})".format(cell_s)
    d = self._get_detector_geometry(el[0].detector)
    print "Detector origin:  ({0:.5f}, {1:.5f}, {2:.5f})".format(*d['origin'].elems)
    print "Detector fast ax: ({0:.5f}, {1:.5f}, {2:.5f})".format(*d['fast'].elems)
    print "Detector slow ax: ({0:.5f}, {1:.5f}, {2:.5f})".format(*d['slow'].elems)
    print "Detector distance: {0}".format(d['distance'])

    # Scan-varying refinement
    cmd = ('dials.refine static_02_refined_experiments.json '
           'static_02_refined_reflections.pickle '
           'scan_varying=true '
           'max_iterations=100 '
           'output.experiments=sv_refined_experiments.json '
           'output.reflections=sv_refined_reflections.pickle '
           'output.log=sv_dials.refine.log')
    if self._run_one_job(cmd) is None: return

    # Plot scan-varying crystal
    cmd = 'dials.plot_scan_varying_crystal sv_refined_experiments.json'
    self._run_one_job(cmd)

    # Print scan-varying cell and new detector geometry
    el = ExperimentListFactory.from_json_file('sv_refined_experiments.json')
    xl = el[0].crystal
    abc = flex.vec3_double()
    angles = flex.vec3_double()
    for n in range(xl.num_scan_points):
      a, b, c, alpha, beta, gamma = xl.get_unit_cell_at_scan_point(n).parameters()
      abc.append((a, b, c))
      angles.append((alpha, beta, gamma))
    a, b, c = abc.mean()
    alpha, beta, gamma = angles.mean()
    cell_s = ', '.join(['{:.3f}'.format(e) for e in [a, b, c, alpha, beta, gamma]])
    print "Average refined unit cell: ({0})".format(cell_s)
    d = self._get_detector_geometry(el[0].detector)
    print "Detector origin:  ({0:.5f}, {1:.5f}, {2:.5f})".format(*d['origin'].elems)
    print "Detector fast ax: ({0:.5f}, {1:.5f}, {2:.5f})".format(*d['fast'].elems)
    print "Detector slow ax: ({0:.5f}, {1:.5f}, {2:.5f})".format(*d['slow'].elems)
    print "Detector distance: {0}".format(d['distance'])

    # integrate
    cmd = ('dials.integrate sv_refined_experiments.json '
           'sv_refined_reflections.pickle nproc={0} '
           'prediction.d_min={1}').format(
              self.params.nproc, self.params.integrate.d_min)
    if self._run_one_job(cmd) is None: return

    # FIXME print out some statistics from integration

    # export
    cmd = ('dials.export integrated_experiments.json integrated.pickle '
           'mtz.ignore_panels=true')
    if self._run_one_job(cmd) is None: return

    # pointless
    cmd = 'pointless hklin integrated.mtz hklout sorted.mtz > pointless.log'
    if self._run_one_job(cmd) is None: return

    # aimless
    cmd = 'aimless hklin sorted.mtz hklout scaled.mtz > aimless.log'
    keywords = ['resolution 80.0 2.0']
    if self._run_one_job(cmd, keywords=keywords) is None: return


    print "DONE"

    return

if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
