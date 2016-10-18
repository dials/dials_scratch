#!/usr/bin/env python
#
#  Copyright (C) 2016 Diamond Light Source and STFC Rutherford Appleton
#  Laboratory, UK.
#
#  Author: David Waterman
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

"""Work through a list of indexed.pickles and perform centroid analysis on
each to produce periodograms of the residuals. Save these plots and also do
scan-varying refinement, followed by a second round of centroid analysis to see
if the varying model accounts for the main features of the residuals.

This will be run on datasets from the Metrix database. The idea is to explore
datasets with many different features, to investigate in which situations it is
possible to automatically determine a sensible interval width for scan-varying
refinement.

Usage: dials.python analyse_periodograms.py filelist

where filelist consists of paths to indexed.pickles. Associated
experiments.jsons will be determined by matching expected filenames."""

from __future__ import division
import os
import datetime
from libtbx.utils import Sorry
from libtbx import easy_run
import libtbx.load_env
from dials.command_line.analyse_output import ensure_directory

class Script(object):
  '''A class for running the script.'''

  def __init__(self):
    '''Initialise the script.'''
    from dials.util.options import OptionParser
    from libtbx.phil import parse
    import libtbx.load_env

    # The phil scope
    phil_scope = parse('''
      force_analysis = False
        .type = bool
        .help = "Set True to force overwriting of previous analysis results."
                "Otherwise subdirectories for jobs that already exist will be"
                "skipped."

      output {
        directory = analyse_periodograms
          .type = path
          .help = "The root directory in which to store results"
      }
    ''', process_includes=True)

    # The script usage
    import __main__
    usage  = "usage: dials.python {0} [options] [param.phil] filelist ".format(
      __main__.__file__)

    # Create the parser
    self.parser = OptionParser(
      usage=usage,
      phil=phil_scope,
      epilog=__doc__)

    return

  def run(self):
    '''Execute the script.'''

    # Parse the command line
    self.params, options, filelist = self.parser.parse_args(show_diff_phil=True,
      return_unhandled=True)

    if len(filelist) != 1:
      raise Sorry('Please provide input in the form:\n  {0}'.format(self.parser.usage))

    # Read the filelist
    try:
      with open(filelist[0] ,'r') as f:
        self.filelist = f.readlines()
    except IOError:
      raise Sorry('Cannot read the specified file list')

    # Create output directory if it does not already exist
    self._directory = os.path.abspath(self.params.output.directory)
    ensure_directory(self._directory)
    print "Analysis will be performed in {0}".format(self._directory)

    # Do analysis for the files in filelist
    for i, f in enumerate(self.filelist):
      s = '{:%Y-%m-%d %H:%M:%S} '.format(datetime.datetime.now())
      s += 'Processing file {0}: '.format(i) + f.rstrip()
      print s
      self.process(f, i)

    return

  def process(self, idx_path, i):
    '''Perform processing tasks for one indexed.pickle'''

    cwd = os.path.abspath(os.curdir)
    status = self._process_core(idx_path, i)
    os.chdir(cwd)
    return status

  def _process_core(self, idx_path, i):

    status = {'skipped':True,
              'spectra':False,
              'refined':False,
              'refined_spectra':False}

    directory = os.path.join(self._directory, 'run%05d' % i)
    try:
      os.mkdir(directory)
    except OSError:
      if not self.params.force_analysis:
        print 'Skipping ' + directory + ' that already exists'
        return status
    status['skipped'] = False
    os.chdir(directory)

    # locate centroid_analysis script
    dials_dir = libtbx.env.find_in_repositories('dials')
    ca = os.path.join(dials_dir, 'algorithms', 'refinement', 'analysis',
      'centroid_analysis.py')

    # run analysis
    cmd = "dials.python {0} {1}".format(ca, idx_path)
    result = easy_run.fully_buffered(command=cmd)
    with open('cmd.txt', 'w') as f:
      f.write(result.command)

    if result.return_code != 0:
      return status
    status['spectra'] = True

    # try to find an associated experiment.json
    exp_path = self.find_experiments_json(idx_path)
    if exp_path is None:
      return status

    # refine, static then scan-varying
    try:
      os.mkdir('refined')
    except OSError:
      if not self.params.force_analysis:
        return status
    os.chdir('refined')
    cmd = 'dials.refine {0} {1}'.format(exp_path, idx_path)
    result = easy_run.fully_buffered(command=cmd)
    with open('cmd.txt', 'w') as f:
      f.write(result.command)
    tst = [os.path.exists(p) for p in ['refined.pickle',
                                       'refined_experiments.json']]
    if tst.count(True) != 2:
      return status
    cmd = 'dials.refine refined.pickle refined_experiments.json scan_varying=True'
    result = easy_run.fully_buffered(command=cmd)
    with open('cmd.txt', 'a') as f:
      f.write(result.command)
    if result.return_code != 0:
      return status
    status['refined'] = True

    # second round of analysis
    cmd = "dials.python {0} {1}".format(ca, 'refined.pickle')
    result = easy_run.fully_buffered(command=cmd)
    with open('cmd.txt', 'a') as f:
      f.write(result.command)

    if result.return_code != 0:
      return status
    status['refined_spectra'] = True

    return status

  def find_experiments_json(self, idx_path):
    '''Given the path to an indexed.pickle, try to identify an associated
    experiments.json in the same directory assuming xia2's naming convention'''

    head, tail = os.path.split(idx_path)
    root, ext = os.path.splitext(tail)
    exp_path = os.path.join(head,
      root.replace('indexed', 'experiments') + '.json')
    if os.path.exists(exp_path):
      return exp_path
    else:
      return None

if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
