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

    # TODO loop through files in filelist, do analysis. Skip existing analyses,
    # unless force_analysis is True
    for i, f in enumerate(self.filelist):
      self.process(f, i)

    return

  def process(self, filename, i):
    '''Perform processing tasks for one indexed.pickle'''

    status = {'skipped':True,
              'spectra':False,
              'refinement':False,
              'refined_spectra':False}

    directory = os.path.join(self._directory, 'run%05d' % i)
    try:
      os.mkdir(directory)
    except OSError:
      if not self.params.force_analysis:
        print 'Skipping ' + directory + ' that already exists'
        return status
    status['skipped'] = False

    cwd = os.path.abspath(os.curdir)
    os.chdir(directory)

    # locate centroid_analysis script
    dials_dir = libtbx.env.find_in_repositories('dials')
    ca = os.path.join(dials_dir, 'algorithms', 'refinement', 'analysis',
      'centroid_analysis.py')

    # run analysis
    cmd = "dials.python {0} {1}".format(ca, filename)
    result = easy_run.fully_buffered(command=cmd)

    if result.return_code != 0:
      return status
    status['spectra'] = True

    # TODO refinement and re-running analysis - in a subdirectory?

    # change back to original directory
    os.chdir(cwd)

    return status

if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
