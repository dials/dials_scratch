#!/usr/bin/env python
# coding: utf-8

"""
xds_scaling.py performs an xds-like parameterisation of scaling and outputs the
calculated inverse scale factors to a integrated_scaled.pickle file(s).

Usage:
  dials_scratch.xds_scaling integrated.pickle integrated_experiments.json
  target_integrated_scaled.pickle [options]

A number of options can be specified, see the phil_scope below.
"""

from __future__ import absolute_import, division
import libtbx.load_env
import logging
logger = logging.getLogger(libtbx.env.dispatcher_name)

import sys
import numpy as np


from dials.util import halraiser
from dials.array_family import flex
from dials.util.options import OptionParser, flatten_reflections, flatten_experiments
from libtbx import phil

phil_scope = phil.parse('''
  debug = False
    .type = bool
    .help = "Output additional debugging information"
  output {
    log = 'dials_scratch.xds_scaling.log'
      .type = str
      .help = "The log filename"

    debug_log = 'dials_scratch.xds_scaling.debug.log'
      .type = str
      .help = "The debug log filename"
  }
  verbosity = 1
    .type = int(value_min=0)
    .help = "The verbosity level"

  n_d_bins = 20
    .type = int
    .help = "Number of bins for resolution gridding"
  rotation_interval = None
    .type = float
    .help = "User specified rotation (phi) interval in degrees for phi binning"
  n_detector_bins = 19
    .type = int
    .help = "Number of bins in each detector dimension for modulation gridding"
  integration_method = 'prf'
    .type = str
    .help = "Option to choose from profile fitted intensities (prf)
             or summation integrated intensities (sum)"
  parameterization = 'standard'
    .type = str
    .help = "Choice of g-value parameterisation - 'standard' (multiplicative) or 'log'"
  decay_correction_rescaling = False
    .type = bool
    .help = "Option to turn on a relative-B factor rescale to the decay scale factors"
  Isigma_min = -5.0
    .type = float
    .help = "Option to use a I/sigma subset of reflections to determine scale factors"
  d_min = 0.0
    .type = float
    .help = "Option to use a d-value subset of reflections to determine scale factors"
''')

from dials_scratch.jbe.scaling_code import minimiser_functions as mf
from dials_scratch.jbe.scaling_code import data_manager_functions as dmf
from dials_scratch.jbe.scaling_code.data_quality_assessment import R_meas, R_pim
from dials_scratch.jbe.scaling_code.data_plotter import (plot_data_decay,
plot_data_absorption, plot_data_modulation, plot_correction_at_multiple_detector_areas)


def main(argv):
  '''main script to run and xds-like scaling algorithm'''
  optionparser = OptionParser(
    usage=__doc__.strip(),
    read_experiments=True,
    read_reflections=True,
    read_datablocks=False,
    phil=phil_scope,
    check_format=False)
  params, options = optionparser.parse_args(argv, show_diff_phil=False)

  from dials.util import log
  log.config(verbosity=1, info=params.output.log, debug=params.output.debug_log)

  if not params.input.experiments or not params.input.reflections:
    optionparser.print_help()
    return

  from dials.util.version import dials_version
  logger.info(dials_version())

  # Log the diff phil
  diff_phil = optionparser.diff_phil.as_str()
  if diff_phil is not '':
    logger.info('The following parameters have been modified:\n')
    logger.info(diff_phil)


  # UNWRAP all of the data objects from the PHIL parser
  reflections = flatten_reflections(params.input.reflections)
  experiments = flatten_experiments(params.input.experiments)

  scaling_options = {'n_d_bins' : None, 'rotation_interval' : None, 'n_detector_bins' : None,
                     'integration_method' : None, 'modulation' : True,
                     'decay' : True, 'absorption' : True, 'Isigma_min' : 3.0,
                     'd_min' : 0.0, 'decay_correction_rescaling': False,
                     'parameterization': 'standard', 'scaling_method' : 'KB'}

  if len(reflections) != 2:
    assert 0, """Incorrect number of reflection files entered
    in the command line (must be 2, one to scale and one integrated_scaled.pickle)"""
  elif len(experiments) != 1:
    assert 0, """Incorrect number of experiments files entered, only need one"""

  phil_parameters = optionparser.phil
  diff_phil_parameters = optionparser.diff_phil

  logger.info("=" * 80)
  logger.info("")
  logger.info("Initialising")
  logger.info("")
  print "Initialising data structures...."

  
  for obj in phil_parameters.objects:
    if obj.name in scaling_options:
      scaling_options[obj.name] = obj.extract()
  for obj in diff_phil_parameters.objects:
    if obj.name in scaling_options:
      scaling_options[obj.name] = obj.extract()
  '''handling of choice of integration method'''
  if scaling_options['integration_method'] not in ['prf', 'sum', 'combine']:
    print 'Invalid integration_method choice, using default profile fitted intensities'
    scaling_options['integration_method'] = 'prf'
  if scaling_options['parameterization'] not in ['standard', 'log']:
    print 'Invalid parameterization choice, using standard g-value parameterisation'
    scaling_options['integration_method'] = 'standard'

  logger.info("Scaling options being used are :")
  for k, v in scaling_options.iteritems():
    logger.info('%s : %s' % (k, v))

  '''do lbfgs minimisation'''
  minimised = scaling_lbfgs(reflections, experiments, scaling_options, logger)

  """'''output plots of scale factors'''
  if scaling_options['multi_mode']:
    if scaling_options['absorption']:
      plot_data_absorption(minimised.dm1, outputfile='g_absorption_multiset1.png')
      n_time_pos = minimised.dm1.g_absorption.ntime_parameters
      plot_correction_at_multiple_detector_areas(minimised.dm1, [0, n_time_pos // 5,
        2 * n_time_pos // 5, 3 * n_time_pos // 5, 4 * n_time_pos // 5, n_time_pos - 2],
        outputfile='g_absorption_surfaces_multiset1.png')
      plot_data_absorption(minimised.dm2, outputfile='g_absorption_multiset2.png')
      n_time_pos = minimised.dm2.g_absorption.ntime_parameters
      plot_correction_at_multiple_detector_areas(minimised.dm2, [0, n_time_pos // 5,
        2 * n_time_pos // 5, 3 * n_time_pos // 5, 4 * n_time_pos // 5, n_time_pos - 2],
        outputfile='g_absorption_surfaces_multiset2.png')
    if scaling_options['decay']:
      plot_data_decay(minimised.dm1, outputfile='g_decay_multiset1.png')
      plot_data_decay(minimised.dm2, outputfile='g_decay_multiset2.png')
    if scaling_options['modulation']:
      plot_data_modulation(minimised.dm1, outputfile='g_modulation_multiset1.png')
      plot_data_modulation(minimised.dm2, outputfile='g_modulation_multiset2.png')
    print "Saved plots of correction factors"
  else:
    if scaling_options['absorption']:
      plot_data_absorption(minimised)
      n_time_pos = minimised.g_absorption.ntime_parameters
      plot_correction_at_multiple_detector_areas(minimised, [0, n_time_pos // 5,
        2 * n_time_pos // 5, 3 * n_time_pos // 5, 4 * n_time_pos // 5, n_time_pos - 2])
    if scaling_options['decay']:
      plot_data_decay(minimised)
    if scaling_options['modulation']:
      plot_data_modulation(minimised)
    print Saved plots of correction factors"""

  '''clean up reflection table for outputting and save data'''
  #minimised.dm1.clean_reflection_table()
  minimised.dm1.save_reflection_table('integrated_targetscaled.pickle')
  print "Saved output to " + str('integrated_targetscaled.pickle')


def scaling_lbfgs(reflections, experiments, scaling_options, logger):
  """This algorithm performs scaling against a target scaled reflection table"""

  loaded_reflections = dmf.targeted_datamanager(reflections[0], 
    experiments[0], reflections[1], scaling_options)

  '''call the optimiser on the Data Manager object'''
  loaded_reflections = mf.LBFGS_optimiser(loaded_reflections,
                                          param_name=None
                                          ).return_data_manager()

  '''the minimisation has only been done on a subset on the data, so apply the
  scale factors to the sorted reflection table.'''
  loaded_reflections.expand_scales_to_all_reflections()
  return loaded_reflections


if __name__ == "__main__":
  try:
    sys.exit(main(sys.argv[1:]))
  except Exception as e:
    halraiser(e)
