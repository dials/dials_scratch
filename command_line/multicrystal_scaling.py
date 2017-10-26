#!/usr/bin/env python
# coding: utf-8

"""
xds_scaling.py performs an xds-like parameterisation of scaling and outputs the
calcualted inverse scale factors to a integrated_scaled.pickle file.
Unfortunately this currently runs quite slowly on large datasets such as the
thaumatin tutorial data, particularly the data reshaping before minimisation.

Usage:
  dials_scratch.multicrystal_scaling integrated.pickle integrated_experiments.json [options]

A number of options can be specified, see the phil_scope below.
"""

from __future__ import absolute_import, division
import libtbx.load_env
import logging
logger = logging.getLogger(libtbx.env.dispatcher_name)

import sys
import numpy as np
import copy as copy

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
  modulation = False
    .type = bool
    .help = "Option to turn off modulation correction"
  decay = True
    .type = bool
    .help = "Option to turn off decay correction"
  absorption = True
    .type = bool
    .help = "Option to turn off absorption correction"
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

  output_path = [s for s in argv if 'integrated.pickle' in s]
  output_path = output_path[0].rstrip('.pickle') + '_scaled.pickle'

  # UNWRAP all of the data objects from the PHIL parser
  reflections = flatten_reflections(params.input.reflections)
  experiments = flatten_experiments(params.input.experiments)

  phil_parameters = optionparser.phil
  diff_phil_parameters = optionparser.diff_phil

  logger.info("=" * 80)
  logger.info("")
  logger.info("Initialising")
  logger.info("")
  print "Initialising data structures...."

  scaling_options = {'n_d_bins' : None, 'rotation_interval' : None, 'n_detector_bins' : None,
                     'integration_method' : None, 'modulation' : True,
                     'decay' : True, 'absorption' : True, 'Isigma_min' : 3.0,
                     'd_min' : 0.0, 'decay_correction_rescaling': False,
                     'parameterization': 'standard'}
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
  minimised = xds_scaling_lbfgs(reflections, experiments, scaling_options, logger)

  '''calculate R metrics'''
  Rmeas = R_meas(minimised)
  Rpim = R_pim(minimised)
  print "R_meas is %s" % (Rmeas)
  print "R_pim is %s" % (Rpim)
  Rmeas = R_meas(minimised.dm1)
  Rpim = R_pim(minimised.dm1)
  print "R_meas is %s" % (Rmeas)
  print "R_pim is %s" % (Rpim)
  Rmeas = R_meas(minimised.dm2)
  Rpim = R_pim(minimised.dm2)
  print "R_meas is %s" % (Rmeas)
  print "R_pim is %s" % (Rpim)

  '''output plots of scale factors'''
  if scaling_options['absorption']:
    plot_data_absorption(minimised.dm1)
    n_time_pos = minimised.dm1.g_absorption.ntime_parameters
    plot_correction_at_multiple_detector_areas(minimised.dm1, [0, n_time_pos // 5,
      2 * n_time_pos // 5, 3 * n_time_pos // 5, 4 * n_time_pos // 5, n_time_pos - 1])
  if scaling_options['decay']:
    plot_data_decay(minimised.dm1)
  if scaling_options['modulation']:
    plot_data_modulation(minimised.dm1)
  print "Saved plots of correction factors"

  '''clean up reflection table for outputting and save data'''
  minimised.dm1.clean_reflection_table()
  minimised.dm1.save_sorted_reflections(output_path)
  print "Saved output to " + str(output_path)




def xds_scaling_lbfgs(reflections, experiments, scaling_options, logger):
  """This algorithm performs an xds-like scaling"""

  '''create a data manager object. Upon creation, negative variance & d-values
  are filtered and the indices are mapped to the asu and sorted. scale factors
  are initialised to unity'''
  #split the reflection table in half
  reflection_table_1 = reflections[0]
  reflection_table_1['z_value'] = reflection_table_1['xyzobs.px.value'].parts()[2]
  reflection_table_2 = copy.deepcopy(reflection_table_1)
  z_max = max(reflection_table_1['z_value'])
  sel = reflection_table_1['z_value'] <= z_max/2.0
  reflection_table_1 = reflection_table_1.select(sel)
  sel = reflection_table_2['z_value'] > z_max/2.0
  reflection_table_2 = reflection_table_2.select(sel)

  loaded_reflections = dmf.multicrystal_datamanager([reflection_table_1], experiments,
                                                    [reflection_table_2], experiments,
                                                    scaling_options)

  '''call the optimiser on the Data Manager object'''
  if scaling_options['absorption']:
    loaded_reflections = mf.LBFGS_optimiser(loaded_reflections,
                                            param_name='g_absorption'
                                           ).return_data_manager()
  if scaling_options['decay']:
    loaded_reflections = mf.LBFGS_optimiser(loaded_reflections,
                                            param_name='g_decay'
                                           ).return_data_manager()
  if scaling_options['modulation']:
    loaded_reflections = mf.LBFGS_optimiser(loaded_reflections,
                                            param_name='g_modulation'
                                           ).return_data_manager()

  '''the minimisation has only been done on a subset on the data, so apply the
  scale factors to the sorted reflection table.'''
  loaded_reflections.expand_scales_to_all_reflections()
  loaded_reflections.join_multiple_datasets()
  return loaded_reflections


if __name__ == "__main__":
  try:
    sys.exit(main(sys.argv[1:]))
  except Exception as e:
    halraiser(e)
