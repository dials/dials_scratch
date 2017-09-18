#!/usr/bin/env python
# coding: utf-8

"""
xds_scaling.py performs an xds-like parameterisation of scaling and outputs the
calcualted inverse scale factors to a integrated_scaled.pickle file.
Unfortunately this currently runs quite slowly on large datasets such as the
thaumatin tutorial data, particularly the data reshaping before minimisation.

Usage:
  dials_scratch.xds_scaling integrated.pickle integrated_experiments.json [options]

A number of options can be specified, see the phil_scope below.
"""

from __future__ import absolute_import, division
import libtbx.load_env
import logging
logger = logging.getLogger(libtbx.env.dispatcher_name)

import sys


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
  n_z_bins = 20
    .type = int
    .help = "Number of bins for phi/time gridding"
  n_detector_bins = 37
    .type = int
    .help = "Number of bins in each detector dimension for modulation gridding"
  integration_method = 'prf'
    .type = str
    .help = "Option to choose from profile fitted intensities (prf) or summation integrated intensities (sum)"
  decay_correction_rescaling = False
    .type = bool
    .help = "Option to turn on a relative-B factor rescale to the decay scale factors"
  modulation = True
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
plot_data_absorption, plot_data_modulation)


def main(argv):
  
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
  output_path =  output_path[0].rstrip('.pickle')+'_scaled.pickle'
  
  # UNWRAP all of the data objects from the PHIL parser
  reflections = flatten_reflections(params.input.reflections)
  experiments = flatten_experiments(params.input.experiments)
  
  apply_scaling(reflections, experiments, optionparser, output_path, logger=logger)

def scaling_lbfgs(reflections, experiments, scaling_options, logger):
  """This algorithm performs an xds-like scaling"""
  '''handling of choice of integration method'''

  if scaling_options['integration_method'] not in ['prf', 'sum']:
    print 'Invalid integration_method choice, using default profile fitted intensities'
    scaling_options['integration_method'] = 'prf'

  '''create a data manager object'''
  loaded_reflections = dmf.XDS_Data_Manager(reflections, experiments, scaling_options)

  """Filter out zero/negative values of d, variance, etc"""
  loaded_reflections.filter_negative_variances()
  loaded_reflections.filter_data('d', -1.0, 0.0)

  '''Map the indices to the asu and also sort the reflection table by miller index'''
  loaded_reflections.map_indices_to_asu()
  loaded_reflections.scale_by_LP_and_dqe()

  #set scaling weightings
  loaded_reflections.scale_weight_Isigma(3.0)
  loaded_reflections.scale_weight_dmin(1.4)

  '''assign a unique reflection index to each group of reflections'''
  '''determine gridding index for scale parameters '''
  loaded_reflections.assign_h_index()
  loaded_reflections.initialise_scale_factors()

  '''call the optimiser on the Data Manager object'''
  decay_correction_rescaling = scaling_options['decay_correction_rescaling']
  
  if scaling_options['absorption']:
    loaded_reflections = mf.LBFGS_optimiser(loaded_reflections,
                                            param_name='g_absorption',
                                            parameters=loaded_reflections.g_absorption
                                           ).return_data_manager()
  if scaling_options['decay']:
    loaded_reflections = mf.LBFGS_optimiser(loaded_reflections,
                                            param_name='g_decay',
                                            parameters=loaded_reflections.g_decay,
                                            decay_correction_rescaling=decay_correction_rescaling
                                           ).return_data_manager()
  if scaling_options['modulation']:
    loaded_reflections = mf.LBFGS_optimiser(loaded_reflections,
                                            param_name='g_modulation',
                                            parameters=loaded_reflections.g_modulation
                                           ).return_data_manager()
  loaded_reflections.sorted_reflections['inverse_scale_factor'] = loaded_reflections.scale_factors
  return loaded_reflections

def apply_scaling(reflections, experiments, optionparser, output_path, logger):
  '''extract parameters from phil'''
  phil_parameters = optionparser.phil
  diff_phil_parameters = optionparser.diff_phil

  logger.info("=" * 80)
  logger.info("")
  logger.info("Initialising")
  logger.info("")
  print "Initialising data structures...."

  scaling_options={'n_d_bins' : None, 'n_z_bins' : None, 'n_detector_bins' : None, 
                  'integration_method' : None, 'decay_correction_rescaling': False,
                  'modulation' : True, 'decay' : True, 'absorption' : True}
  for obj in phil_parameters.objects:
    if obj.name in scaling_options:
      scaling_options[obj.name] = obj.extract()
  for obj in diff_phil_parameters.objects:
    if obj.name in scaling_options:
      scaling_options[obj.name] = obj.extract()

  logger.info("Scaling options being used are :") 
  for k,v in scaling_options.iteritems():
    logger.info('%s : %s' % (k, v))

  '''do lbfgs minimisation'''
  minimised = scaling_lbfgs(reflections, experiments, scaling_options, logger=logger)

  '''calculate R metrics'''
  Rmeas = R_meas(minimised)
  Rpim = R_pim(minimised)
  print "R_meas is %s" % (Rmeas)
  print "R_pim is %s" % (Rpim)

  '''clean up reflection table for outputting and save data'''
  minimised.clean_reflection_table()
  minimised.save_sorted_reflections(output_path)
  print "Saved output to " + str(output_path)

  '''output plots of scale factors'''
  if scaling_options['absorption']:
    plot_data_absorption(minimised)
  if scaling_options['decay']:
    plot_data_decay(minimised)
  if scaling_options['modulation']:
    plot_data_modulation(minimised)
  print "Saved plots of correction factors"


if __name__ == "__main__":
  try:
    sys.exit(main(sys.argv[1:]))
  except Exception as e:
    halraiser(e)
