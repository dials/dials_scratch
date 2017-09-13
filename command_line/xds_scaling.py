#!/usr/bin/env python
# coding: utf-8

"""
xds_scaling.py performs an xds-like parameterisation of scaling and outputs the
calcualted inverse scale factors to a integrated_scaled.pickle file.
Unfortunately this currently runs quite slowly on large datasets such as the
thaumatin tutorial data, particularly the data reshaping before minimisation.

Usage:
  dials_scratch.xds_scaling integrated.pickle integrated_experiments.json
"""

from __future__ import absolute_import, division

import sys
import logging

import dials.util.log
from dials.util import halraiser
from dials.array_family import flex
from dials.util.options import OptionParser, flatten_reflections, flatten_experiments
import libtbx.load_env
from libtbx import phil

logger = logging.getLogger(libtbx.env.dispatcher_name)

#phil_scope = phil.parse('''
#  debug = False
#    .type = bool
#    .help = "Output additional debugging information"
#  output {
#    log = 'dials_scratch.xds_scaling.log'
#      .type = str
#      .help = "The log filename"
#
#    debug_log = 'dials_scratch.xds_scaling.debug.log'
#      .type = str
#      .help = "The debug log filename"
#  }
#''')

phil_scope = phil.parse('')

from dials_scratch.jbe.scaling_code import minimiser_functions as mf
from dials_scratch.jbe.scaling_code import data_manager_functions as dmf
from dials_scratch.jbe.scaling_code.data_quality_assessment import R_meas, R_pim

def scaling_lbfgs(reflections, experiments, gridding_parameters, scaling_options):
  """This algorithm performs an xds-like scaling"""
  '''handling of choice of integration method'''
  if scaling_options['integration_method'] == 'profile_fitting':
    int_method = (str('intensity.prf.value'), str('intensity.prf.variance'))
  elif scaling_options['integration_method'] == 'summation_integration':
      int_method = (str('intensity.sum.value'), str('intensity.sum.variance'))
  else:
      raise ValueError('Incorrect choice of integration_method')

  '''create a data manager object'''
  loaded_reflections = dmf.XDS_Data_Manager(reflections, experiments, int_method,
                                                 gridding_parameters)

  """Filter out zero/negative values of d, variance, etc"""
  loaded_reflections.filter_data(int_method[1], -10.0, 0.0)
  loaded_reflections.filter_data('d', -1.0, 0.0)

  '''Map the indices to the asu and also sort the reflection table by miller index'''
  loaded_reflections.map_indices_to_asu()

  '''determine gridding index for scale parameters '''
  loaded_reflections.initialise_scale_factors()

  '''assign a unique reflection index to each group of reflections'''
  loaded_reflections.assign_h_index()
  loaded_reflections.scale_by_LP_and_dqe()
    
  '''call the optimiser on the Data Manager object'''
  decay_correction_rescaling = scaling_options['decay_correction_rescaling']
  minimised = mf.LBFGS_optimiser(loaded_reflections, param_name='g_absorption',
                                 parameters=loaded_reflections.g_absorption)
  minimised = mf.LBFGS_optimiser(minimised.data_manager, param_name='g_decay',
                                 parameters=minimised.data_manager.g_decay,
                                 decay_correction_rescaling=decay_correction_rescaling)
  minimised = mf.LBFGS_optimiser(minimised.data_manager, param_name='g_modulation',
                                 parameters=minimised.data_manager.g_modulation)
  return minimised

def apply_scaling(reflections, experiments):
  #default parameters
  gridding_parameters = {'ndbins':20, 'nzbins':20, 'n_detector_bins':37}
  scaling_options = {'integration_method':'profile_fitting',
                     'decay_correction_rescaling':True}
  minimised = scaling_lbfgs(reflections, experiments, gridding_parameters, 
                            scaling_options)

  '''clean up reflection table for outputting'''
  minimised.data_manager.sorted_reflections['inverse_scale_factor'] = (
    minimised.data_manager.scale_factors)
  minimised.data_manager.initial_keys.append('inverse_scale_factor')
  for key in minimised.data_manager.reflection_table.keys():
    if not key in minimised.data_manager.initial_keys:
      del minimised.data_manager.sorted_reflections[key]

  '''save data'''
  filename = sys.argv[1].strip('.pickle')+str('_scaled.pickle')
  minimised.data_manager.save_sorted_reflections(filename)
  print "Saved output to " + str(filename)

  '''calculate R metrics'''
  Rmeas = R_meas(minimised.data_manager)
  Rpim = R_pim(minimised.data_manager)
  print "R_meas is %s" % (Rmeas)
  print "R_pim is %s" % (Rpim)


def main(argv):
  optionparser = OptionParser(
    usage=__doc__.strip(),
    read_experiments=True,
    read_reflections=True,
    read_datablocks=False,
    phil=phil_scope,
    check_format=False)
  params, options = optionparser.parse_args(argv)

  #dials.util.log.config(
  #  verbosity=options.verbose,
  #  info=params.output.log,
  #  debug=params.output.debug_log)

  if not params.input.experiments or not params.input.reflections:
    optionparser.print_help()
    return

  # UNWRAP all of the data objects from the PHIL parser
  reflections = flatten_reflections(params.input.reflections)
  experiments = flatten_experiments(params.input.experiments)
  
  apply_scaling(reflections, experiments)

if __name__ == "__main__":
  try:
    sys.exit(main(sys.argv[1:]))
  except Exception as e:
    halraiser(e)
