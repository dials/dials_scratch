#!/usr/bin/env python
# coding: utf-8

"""
Usage: dials_scratch.aimless_scaling integrated.pickle integrated_experiments.json
[integrated.pickle(2) integrated_experiments.json(2) ....] [options]

This program performs scaling on the input datasets, using a physical
parameterisation based on that used in the program Aimless. If multiple input
files have been specified, the datasets will be jointly scaled against a common
target of unique reflection intensities.

By default, a scale, decay and absorption correction parameterisation for each
dataset is used. An 'inverse_scale_factor' column is added to the integrated pickle
file for each dataset (output to an integrated_scaled pickle file), and if there
are multiple datasets a combined reflection table is also saved.

Scaling against a target can also be performed, by inputting one dataset
followed by the option target=integrated_scaled.pickle, where integrated_scaled.pickle
is a reflection table with an 'inverse_scale_factor' column from a previous
scaling run.
It is necessary to use the command scaling_method=KB in order to perform simple
KB scaling of a small-wedge dataset against a target dataset.
"""

from __future__ import absolute_import, division, print_function
import libtbx.load_env
import time
start_time=time.time()
import logging
logger = logging.getLogger('dials.scale')

import sys

from dials.util import halraiser
from dials.util.options import OptionParser, flatten_reflections, flatten_experiments
from libtbx import phil

phil_scope = phil.parse('''
  debug = False
    .type = bool
    .help = "Output additional debugging information"
  scaling_method = 'aimless'
      .type = str
      .help = "Set method for scaling - 'aimless' or 'KB for simple KB scaling"
  output {
    log = dials_scratch.scaling.log
      .type = str
      .help = "The log filename"
    debug_log = dials_scratch.scaling.debug.log
      .type = str
      .help = "The debug log filename"
    plot_scalefactors = True
      .type = bool
      .help = "Option to switch off scalefactor plotting."
  }
  include scope dials_scratch.jbe.scaling_code.scaling_options.phil_scope
''', process_includes=True)


from dials_scratch.jbe.scaling_code import minimiser_functions as mf
from dials_scratch.jbe.scaling_code import data_manager_functions as dmf
from dials_scratch.jbe.scaling_code.data_quality_assessment import R_pim_meas

def main(argv):
  '''main script to run the scaling algorithm'''

  optionparser = OptionParser(usage=__doc__.strip(), read_experiments=True,
    read_reflections=True, read_datablocks=False, phil=phil_scope,
    check_format=False)
  params, options = optionparser.parse_args(argv, show_diff_phil=False)

  import libtbx.load_env
  from libtbx.utils import Sorry
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

  # Unwrap all of the data objects from the PHIL parser
  reflections = flatten_reflections(params.input.reflections)
  experiments = flatten_experiments(params.input.experiments)

  target_reflections = None
  if params.scaling_options.target:
    target_params, options = optionparser.parse_args([params.scaling_options.target])
    target_reflections = flatten_reflections(target_params.input.reflections)[0]
    logger.info(('{0} has been selected as a target dataset for scaling. \n').format(
      params.scaling_options.target))
  # how to make sure correct experiments/reflections are matched for multiple datasets?

  # Check number of input files and determine whether multi mode or not.
  len_refl = len(reflections)
  len_exp = len(experiments)
  if len_refl > 1 and len_exp > 1 and (len_exp == len_refl):
    params.scaling_options.__inject__('multi_mode', True)
  elif len_refl == 1 and len_exp == 1:
    params.scaling_options.__inject__('multi_mode', False)
  else:
    assert 0, """Incorrect number of reflection and/or experiment files entered
    in the command line: must be an equal number of each"""

  if params.scaling_options.minimisation_parameterisation not in ['standard', 'log']:
    msg = ('Invalid minimisation parameterisation choice, proceeding using {sep}'
           'using standard g-value parameterisation').format(sep='\n')
    logger.info(msg)
    params.scaling_options.minimisation_parameterisation = 'standard'

  # do the main scaling
  if target_reflections:
    minimised = scale_against_target(reflections, experiments, target_reflections,
      params)
  else:
    #params.__inject__('scaling_method', 'aimless') #add this as a command line
    #parameter when XSCALE choice available?
    minimised = aimless_scaling_lbfgs(reflections, experiments, params)

  # calculate R metrics
  if not params.scaling_options.target:
    logger.info('Calculating metrics for scaling quality assessment.')
    '''calculate R metrics'''
    if params.scaling_options.multi_mode:
      Rpim, Rmeas = R_pim_meas(minimised)
      logger.info(("R_meas of the combined scaled dataset is {0:.6f}").format(Rmeas))
      logger.info(("R_pim of the combined scaled dataset is {0:.6f} {sep}").format(
        Rpim, sep='\n'))
      for j, dm in enumerate(minimised.data_managers):
        Rpim, Rmeas = R_pim_meas(dm)
        logger.info(("R_meas of the scaled dataset {0} is {1:.6f}").format(j+1, Rmeas))
        logger.info(("R_pim of the scaled dataset {0} is {1:.6f} {sep}").format(
          j+1, Rpim, sep='\n'))
    else:
      Rpim, Rmeas = R_pim_meas(minimised)
      logger.info(("R_meas of the scaled dataset is {0:.6f}").format(Rmeas))
      logger.info(("R_pim of the scaled dataset is {0:.6f} {sep}").format(Rpim, sep='\n'))

  # Plot scalefactors
  if params.output.plot_scalefactors and not params.scaling_options.target:
    from dials_scratch.jbe.scaling_code.data_plotter import (plot_smooth_scales,
      plot_absorption_surface)
    logger.info('\nPlotting graphs of scale factors. \n')
    if params.scaling_options.multi_mode:
      for j, dm in enumerate(minimised.data_managers):
        plot_smooth_scales(dm, outputfile='smooth_scale_factors_'+str(j+1)+'.png')
        if params.parameterisation.absorption_term:
          plot_absorption_surface(dm, outputfile='absorption_surface_'+str(j+1)+'.png')
    else:
      plot_smooth_scales(minimised, outputfile='smooth_scale_factors.png')
      if params.parameterisation.absorption_term:
        plot_absorption_surface(minimised)
    logger.info('Saved plots of correction factors. \n')

  # Clean up reflection table for outputting and save data
  if params.scaling_options.multi_mode:
    for j, dm in enumerate(minimised.data_managers):
      dm.clean_reflection_table()
      dm.save_reflection_table('integrated_scaled_'+str(j+1)+'.pickle')
      logger.info(('Saved output to {0}').format('integrated_scaled_'+str(j+1)+'.pickle'))
  elif params.scaling_options.target:
    minimised.dm1.save_reflection_table('integrated_targetscaled.pickle')
    logger.info("Saved output to " + str('integrated_targetscaled.pickle'))
  else:
    minimised.clean_reflection_table()
    minimised.save_reflection_table('integrated_scaled.pickle')
    logger.info(('\nSaved output to {0}').format('integrated_scaled.pickle'))

  # All done!
  finish_time=time.time()
  logger.info("\nTime taken: {0:.4f}s ".format((finish_time - start_time)))
  logger.info('\n'+'*'*40+'\n')

def scale_against_target(reflections, experiments, target_reflections, params):
  """This algorithm performs scaling against a target scaled reflection table"""
  logger.info('\n'+'*'*40+'\n')
  loaded_reflections = dmf.TargetedDataManager(reflections[0], experiments[0],
    target_reflections, params)

  '''call the optimiser on the Data Manager object'''
  param_name = []
  if params.parameterisation.scale_term:
    param_name.append('g_scale')
  if params.parameterisation.decay_term:
    param_name.append('g_decay')
  if not param_name:
      assert 0, 'no parameters have been chosen for scaling, aborting process'
  loaded_reflections = mf.LBFGS_optimiser(loaded_reflections,
    param_name=param_name).return_data_manager()

  '''the minimisation has only been done on a subset on the data, so apply the
  scale factors to the sorted reflection table.'''
  loaded_reflections.expand_scales_to_all_reflections()
  return loaded_reflections


def aimless_scaling_lbfgs(reflections, experiments, params):
  """This algorithm performs an aimless-like scaling"""
  logger.info('\n'+'*'*40+'\n')

  # Initialise the datamanager.
  if params.scaling_options.multi_mode:
    loaded_reflections = dmf.MultiCrystalDataManager(reflections, experiments,
      params)
  else:
    loaded_reflections = dmf.AimlessDataManager(reflections[0], experiments[0],
      params)

  # for now, assume you always want a scale, therefore option of doing decay also ////fix this?

  # Build a param_name list based on the scaling_options
  if params.scaling_options.concurrent_scaling:
    param_name = ['g_scale'] # force a scale term for now
    if params.parameterisation.decay_term:
      param_name.append('g_decay')
    if params.parameterisation.absorption_term:
      param_name.append('g_absorption')
    if not param_name:
      assert 0, 'no parameters have been chosen for scaling, aborting process'
    # Call the optimiser on the Data Manager object
    loaded_reflections = mf.LBFGS_optimiser(loaded_reflections,
      param_name=param_name).return_data_manager()
    # Optimise the error model and then do another minimisation
    if params.weighting.optimise_error_model:
      loaded_reflections.update_error_model()
      # Second minimisation with new weights
      loaded_reflections = mf.LBFGS_optimiser(loaded_reflections,
        param_name=param_name).return_data_manager()

  else: # Not concurrent_scaling, so do scale/decay term first then absorption
    if params.parameterisation.decay_term:
      loaded_reflections = mf.LBFGS_optimiser(loaded_reflections,
        param_name=['g_scale', 'g_decay']).return_data_manager()
    else: #just do scale factor if you don't want decay.
      loaded_reflections = mf.LBFGS_optimiser(loaded_reflections,
        param_name=['g_scale']).return_data_manager()
    if params.parameterisation.absorption_term:
      loaded_reflections = mf.LBFGS_optimiser(loaded_reflections,
        param_name=['g_absorption']).return_data_manager()
    # Optimise the error model and then do another minimisation
    if params.weighting.optimise_error_model:
      loaded_reflections.update_error_model()
      # Second pass
      if params.parameterisation.decay_term:
        loaded_reflections = mf.LBFGS_optimiser(loaded_reflections,
          param_name=['g_scale', 'g_decay']).return_data_manager()
      else: #just do scale factor if you don't want decay.
        loaded_reflections = mf.LBFGS_optimiser(loaded_reflections,
          param_name=['g_scale']).return_data_manager()
      if params.parameterisation.absorption_term:
        loaded_reflections = mf.LBFGS_optimiser(loaded_reflections,
          param_name=['g_absorption']).return_data_manager()

  # The minimisation has only been done on a subset on the data, so apply the
  # scale factors to the whole reflection table.
  loaded_reflections.expand_scales_to_all_reflections()
  if params.scaling_options.multi_mode:
    loaded_reflections.join_multiple_datasets()
  return loaded_reflections


if __name__ == "__main__":
  try:
    sys.exit(main(sys.argv[1:]))
  except Exception as e:
    halraiser(e)
