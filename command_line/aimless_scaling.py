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
It is necessary to use the command scaling_model=KB in order to perform simple
KB scaling of a small-wedge dataset against a target dataset.
"""
from __future__ import absolute_import, division, print_function
#import libtbx.load_env
import time
import logging
import sys
from dials.util import halraiser
from dials.util.options import OptionParser, flatten_reflections, flatten_experiments
from libtbx import phil
from dials_scratch.jbe.scaling_code import minimiser_functions as mf
from dials_scratch.jbe.scaling_code import ScalingModelFactory
from dials_scratch.jbe.scaling_code import ScalerFactory

start_time = time.time()
logger = logging.getLogger('dials.scale')

phil_scope = phil.parse('''
  debug = False
    .type = bool
    .help = "Output additional debugging information"
  scaling_model = 'aimless'
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
    plot_merging_stats = False
      .type = bool
      .help = "Option to switch on plotting of merging stats."
    experiments_out = "scaled_experiments.json"
      .type = str
      .help = "Option to set filepath for output json."
    scaled_out = "scaled.pickle"
      .type = str
      .help = "Option to set filepath for output pickle file of scaled intensities."
  }
  include scope dials_scratch.jbe.scaling_code.scaling_options.phil_scope
''', process_includes=True)

def main(argv):
  '''main script to run the scaling algorithm'''

  optionparser = OptionParser(usage=__doc__.strip(), read_experiments=True,
    read_reflections=True, read_datablocks=False, phil=phil_scope,
    check_format=False)
  params, _ = optionparser.parse_args(argv, show_diff_phil=False)

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

  #put this bit somewhere else?
  if params.scaling_options.minimisation_parameterisation not in ['standard', 'log']:
    msg = ('Invalid minimisation parameterisation choice, proceeding using {sep}'
           'using standard g-value parameterisation').format(sep='\n')
    logger.info(msg)
    params.scaling_options.minimisation_parameterisation = 'standard'

  #first create the scaling model
  experiments = ScalingModelFactory.Factory.create(params, experiments, reflections)

  #now create the scaler
  scaler = ScalerFactory.Factory.create(params, experiments, reflections)

  #do the scaling
  minimised = lbfgs_scaling(scaler)

  #calculate merging stats
  results = minimised.calc_merging_statistics()
  logger.info('*'*40)
  logger.info("Dataset statistics")
  plot_labels = []
  #result.overall.show_summary(out=log.info_handle(logger))
  for i, result in enumerate(results):
    if len(results) == 1:
      logger.info("")
      result.overall.show_summary()
      plot_labels.append('Single dataset ')
    else:
      if i < len(results) - 1:
        logger.info("\nStatistics for dataset " + str(i+1))
        result.overall.show_summary()
        plot_labels.append('Dataset ' + str(i+1))
      else:
        logger.info("\nStatistics for combined datasets")
        result.overall.show_summary()
        plot_labels.append('Combined datasets')
  if params.output.plot_merging_stats:
    from xia2.command_line.compare_merging_stats import plot_merging_stats
    plot_merging_stats(results, labels=plot_labels)


  #save scaled_experiments.json file
  save_experiments(experiments, params.output.experiments_out)

  # Clean up reflection table for outputting and save data
  minimised.clean_reflection_table()
  minimised.save_reflection_table(params.output.scaled_out)
  logger.info(('\nSaved reflection table to {0}').format(params.output.scaled_out))

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

  # All done!
  finish_time = time.time()
  logger.info("\nTotal time taken: {0:.4f}s ".format((finish_time - start_time)))
  logger.info('\n'+'*'*40+'\n')

def lbfgs_scaling(scaler):
  '''main function to do lbfbs scaling'''
  logger.info('\n'+'*'*40+'\n')

  scaler.params.scaling_options.__inject__('multi_mode', False)

  if isinstance(scaler, ScalerFactory.TargetScaler):
    #do a prescaling round against a target
    param_name = []
    if scaler.params.parameterisation.scale_term:
      param_name.append('g_scale')
    if scaler.params.parameterisation.decay_term:
      param_name.append('g_decay')
    if scaler.params.parameterisation.absorption_term:
      param_name.append('g_absorption')
    if not param_name:
        assert 0, 'no parameters have been chosen for scaling, aborting process'
    scaler = mf.LBFGS_optimiser(scaler,
      param_name=param_name).return_data_manager()

    '''the minimisation has only been done on a subset on the data, so apply the
    scale factors to the sorted reflection table.'''
    scaler.expand_scales_to_all_reflections()
    #now pass to a multiscaler
    scaler = ScalerFactory.MultiScalerFactory.create_from_targetscaler(scaler)

  #from here onwards, scaler should only be a SingleScaler
  #or MultiScaler (not TargetScaler)
  if isinstance(scaler, ScalerFactory.MultiScaler):
    scaler.params.scaling_options.multi_mode = True

  if scaler.params.scaling_options.concurrent_scaling:
    param_name = ['g_scale'] # force a scale term for now
    if scaler.params.parameterisation.decay_term:
      param_name.append('g_decay')
    if scaler.params.parameterisation.absorption_term:
      param_name.append('g_absorption')
    if not param_name:
      assert 0, 'no parameters have been chosen for scaling, aborting process'
    # Call the optimiser on the Data Manager object
    scaler = mf.LBFGS_optimiser(scaler,
      param_name=param_name).return_data_manager()
    # Optimise the error model and then do another minimisation
    if scaler.params.weighting.optimise_error_model:
      scaler.update_error_model()
      # Second minimisation with new weights
      scaler = mf.LBFGS_optimiser(scaler,
        param_name=param_name).return_data_manager()

  else: # Not concurrent_scaling, so do scale/decay term first then absorption
    if scaler.params.parameterisation.decay_term:
      scaler = mf.LBFGS_optimiser(scaler,
        param_name=['g_scale', 'g_decay']).return_data_manager()
    else: #just do scale factor if you don't want decay.
      scaler = mf.LBFGS_optimiser(scaler,
        param_name=['g_scale']).return_data_manager()
    if scaler.params.parameterisation.absorption_term:
      scaler = mf.LBFGS_optimiser(scaler,
        param_name=['g_absorption']).return_data_manager()
    # Optimise the error model and then do another minimisation
    if scaler.params.weighting.optimise_error_model:
      scaler.update_error_model()
      # Second pass
      if scaler.params.parameterisation.decay_term:
        scaler = mf.LBFGS_optimiser(scaler,
          param_name=['g_scale', 'g_decay']).return_data_manager()
      else: #just do scale factor if you don't want decay.
        scaler = mf.LBFGS_optimiser(scaler,
          param_name=['g_scale']).return_data_manager()
      if scaler.params.parameterisation.absorption_term:
        scaler = mf.LBFGS_optimiser(scaler,
          param_name=['g_absorption']).return_data_manager()

  # The minimisation has only been done on a subset on the data, so apply the
  # scale factors to the whole reflection table.
  scaler.expand_scales_to_all_reflections()
  if isinstance(scaler, ScalerFactory.MultiScaler):
    scaler.join_multiple_datasets()
  return scaler

def save_experiments(experiments, filename):
  ''' Save the profile model parameters. '''
  from time import time
  from dxtbx.model.experiment_list import ExperimentListDumper
  st = time()
  logger.info('\nSaving the experiments to %s' % filename)
  dump = ExperimentListDumper(experiments)
  with open(filename, "w") as outfile:
    outfile.write(dump.as_json(split=True))
  logger.info('Time taken: %g' % (time() - st))


if __name__ == "__main__":
  try:
    sys.exit(main(sys.argv[1:]))
  except Exception as e:
    halraiser(e)
