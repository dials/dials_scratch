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
from dials_scratch.jbe.scaling_code import ParameterHandler

start_time = time.time()
logger = logging.getLogger('dials.scale')

phil_scope = phil.parse('''
  debug = False
    .type = bool
    .help = "Output additional debugging information"
  scaling_model = aimless
      .type = choice
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
  include scope dials_scratch.jbe.scaling_code.ScalingModelFactory.phil_scope
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

  diff_phil = optionparser.diff_phil.as_str()
  if diff_phil is not '':
    logger.info('The following parameters have been modified:\n')
    logger.info(diff_phil)

  reflections = flatten_reflections(params.input.reflections)
  experiments = flatten_experiments(params.input.experiments)

  def parse_multiple_datasets(reflections):
    'method to parse multiple datasets to see if any have already been scaled.'
    single_reflection_tables = []
    for refl_table in reflections:
      if 'dataset_id' in refl_table.keys():
        dataset_ids = set(refl_table['dataset_id'])
        if len(dataset_ids) > 1: #more than one scaled dataset in refl_tab
          logger.info(('\nDetected existence of a multi-dataset scaled reflection table {sep}'
            'containing {0} datasets. {sep}').format(len(dataset_ids), sep='\n'))
          for dataset_id in dataset_ids:
            single_refl_table = refl_table.select(refl_table['dataset_id'] == dataset_id)
            single_reflection_tables.append(single_refl_table)
        else: #only one already scaled dataset in refl_tab
          single_reflection_tables.append(refl_table)
      else: #refl_table has not previously been scaled
        single_reflection_tables.append(refl_table)
    logger.info("Found %s reflection tables in total." % len(single_reflection_tables))
    return single_reflection_tables

  if len(experiments) != 1:
    logger.info(('Checking for the existence of a reflection table containing {sep}'
    'multiple scaled datasets. {sep}').format(sep='\n'))
    reflections = parse_multiple_datasets(reflections)
    logger.info("Found %s experiments in total." % len(experiments))
    logger.info('\n'+'*'*40)

  assert len(experiments) == len(reflections), ''' mismatched number of
  experiments and reflection tables found'''

  '''first create the scaling model'''
  experiments = ScalingModelFactory.Factory.create(params, experiments, reflections)

  '''now create the scaler and do the scaling'''
  scaler = ScalerFactory.Factory.create(params, experiments, reflections)
  minimised = lbfgs_scaling(scaler)

  '''calculate merging stats'''
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
  
  '''plot merging stats if requested'''
  if params.output.plot_merging_stats:
    from xia2.command_line.compare_merging_stats import plot_merging_stats
    plot_merging_stats(results, labels=plot_labels)

  '''save scaled_experiments.json file'''
  save_experiments(experiments, params.output.experiments_out)

  '''plot scalefactors if requested'''
  if params.output.plot_scalefactors:# and not params.scaling_options.target:
    from dials_scratch.jbe.scaling_code.data_plotter import (plot_smooth_scales,
      plot_absorption_surface)
    logger.info('\nPlotting graphs of scale factors. \n')
    if isinstance(scaler, ScalerFactory.MultiScaler):
      for j, dm in enumerate(minimised.data_managers):
        plot_smooth_scales(dm, outputfile='smooth_scale_factors_'+str(j+1)+'.png')
        if params.parameterisation.absorption_term:
          plot_absorption_surface(dm, outputfile='absorption_surface_'+str(j+1)+'.png')
    else:
      plot_smooth_scales(minimised, outputfile='smooth_scale_factors.png')
      if params.parameterisation.absorption_term:
        plot_absorption_surface(minimised)
    logger.info('Saved plots of correction factors. \n')

    '''Clean up reflection table for outputting and save data'''
  minimised.clean_reflection_table()
  minimised.save_reflection_table(params.output.scaled_out)
  logger.info(('\nSaved reflection table to {0}').format(params.output.scaled_out))

  '''All done!'''
  finish_time = time.time()
  logger.info("\nTotal time taken: {0:.4f}s ".format((finish_time - start_time)))
  logger.info('\n'+'*'*40+'\n')

def lbfgs_scaling(scaler):
  '''main function to do lbfbs scaling'''

  if isinstance(scaler, ScalerFactory.TargetScaler):
    '''do a scaling round against a target of already scaled datasets'''
    param_name = ParameterHandler.ParameterlistFactory.full_active_list(scaler.params)
    '''Pass the scaler to the optimiser'''
    scaler = mf.LBFGS_optimiser(scaler,
      param_name=param_name[0]).return_scaler()

    '''the minimisation has only been done on a subset on the data, so apply the
    scale factors to the sorted reflection table.'''
    scaler.expand_scales_to_all_reflections()
    if scaler.params.scaling_options.only_target is True:
      scaler.join_multiple_datasets()
      return scaler
    '''now pass to a multiscaler ready for next round of scaling.'''
    scaler = ScalerFactory.MultiScalerFactory.create_from_targetscaler(scaler)

  '''from here onwards, scaler should only be a SingleScaler
  or MultiScaler (not TargetScaler)'''
  if scaler.params.scaling_options.concurrent_scaling:
    param_name = ParameterHandler.ParameterlistFactory.full_active_list(scaler.params)
  else:
    param_name = ParameterHandler.ParameterlistFactory.consecutive_list(scaler.params)
  for param in param_name:
    '''Pass the scaler to the optimiser'''
    scaler = mf.LBFGS_optimiser(scaler,
      param_name=param).return_scaler()
  '''Optimise the error model and then do another minimisation'''
  if scaler.params.weighting.optimise_error_model:
    scaler.update_error_model()
    for param in param_name:
      scaler = mf.LBFGS_optimiser(scaler,
        param_name=param).return_scaler()

  '''The minimisation has only been done on a subset on the data, so apply the
  scale factors to the whole reflection table.'''
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
