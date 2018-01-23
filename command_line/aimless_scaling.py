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
start_time = time.time()
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
    experiments = "scaled_experiments.json"
      .type = str
      .help = "Option to set filepath for output json."
  }
  include scope dials_scratch.jbe.scaling_code.scaling_options.phil_scope
''', process_includes=True)


from dials_scratch.jbe.scaling_code import minimiser_functions as mf
from dials_scratch.jbe.scaling_code import data_manager_functions as dmf
from dials_scratch.jbe.scaling_code import ScalingModelFactory

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
  
  if params.scaling_options.minimisation_parameterisation not in ['standard', 'log']:
    msg = ('Invalid minimisation parameterisation choice, proceeding using {sep}'
           'using standard g-value parameterisation').format(sep='\n')
    logger.info(msg)
    params.scaling_options.minimisation_parameterisation = 'standard'

  #additional parsing to determine the type of reflection files.
  n_already_scaled = 0
  if len(reflections) > 1:
    single_reflection_tables = []
    is_already_scaled = []
    for refl_table in reflections:
      if 'dataset_id' in refl_table.keys():
        if refl_table['dataset_id'].count(0) < len(refl_table):
          dataset_ids = set(refl_table['dataset_id'])
          logger.info(('Detected existence of a multi-dataset scaled reflection table, {sep}'
            'containing {0} datasets. {sep}').format(len(dataset_ids), sep='\n'))
          for dataset_id in dataset_ids:
            single_refl_table = refl_table.select(refl_table['dataset_id'] == dataset_id)
            single_reflection_tables.append(single_refl_table)
            n_already_scaled += 1
            is_already_scaled.append(True)
          logger.info("Successfully parsed multiple scaled reflection tables. \n")
        else:
          single_reflection_tables.append(refl_table)
          is_already_scaled.append(True)
      else:
        single_reflection_tables.append(refl_table)
        is_already_scaled.append(False)
    reflections = single_reflection_tables

  len_refl = len(reflections)
  len_exp = len(experiments)

  if len_refl > 1 and len_exp > 1 and (len_exp == len_refl):
    params.scaling_options.__inject__('multi_mode', True)
  elif len_refl == 1 and len_exp == 1:
    params.scaling_options.__inject__('multi_mode', False)
  else:
    assert 0, """Incorrect number of reflection and/or experiment files entered
    in the command line: must be an equal number of each"""

  #first create the scaling model
  experiments = ScalingModelFactory.Factory.create(params, experiments, reflections)

  #do initial targeted scaling if some already scaled.
  if n_already_scaled > 0 and n_already_scaled != len_refl:
    minimised = scale_against_target(reflections, experiments, is_already_scaled, params)
  else:
    minimised = aimless_scaling_lbfgs(reflections, experiments, params)

  #calculate merging stats
  results = minimised.calc_merging_statistics()
  logger.info('*'*40)
  logger.info("Dataset statistics")
  labels = []
  #result.overall.show_summary(out=log.info_handle(logger))
  for i, result in enumerate(results):
    if len(results) == 1:
      logger.info("")
      result.overall.show_summary()
      labels.append('Single dataset ')
    else:
      if i < len(results) - 1:
        logger.info("\nStatistics for dataset " + str(i+1))
        result.overall.show_summary()
        labels.append('Dataset ' + str(i+1))
      else:
        logger.info("\nStatistics for combined datasets")
        result.overall.show_summary()
        labels.append('Combined datasets')
  if params.output.plot_merging_stats:
    from xia2.command_line.compare_merging_stats import plot_merging_stats
    plot_merging_stats(results, labels=labels)

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

  #save scaled_experiments.json file
  save_experiments(experiments, params.output.experiments)
  # Clean up reflection table for outputting and save data
  minimised.clean_reflection_table()
  minimised.save_reflection_table('integrated_scaled.pickle')
  logger.info(('\nSaved reflection table to {0}').format('integrated_scaled.pickle'))

  # All done!
  finish_time=time.time()
  logger.info("\nTotal time taken: {0:.4f}s ".format((finish_time - start_time)))
  logger.info('\n'+'*'*40+'\n')

  

def scale_against_target(reflections, experiments, is_already_scaled, params):
  """This algorithm performs scaling against a target scaled reflection table"""
  logger.info('\n'+'*'*40+'\n')
  params.scaling_options.multi_mode = False
  loaded_reflections = dmf.TargetedDataManager(reflections, experiments,
    is_already_scaled, params)

  '''call the optimiser on the Data Manager object'''
  param_name = []
  if params.parameterisation.scale_term:
    param_name.append('g_scale')
  if params.parameterisation.decay_term:
    param_name.append('g_decay')
  if params.parameterisation.absorption_term:
    param_name.append('g_absorption')
  if not param_name:
      assert 0, 'no parameters have been chosen for scaling, aborting process'
  loaded_reflections = mf.LBFGS_optimiser(loaded_reflections,
    param_name=param_name).return_data_manager()

  '''the minimisation has only been done on a subset on the data, so apply the
  scale factors to the sorted reflection table.'''
  loaded_reflections.expand_scales_to_all_reflections()

  #reflections = [i.reflection_table for i in loaded_reflections.data_managers]
  #experiments = [i.experiments for i in loaded_reflections.data_managers]
  #reflections.append(loaded_reflections.dm1.reflection_table)
  #experiments.append(loaded_reflections.dm1.experiments)
  loaded_reflections.params.scaling_options.multi_mode = True
  return aimless_scaling_lbfgs(multi_data_manager=loaded_reflections,
    reflections=None, experiments=None, params=loaded_reflections.params)#reflections, experiments, params)


def aimless_scaling_lbfgs(reflections, experiments, params, multi_data_manager=None):
  """This algorithm performs an aimless-like scaling"""
  logger.info('\n'+'*'*40+'\n')

  # Initialise the datamanager.
  if params.scaling_options.multi_mode:
    if multi_data_manager is not None:
      loaded_reflections = dmf.MultiCrystalDataManager()
      loaded_reflections.init_from_datamanager(
        multi_data_manager)
    else:
      loaded_reflections = dmf.MultiCrystalDataManager()
      loaded_reflections.init_from_refl_tables(
        reflections, experiments, params)
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
