#!/usr/bin/env python
# coding: utf-8
"""
Usage: dials_scratch.scale integrated.pickle integrated_experiments.json
[integrated.pickle(2) integrated_experiments.json(2) ....] [options]

This program performs scaling on the input datasets. The default parameterisation
is a physical parameterisation based on that used in the program Aimless. If
multiple input files have been specified, the datasets will be jointly scaled
against a common target of unique reflection intensities.

By default, a scale, decay and absorption correction parameterisation for each
dataset is used. One scaled.pickle and scaled_experiments.json files are output,
which may contain data and scale models from multiple experiments. The reflection
intensities are left unscaled and unmerged in the output, but an
'inverse_scale_factor' column is added to the pickle file.

To plot the scale factors determined by this program, one should subsequently run:
dials_scratch.plot_scaling_models scaled.pickle scaled_experiments.json
"""
from __future__ import absolute_import, division, print_function
import time
import logging
import sys
import libtbx.load_env
from libtbx import phil
from dials.util import halraiser
from dials.util.options import OptionParser, flatten_reflections, flatten_experiments

from dials_scratch.jbe.scaling_code.minimiser_functions import LBFGS_optimiser
from dials_scratch.jbe.scaling_code.ScalingRefiner import \
  AdaptLbfgs, GaussNewtonIterations, LevenbergMarquardtIterations
from dials_scratch.jbe.scaling_code.model import ScalingModelFactory
from dials_scratch.jbe.scaling_code import ScalerFactory
from dials_scratch.jbe.scaling_code import Scaler
from dials_scratch.jbe.scaling_code import ParameterHandler
from dials_scratch.jbe.scaling_code.target_function import ScalingTarget, ScalingTargetFixedIH
from dials_scratch.jbe.scaling_code.scaling_utilities import (
  parse_multiple_datasets, save_experiments)


start_time = time.time()
logger = logging.getLogger('dials')

phil_scope = phil.parse('''
  debug = False
    .type = bool
    .help = "Output additional debugging information"
  scaling_model = aimless
      .type = str
      .help = "Set method for scaling - 'aimless', 'xscale' or 'KB'. "
  output {
    log = dials_scratch.scaling.log
      .type = str
      .help = "The log filename"
    debug_log = dials_scratch.scaling.debug.log
      .type = str
      .help = "The debug log filename"
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

  diff_phil = optionparser.diff_phil.as_str()
  if diff_phil is not '':
    logger.info('The following parameters have been modified:\n')
    logger.info(diff_phil)

  reflections = flatten_reflections(params.input.reflections)
  experiments = flatten_experiments(params.input.experiments)

  if len(experiments) != 1:
    logger.info(('Checking for the existence of a reflection table containing {sep}'
    'multiple scaled datasets. {sep}').format(sep='\n'))
    reflections = parse_multiple_datasets(reflections)
    logger.info("Found %s experiments in total." % len(experiments))
    logger.info('\n'+'*'*40)

  assert len(experiments) == len(reflections), ''' mismatched number of
  experiments and reflection tables found'''

  '''first create the scaling model if it didn't already exist in the
  experiments files'''
  experiments = ScalingModelFactory.Factory.create(params, experiments, reflections)

  '''now create the scaler and do the scaling'''
  scaler = ScalerFactory.Factory.create(params, experiments, reflections)
  minimised = scaling_algorithm(scaler)

  for experiment in experiments:
    experiment.scaling_model.set_scaling_model_as_scaled()

  if minimised.outlier_table:
    minimised.save_outlier_table('outliers.pickle')
    logger.info("%s outliers in total were removed from the dataset and saved to %s"
      % (len(minimised.outlier_table), 'outliers.pickle'))

  '''calculate merging stats'''
  results, scaled_ids = minimised.calc_merging_statistics()
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
        logger.info("\nStatistics for dataset " + str(scaled_ids[i]))
        result.overall.show_summary()
        plot_labels.append('Dataset ' + str(scaled_ids[i]))
      else:
        logger.info("\nStatistics for combined datasets")
        result.overall.show_summary()
        plot_labels.append('Combined datasets')

  '''plot merging stats if requested'''
  if params.output.plot_merging_stats:
    from xia2.command_line.compare_merging_stats import plot_merging_stats
    plot_merging_stats(results, labels=plot_labels)

  #correl_list = minimised.calc_correlation()
  #if correl_list:
  #  n = len(minimised.single_scalers)
  #  print(correl_list)
  #  for i in range(0, n*n, n):
  #    print(correl_list[i:i+n])

  '''save scaled_experiments.json file'''
  save_experiments(experiments, params.output.experiments_out)

  '''Save scaled.pickle datafile'''
  minimised.clean_reflection_table()
  minimised.save_reflection_table(params.output.scaled_out)
  logger.info(('\nSaved reflection table to {0}').format(params.output.scaled_out))

  '''All done!'''
  finish_time = time.time()
  logger.info("\nTotal time taken: {0:.4f}s ".format((finish_time - start_time)))
  logger.info('\n'+'*'*40+'\n')

def perform_scaling(scaler, target_type=ScalingTarget):
  '''a function to call to do a complete minimisation based on the current state'''
  apm_factory = ParameterHandler.ActiveParameterFactory.create(scaler)
  for _ in range(apm_factory.n_cycles):
    apm = apm_factory.make_next_apm()
    refinery = AdaptLbfgs(scaler, target=target_type(scaler, apm),
      prediction_parameterisation=apm, max_iterations=25)
    refinery.run()
    scaler = refinery.return_scaler()
  return scaler

def scaling_algorithm(scaler):
  '''main function to do lbfbs scaling'''
  if scaler.id_ == 'target':
    '''do a scaling round against a target of already scaled datasets'''
    scaler = perform_scaling(scaler, target_type=ScalingTargetFixedIH)

    '''the minimisation has only been done on a subset on the data, so apply the
    scale factors to the whole reflection table.'''
    scaler.expand_scales_to_all_reflections()
    if scaler.params.scaling_options.only_target is True:
      scaler.join_multiple_datasets()
      return scaler
    '''now pass to a multiscaler ready for next round of scaling.'''
    scaler = ScalerFactory.MultiScalerFactory.create_from_targetscaler(scaler)

  '''from here onwards, scaler should only be a SingleScaler
  or MultiScaler (not TargetScaler)'''
  scaler = perform_scaling(scaler)
  '''Optimise the error model and then do another minimisation'''
  if scaler.params.weighting.optimise_error_model:
    scaler.update_error_model()
    scaler = perform_scaling(scaler)

  apm_factory = ParameterHandler.ActiveParameterFactory.create(scaler)
  for _ in range(apm_factory.n_cycles):
    apm = apm_factory.make_next_apm()
    refinery = GaussNewtonIterations(scaler, target=ScalingTarget(scaler, apm),
      prediction_parameterisation=apm, max_iterations=1)
    refinery.run()
    scaler = refinery.return_scaler()
  
  '''The minimisation has only been done on a subset on the data, so apply the
  scale factors to the whole reflection table.'''
  scaler.expand_scales_to_all_reflections()
  if isinstance(scaler, Scaler.MultiScalerBase):
    scaler.join_multiple_datasets()
  return scaler

if __name__ == "__main__":
  try:
    sys.exit(main(sys.argv[1:]))
  except Exception as e:
    halraiser(e)
