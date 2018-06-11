
"""
A cross validation program for scaling.

The program runs dials.scale with each option in turn, using a free set to
score the model - the results are printed in a table and the model with the
lowest free set rmsd is indicated. For each option, the analysis will be
repeated nfolds times, with a different free set chosen each time, and the
final rmsds averaged. For full k-fold cross validation, nfolds should be set to
100/free_set_percentage, which would be nfolds=10 for the default
free_set_percentage=10.0.

There are currently three types of options to use;
1) Toggling a list of options on and off with the command toggle_bools="a b c",
  where a, b, c etc are boolean command line options.
2) optimising a single parameter with the command optimise_parameter=x, where
  x is the name of the command line parameter for dials.scale. A list of
  parameter_values must also be supplied.
3) Optimising a choice with the command optimise_choice=x, where
  x is the name of the command line parameter for dials.scale which has a
  type = choice. A list of choice_values=a,b etc must also be supplied.
"""

from __future__ import absolute_import, division, print_function
import logging
import itertools
import time
from copy import deepcopy
from libtbx import phil
from libtbx.table_utils import simple_table
from dials.util import halraiser, log
from dials.util.options import OptionParser, flatten_reflections,\
  flatten_experiments
from dials.command_line.scale import Script, phil_scope

phil_scope = phil.parse('''
  log = dials.cross_validate.log
      .type = str
      .help = "The log filename"
  debug_log = dials.cross_validate.debug.log
      .type = str
      .help = "The debug log filename"
  cross_validation {
    toggle_bools = None
      .type = strings
      .help = "Command-line options to turn on and off for cross validation e.g.
               decay_term"
    optimise_choice = None
      .type = str
      .help = "Command-line option to turn on and off for cross validation e.g.
               decay_term"
    choice_values = None
      .type = strings
      .help = "Choice values to try."
    optimise_parameter = None
      .type = str
      .help = "Optimise a command-line parameter. sample_values must also be
               specified."
    parameter_values = None
      .type = floats
      .help = "Parameter values to compare."
    nfolds = 1
      .type = int
      .help = "Number of cross-validation folds to perform. If nfolds > 1, the
               minimisation for each option is repeated nfolds times, with an
               incremental offset for the free set. The max number of folds
               allowed is 1/free_set_percentage; if set greater than this then
               the repetition will finish afer 1/free_set_percentage folds."
  }
  include scope dials.command_line.scale.phil_scope
''', process_includes=True)

logger = logging.getLogger('dials')
info_handle = log.info_handle(logger)

def cross_validate():
  """Run cross validation script."""
  optionparser = OptionParser(usage=__doc__.strip(), read_experiments=True,
      read_reflections=True, read_datablocks=False, phil=phil_scope,
      check_format=False)
  params, _ = optionparser.parse_args(show_diff_phil=False)

  diff_phil = optionparser.diff_phil
  diff_phil.objects = [obj for obj in diff_phil.objects if not (
    obj.name == 'input' or obj.name == 'cross_validation')]

  reflections = flatten_reflections(params.input.reflections)
  experiments = flatten_experiments(params.input.experiments)

  log.config(verbosity=1, info=params.log,
      debug=params.debug_log)

  options_dict = {}

  if [bool(params.cross_validation.toggle_bools),
    bool(params.cross_validation.optimise_choice),
    bool(params.cross_validation.optimise_parameter)].count(True) > 1:
    assert 0, """Cannot set more than one of toggle_bools, optimise_choice and
      optimise_parameter simultaneously, please choose only one of these and retry."""

  start_time = time.time()

  if params.cross_validation.toggle_bools:
    for option in params.cross_validation.toggle_bools:
      options_dict[option] = [True, False]
    keys, values = zip(*options_dict.items())
    results_dict = {}
    for i, v in enumerate(itertools.product(*values)):
      e = dict(zip(keys, v))
      results_dict[i] = {"configuration": [], "work_rmsd": [], "free_rmsd": []}
      for k, v in e.iteritems():
        params = set_parameter(params, k, v)
        results_dict[i]["configuration"].append(str(k)+'='+str(v))
      for n in range(params.cross_validation.nfolds):
        if n < 100.0/params.scaling_options.free_set_percentage:
          params.scaling_options.free_set_offset = n
          results_dict[i] = run_script(params, experiments, reflections,
            results_dict[i])

  elif params.cross_validation.optimise_choice:
    assert params.cross_validation.choice_values, """choice_values must be
      specified."""
    results_dict = {}
    for i, value in enumerate(params.cross_validation.choice_values):
      k = params.cross_validation.optimise_choice
      params = set_parameter(params, k, value)
      results_dict[i] = {"configuration": [str(k)+'='+str(value)],
        "work_rmsd": [], "free_rmsd": []}
      for n in range(params.cross_validation.nfolds):
        if n < 100.0/params.scaling_options.free_set_percentage:
          params.scaling_options.free_set_offset = n
          results_dict[i] = run_script(params, experiments, reflections,
            results_dict[i])

  elif params.cross_validation.optimise_parameter:
    assert params.cross_validation.parameter_values, """parameter_values must be
      specified."""
    results_dict = {}
    for i, value in enumerate(params.cross_validation.parameter_values):
      k = params.cross_validation.optimise_parameter
      params = set_parameter(params, k, value)
      results_dict[i] = {"configuration": [str(k)+'='+str(value)],
        "work_rmsd": [], "free_rmsd": []}
      for n in range(params.cross_validation.nfolds):
        if n < 100.0/params.scaling_options.free_set_percentage:
          params.scaling_options.free_set_offset = n
          results_dict[i] = run_script(params, experiments, reflections,
            results_dict[i])

  interpret_results(results_dict)
  if diff_phil.objects:
    logger.info("\nAdditional configuration for all runs: \n")
    logger.info(diff_phil.as_str())
  logger.info("\nCross-validation finished.\n")

  finish_time = time.time()
  logger.info("\nTotal time taken: {0:.4f}s ".format(finish_time - start_time))
  logger.info('\n'+'='*80+'\n')

def set_parameter(params, name, val):
  """Find the name in the params scope extract and set it to the val."""
  #Note: must be a better way to do this?
  if name in ['lmax', 'n_modulation_bins', 'n_resolution_bins',
    'n_absorption_bins']:
    params.parameterisation.__setattr__(name, int(val)) #convert float to int
  elif name in ['scale_term', 'scale_interval', 'decay_term', 'decay_interval',
    'absorption_term', 'surface_weight', 'modulation_term']:
    params.parameterisation.__setattr__(name, val)
  elif name in ['optimise_errors']:
    params.weighting.__setattr__(name, val)
  elif name in ['d_min', 'd_max']: #But what about biasing by n_refl?
    params.cut_data.__setattr__(name, val)
  elif name in ['target_cycle', 'concurrent', 'full_matrix', 'outlier_zmax',
    'outlier_rejection']:
    params.scaling_options.__setattr__(name, val)
  else:
    assert 0, "Unable to set chosen attribute " + str(name)+"="+str(val)
  return params

def run_script(params, experiments, reflections, results_dict):
  """Run the scaling script with the params and append to results dict."""
  params.scaling_options.__setattr__("use_free_set", True)
  script = Script(params, experiments=deepcopy(experiments),
    reflections=deepcopy(reflections))
  script.run(save_data=False)
  results_dict["work_rmsd"].append(script.minimised.final_rmsds[1])
  results_dict["free_rmsd"].append(script.minimised.final_rmsds[2])
  return results_dict

def interpret_results(results_dict):
  """Pass in a dict of results. Each item is a different attempt.
  Expect a configuration and final_rmsds columns. Score the data and make a
  nice table."""
  rows = []
  headers = ['option', 'work_rmsd', 'free_rmsd']
  free_rmsds = []
  for v in results_dict.itervalues():
    config_str = ' '.join(v['configuration'])
    avg_work = round(sum(v['work_rmsd'])/len(v['work_rmsd']), 5)
    avg_free = round(sum(v['free_rmsd'])/len(v['free_rmsd']), 5)
    rows.append([config_str, str(avg_work), str(avg_free)])
    free_rmsds.append(avg_free)
  #find lowest free rmsd
  low_rmsd_idx = free_rmsds.index(min(free_rmsds))
  rows[low_rmsd_idx][2] += '*'
  st = simple_table(rows, headers)
  logger.info('Summary of the cross validation analysis: \n')
  logger.info(st.format())

if __name__ == "__main__":
  try:
    cross_validate()
  except Exception as e:
    halraiser(e)
