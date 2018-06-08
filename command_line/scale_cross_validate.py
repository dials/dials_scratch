
"""
A cross validation program for scaling.
"""
from __future__ import absolute_import, division, print_function
import logging
import itertools
from copy import deepcopy
from libtbx import phil
from libtbx.table_utils import simple_table
from dials.util import halraiser, log
from dials.util.options import OptionParser, flatten_reflections,\
  flatten_experiments
from dials.command_line.scale import Script


phil_scope = phil.parse('''
  cross_validation {
    toggle_components = None
      .type = strings
      .help = "Phil options to turn on and off for cross validation e.g.
               decay_term"
  }
  include scope dials.command_line.scale.phil_scope
''', process_includes=True)

logger = logging.getLogger('dials')
info_handle = log.info_handle(logger)

def cross_validate():
  "run cross valudation"
  optionparser = OptionParser(usage=__doc__.strip(), read_experiments=True,
      read_reflections=True, read_datablocks=False, phil=phil_scope,
      check_format=False)
  params, _ = optionparser.parse_args(show_diff_phil=False)

  diff_phil = optionparser.diff_phil
  diff_phil.objects = [obj for obj in diff_phil.objects if not (
    obj.name == 'input' or obj.name == 'cross_validation')]

  reflections = flatten_reflections(params.input.reflections)
  experiments = flatten_experiments(params.input.experiments)

  options_dict = {}
  for option in params.cross_validation.toggle_components:
    options_dict[option] = [True, False]
  keys, values = zip(*options_dict.items())
  results_dict = {}
  for i, v in enumerate(itertools.product(*values)):
    e = dict(zip(keys, v))
    results_dict[i] = {"configuration": []}
    for k, v in e.iteritems():
      params.parameterisation.__setattr__(k, v)
      results_dict[i]["configuration"].append(str(k)+'='+str(v))

    #code specific to scaling
    params.scaling_options.__setattr__("use_free_set", True)
    script = Script(params, experiments=deepcopy(experiments),
      reflections=deepcopy(reflections))
    script.run(save_data=False)
    results_dict[i]["final_rmsds"] = script.minimised.final_rmsds

  interpret_results(results_dict)
  if diff_phil.objects:
    logger.info("\nAdditional configuration for all runs: \n")
    logger.info(diff_phil.as_str())
  logger.info("\nCross-validation finished.\n")

def interpret_results(results_dict):
  """Pass in a dict of results. Each item is a different attempt.
  Expect a configuration and final_rmsds columns. Score the data and make a
  nice table."""

  rows = []
  headers = ['option', 'work_rmsd', 'free_rmsd']
  free_rmsds = []
  for v in results_dict.itervalues():
    free_rmsds.append(v['final_rmsds'][2])
    config_str = ' '.join(v['configuration'])
    rows.append([config_str, str(round(v['final_rmsds'][1], 5)),
      str(round(v['final_rmsds'][2], 5))])
  #find lowest score and free rmsd
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
