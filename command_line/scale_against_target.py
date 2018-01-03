#!/usr/bin/env python
# coding: utf-8

"""
Usage:
  dials_scratch.xds_scaling integrated.pickle integrated_experiments.json
  target_integrated_scaled.pickle [options]

A number of options can be specified, see the phil_scope below.
"""

from __future__ import absolute_import, division, print_function
import libtbx.load_env
import logging
logger = logging.getLogger('dials.scale')

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
    log = 'dials_scratch.aimless_scaling.log'
      .type = str
      .help = "The log filename"
    debug_log = 'dials_scratch.aimless_scaling.debug.log'
      .type = str
      .help = "The debug log filename"
    plot_scalefactors = True
      .type = bool
      .help = "Option to switch off scalefactor plotting."
  }
  parameterisation {
    scale_term = True
      .type = bool
      .help = "Option to turn off decay correction (only for KB scaling)"
    rotation_interval = 15.0
      .type = float
      .help = "User specified rotation (phi) interval in degrees for phi binning
              for the scale term"
    decay_term = True
      .type = bool
      .help = "Option to turn off decay correction"
    B_factor_interval = 20.0
      .type = float
      .help = "User specified rotation (phi) interval in degrees for phi binning
              for the decay term"
    absorption_term = True
      .type = bool
      .help = "Option to turn off absorption correction"
    lmax = 4
      .type = int
      .help = "Number of spherical harmonics to include for absorption correction,
              recommended to be no more than 6."
  }
  reflection_selection {
    E2min = 0.8
      .type = float
      .help = "Minimum normalised E^2 value to select reflections for scaling"
    E2max = 5.0
      .type = float
      .help = "Maximum normalised E^2 value to select reflections for scaling"
    Isigma_min = -5.0
      .type = float
      .help = "Option to use a I/sigma subset of reflections to determine scale factors"
    d_min = 0.0
      .type = float
      .help = "Option to use a d-value subset of reflections to determine scale factors"
  }
  scaling_options {
    force_space_group = None
      .type = str
      .help = "Option to specify space group for scaling"
    concurrent_scaling = True
      .type = bool
      .help = "Option to allow absorption correction after decay/scale,
              if concurrent_scaling is set to False"
    optimise_error_model = True
      .type = bool
      .help = "Option to allow optimisation of weights for scaling. Performs
               and additional scale factor minimisation after adjusting weights."
    error_model_params = None
      .type = floats(size=2)
      .help = "Ability to force an error model adjustment, using the model
              in aimless - factors are called SDFac, SDadd in aimless."
    reject_outliers = True
      .type = bool
      .help = "Option to turn on outlier rejection"
    verbosity = 1
      .type = int(value_min=0)
      .help = "The verbosity level"
    integration_method = 'prf'
      .type = str
      .help = "Option to choose from profile fitted intensities (prf)
              or summation integrated intensities (sum)"
    minimisation_parameterisation = 'standard'
      .type = str
      .help = "Choice of 'standard' (multiplicative) or 'log' g-value
               minimisation parameterisation"
    target = None
      .type = str
      .help = "Choice to specify a target dataset for scaling"
  }
''')

from dials_scratch.jbe.scaling_code import minimiser_functions as mf
from dials_scratch.jbe.scaling_code import data_manager_functions as dmf
from dials_scratch.jbe.scaling_code.data_quality_assessment import R_pim_meas
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


  # UNWRAP all of the data objects from the PHIL parser
  reflections = flatten_reflections(params.input.reflections)
  experiments = flatten_experiments(params.input.experiments)

  if len(reflections) != 2:
    assert 0, """Incorrect number of reflection files entered
    in the command line (must be 2, one to scale and one integrated_scaled.pickle)"""
  elif len(experiments) != 1:
    assert 0, """Incorrect number of experiments files entered, only need one"""

  phil_parameters = optionparser.phil
  diff_phil_parameters = optionparser.diff_phil

  params.scaling_options.__inject__('multi_mode', False)
  params.__inject__('scaling_method', 'KB')


  '''do lbfgs minimisation'''
  minimised = scaling_lbfgs(reflections, experiments, params, logger)


  '''clean up reflection table for outputting and save data'''
  #minimised.dm1.clean_reflection_table()
  minimised.dm1.save_reflection_table('integrated_targetscaled.pickle')
  logger.info("Saved output to " + str('integrated_targetscaled.pickle'))

  logger.info('\n'+'*'*40+'\n')

def scaling_lbfgs(reflections, experiments, params, logger):
  """This algorithm performs scaling against a target scaled reflection table"""
  logger.info('\n'+'*'*40+'\n')
  loaded_reflections = dmf.targeted_datamanager(reflections[0],
    experiments[0], reflections[1], params)

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


if __name__ == "__main__":
  try:
    sys.exit(main(sys.argv[1:]))
  except Exception as e:
    halraiser(e)
