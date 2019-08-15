"""
A cross validation program for scaling.

The program runs dials.scale with each option in turn, using a free set to
score the model - the results are printed in a table and the model with the
lowest free set rmsd is indicated. For each option, the analysis will be
repeated nfolds times, with a different free set chosen each time, and the
final rmsds averaged. For full k-fold cross validation, nfolds should be set to
100/free_set_percentage, which would be nfolds=10 for the default
free_set_percentage=10.0.

Two different modes are currently supported, controlled by mode=;
1) mode=single
   dials.scale in run nfolds times for the user specified dials.scale options
2) mode=multi
   optimise a dials.scale parameter, specified by parameter= .
   parameter_values must also be specified as a string of space separated values,
   unless the dials.scale parameter type is bool.

Therefore one must choose:
  mode=              (single or multi)
  parameter=         (a supported dials.scale command line option, optional
                      if mode=single)
  parameter_values=  (values to test, only optional if parameter= selectes a
                      boolean command-line parameter)

For example
mode=multi parameter=absorption_term
mode=multi parameter=decay_interval parameter_values="5.0 10.0 15.0"
mode=multi parameter=model parameter_values="array physical"
"""

from __future__ import absolute_import, division, print_function
import logging
import itertools
import time
import sys
from copy import deepcopy
from libtbx import phil
from libtbx.table_utils import simple_table
from scitbx.array_family import flex
from dials.util import log, show_mail_on_error
from dials.util.options import OptionParser, flatten_reflections, flatten_experiments
from dials.util.version import dials_version
from dials.command_line.scale import Script

phil_scope = phil.parse(
    """
  log = dials.cross_validate.log
      .type = str
      .help = "The log filename"
  debug_log = dials.cross_validate.debug.log
      .type = str
      .help = "The debug log filename"
  cross_validation {
    mode = multi *single
      .type = choice
      .help = "Choose the cross validation running mode, for a full description"
              "see the module docstring. Choice is used for testing a parameter"
              "that can only have discreet values (a choice or bool phil parameter)."
              "Variable is used for testing a parameter that can have a float or"
              "int value (that is also not a 'choice' type). Single just performs"
              "cross validation on one parameter configuration."
    parameter = None
      .type = str
      .help = "Optimise a command-line parameter. parameter_values must also be"
              "specified, unless the parameter is a True/False option."
    parameter_values = None
      .type = strings
      .help = "Parameter values to compare, entered as a string of space"
              "separated values."
    nfolds = 1
      .type = int(value_min=1)
      .help = "Number of cross-validation folds to perform. If nfolds > 1, the"
              "minimisation for each option is repeated nfolds times, with an"
              "incremental offset for the free set. The max number of folds"
              "allowed is 1/free_set_percentage; if set greater than this then"
              "the repetition will finish afer 1/free_set_percentage folds."
  }
  include scope dials.command_line.scale.phil_scope
""",
    process_includes=True,
)

logger = logging.getLogger("dials")
info_handle = log.info_handle(logger)


def cross_validate():
    """Run cross validation script."""
    optionparser = OptionParser(
        usage=__doc__.strip(),
        read_experiments=True,
        read_reflections=True,
        read_datablocks=False,
        phil=phil_scope,
        check_format=False,
    )
    params, _ = optionparser.parse_args(show_diff_phil=False)
    if not params.input.experiments or not params.input.reflections:
        optionparser.print_help()
        sys.exit()

    log.config(verbosity=1, info=params.log, debug=params.debug_log)
    logger.info(dials_version())
    diff_phil = optionparser.diff_phil
    if diff_phil.as_str() is not "":
        logger.info("The following parameters have been modified:\n")
        logger.info(diff_phil.as_str())

    diff_phil.objects = [
        obj
        for obj in diff_phil.objects
        if not (obj.name == "input" or obj.name == "cross_validation")
    ]

    reflections = flatten_reflections(params.input.reflections)
    experiments = flatten_experiments(params.input.experiments)

    options_dict = {}

    start_time = time.time()

    if params.cross_validation.mode == "single":
        # just run the setup nfolds times
        results_dict = {}
        results_dict[0] = {
            "configuration": ["user"],
            "Rwork": [],
            "Rfree": [],
            "CCwork": [],
            "CCfree": [],
        }
        for n in range(params.cross_validation.nfolds):
            if n < 100.0 / params.scaling_options.free_set_percentage:
                params.scaling_options.free_set_offset = n
                results_dict[0] = run_script(
                    params, experiments, reflections, results_dict[0]
                )

    elif params.cross_validation.mode == "multi":
        # run each option nfolds times
        if params.cross_validation.parameter is None:
            assert (
                0
            ), "parameter= must be set to specify what command line option should be optimised"

        choice = params.cross_validation.parameter
        # inspect the phil scope to see what the parameter type is - bool, choice,
        # int or float.
        typ, _ = get_parameter_type_and_value(choice, 0)

        if typ == "bool" and not params.cross_validation.parameter_values:
            # values not specified, implied that should test both True and False
            options_dict[choice] = [True, False]
        else:
            if not params.cross_validation.parameter_values:
                assert (
                    0
                ), "parameter_values= must be set to specify what options should be tested"
            options_dict[choice] = []
            if typ == "bool":
                if (
                    "true" in params.cross_validation.parameter_values
                    or "True" in params.cross_validation.parameter_values
                ):
                    options_dict[choice].append(True)
                if (
                    "false" in params.cross_validation.parameter_values
                    or "False" in params.cross_validation.parameter_values
                ):
                    options_dict[choice].append(False)
            elif typ == "choice":
                for option in params.cross_validation.parameter_values:
                    options_dict[choice].append(option)
            elif typ == "int" or typ == "float":
                for value in params.cross_validation.parameter_values:
                    # first convert str to float or int
                    _, val = get_parameter_type_and_value(choice, value)
                    options_dict[choice].append(val)
            else:
                assert 0, "Error in interpreting parameter and parameter_values"

        # this code below should work for more than one parameter to be optimised,
        # but one cannot specify this yet from the command line
        keys, values = zip(*options_dict.items())
        results_dict = {}
        for i, v in enumerate(itertools.product(*values)):
            e = dict(zip(keys, v))
            results_dict[i] = {
                "configuration": [],
                "Rwork": [],
                "Rfree": [],
                "CCwork": [],
                "CCfree": [],
            }
            for k, v in e.iteritems():
                params = set_parameter(params, k, v)
                results_dict[i]["configuration"].append(str(k) + "=" + str(v))
            for n in range(params.cross_validation.nfolds):
                if n < 100.0 / params.scaling_options.free_set_percentage:
                    params.scaling_options.free_set_offset = n
                    results_dict[i] = run_script(
                        params, experiments, reflections, results_dict[i]
                    )

    else:
        assert 0, "Error in interpreting mode and options."

    interpret_results(results_dict)
    if diff_phil.objects:
        logger.info("\nAdditional configuration for all runs: \n")
        logger.info(diff_phil.as_str())
    logger.info("\nCross-validation finished.\n")

    finish_time = time.time()
    logger.info("\nTotal time taken: {0:.4f}s ".format(finish_time - start_time))
    logger.info("\n" + "=" * 80 + "\n")


def set_parameter(params, name, val):
    """Find the name in the params scope extract and set it to the val."""
    # Note: must be a better way to do this?
    if name in ["lmax", "n_modulation_bins", "n_resolution_bins", "n_absorption_bins"]:
        params.parameterisation.__setattr__(name, int(val))  # convert float to int
    elif name in [
        "scale_term",
        "scale_interval",
        "decay_term",
        "decay_interval",
        "absorption_term",
        "surface_weight",
        "modulation_term",
    ]:
        params.parameterisation.__setattr__(name, val)
    elif name in ["optimise_errors"]:
        params.weighting.__setattr__(name, val)
    elif name in ["d_min", "d_max"]:  # But what about biasing by n_refl?
        params.cut_data.__setattr__(name, val)
    elif name in [
        "target_cycle",
        "concurrent",
        "full_matrix",
        "outlier_zmax",
        "outlier_rejection",
    ]:
        params.scaling_options.__setattr__(name, val)
    elif name in ["model"]:
        params.__setattr__(name, val)
    else:
        assert 0, "Unable to set chosen attribute " + str(name) + "=" + str(val)
    return params


def get_parameter_type_and_value(name, value):
    """Find the parameter type for a discreet phil option - bool or choice."""
    # Note - ideally could inspect the phil_scope
    if name in ["outlier_rejection", "model"]:
        return "choice", value
    elif name in [
        "absorption_term",
        "decay_term",
        "scale_term",
        "modulation_term",
        "optimise_errors",
        "full_matrix",
        "concurrent",
        "target_cycle",
    ]:
        return "bool", value
    elif name in [
        "lmax",
        "n_modulation_bins",
        "n_resolution_bins",
        "n_absorption_bins",
    ]:
        return "int", int(value)
    elif name in [
        "surface_weight",
        "scale_interval",
        "decay_interval",
        "d_min",
        "d_max",
        "outlier_zmax",
    ]:
        return "float", float(value)


def run_script(params, experiments, reflections, results_dict):
    """Run the scaling script with the params and append to results dict."""
    params.scaling_options.__setattr__("use_free_set", True)
    script = Script(
        params, experiments=deepcopy(experiments), reflections=deepcopy(reflections)
    )
    script.run(save_data=False)
    results_dict["Rwork"].append(script.scaler.final_rmsds[0])
    results_dict["Rfree"].append(script.scaler.final_rmsds[1])
    results_dict["CCwork"].append(script.scaler.final_rmsds[2])
    results_dict["CCfree"].append(script.scaler.final_rmsds[3])
    return results_dict


def interpret_results(results_dict):
    """Pass in a dict of results. Each item is a different attempt.
  Expect a configuration and final_rmsds columns. Score the data and make a
  nice table."""
    rows = []
    headers = ["option", "", "Rwork", "Rfree", "CCwork", "CCfree"]
    free_rmsds = []
    free_cc12s = []

    def avg_sd_from_list(lst):
        """simple function to get average and standard deviation"""
        arr = flex.double(lst)
        avg = round(flex.mean(arr), 5)
        std = round(arr.standard_deviation_of_the_sample(), 5)
        return avg, std

    for v in results_dict.itervalues():
        config_str = " ".join(v["configuration"])
        avg_work, std_work = avg_sd_from_list(v["Rwork"])
        avg_free, std_free = avg_sd_from_list(v["Rfree"])
        avg_ccwork, std_ccwork = avg_sd_from_list(v["CCwork"])
        avg_ccfree, std_ccfree = avg_sd_from_list(v["CCfree"])
        rows.append(
            [
                config_str,
                "mean",
                str(avg_work),
                str(avg_free),
                str(avg_ccwork),
                str(avg_ccfree),
            ]
        )
        rows.append(
            [
                "",
                "std dev",
                str(std_work),
                str(std_free),
                str(std_ccwork),
                str(std_ccfree),
            ]
        )
        free_rmsds.append(avg_free)
        free_cc12s.append(avg_ccfree)
    # find lowest free rmsd
    low_rmsd_idx = free_rmsds.index(min(free_rmsds)) * 2  # *2 to skip std rows
    high_cc12_idx = free_cc12s.index(max(free_cc12s)) * 2
    rows[low_rmsd_idx][3] += "*"
    rows[high_cc12_idx][5] += "*"
    st = simple_table(rows, headers)
    logger.info("Summary of the cross validation analysis: \n")
    logger.info(st.format())


if __name__ == "__main__":
    with show_mail_on_error():
        cross_validate()
