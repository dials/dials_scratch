from __future__ import print_function
import copy
from dials.array_family import flex
from cctbx import miller, crystal
import dials_scratch.jbe.scaling_code.minimiser_functions as mf
from dials.util.options import flatten_experiments, flatten_reflections
import numpy as np
import cPickle as pickle
from dials_scratch.jbe.scaling_code.target_function import *
from dials_scratch.jbe.scaling_code.basis_functions import *
from dials_scratch.jbe.scaling_code.scaling_utilities import sph_harm_table
import dials_scratch.jbe.scaling_code.scale_factor as SF
from dials_scratch.jbe.scaling_code.reflection_weighting import *
from dials_scratch.jbe.scaling_code.data_quality_assessment import R_pim_meas
from dials_scratch.jbe.scaling_code.target_Ih import SingleIhTable
import matplotlib.pyplot as plt
from dials_scratch.jbe.scaling_code.data_manager_functions import AimlessDataManager
from dials_scratch.jbe.scaling_code.data_plotter import (
    plot_data_decay,
    plot_data_absorption,
    plot_data_modulation,
    plot_smooth_scales,
    plot_absorption_surface,
)
from dials.util.options import OptionParser
from libtbx import phil
from dxtbx.model.experiment_list import ExperimentList
from dxtbx.model import Crystal, Scan, Beam, Goniometer


class test_data_manager(AimlessDataManager):
    def __init__(self, reflections, experiments, params):
        super(test_data_manager, self).__init__(reflections, experiments, params)
        """self._reflection_table = reflections
    self.miller_set = miller_set
    self._initial_keys = [key for key in self.reflection_table.keys()]
    self._reflection_table['inverse_scale_factor'] = flex.double(
      [1.0] * len(self.reflection_table))
    self._reflection_table['Ih_values'] = flex.double([0.0] * len(self.reflection_table))
    self._reflection_table = self.map_indices_to_asu(self.reflection_table)
    'assign initial weights (will be statistical weights at this point)'
    self.weights_for_scaling = self._update_weights_for_scaling(
      self.reflection_table, weights_filter=False)
    #aimless initialisation
    self.g_absorption = None
    self.g_scale = None
    self.g_decay = None
    #self.n_active_params = 0
    self.g_parameterisation = {}
    #self.active_parameters = flex.double([])
    self.scaling_options = {}
    self.scaling_options['lmax'] = 6
    self.scaling_options['Isigma_min'] = 0.0
    self.scaling_options['d_min'] = 0.0
    self.scaling_options['scaling_method'] = 'aimless'
    self.scaling_options['multi_mode'] = False
    self.scaling_options['decay_term'] = False
    self.scaling_options['scale_term'] = True
    self.scaling_options['absorption_term'] = True
    '''bin reflections, determine outliers, extract reflections and weights for
    scaling and set normalised values.'''
    self.initialise_scale_factors(self.reflection_table)
    (reflections_for_scaling, weights_for_scaling) = (
      self.extract_reflections_for_scaling(self.reflection_table))
    self.Ih_table = SingleIhTable(reflections_for_scaling, weights_for_scaling.get_weights())
    '''refactor the next two operations into extract_reflections?
    reset the normalised values within the scale_factor object to current'''
    self.g_scale.set_normalised_values(reflections_for_scaling[
      'norm_rot_angle'])
    self.g_decay.set_normalised_values(reflections_for_scaling[
      'normalised_time_values'])
    self.g_decay.set_d_values(reflections_for_scaling['d'])
    self.g_absorption.set_values(sph_harm_table(reflections_for_scaling,
                                                self.scaling_options['lmax']))"""

    """def map_indices_to_asu(self, reflection_table):
    '''Create a miller_set object, map to the asu and create a sorted
       reflection table, sorted by asu miller index'''
    reflection_table["asu_miller_index"] = self.miller_set.map_to_asu().indices()
    permuted = (self.miller_set.map_to_asu()).sort_permutation(by_value='packed_indices')
    reflection_table = reflection_table.select(permuted)
    return reflection_table

  def initialise_scale_factors(self, reflection_table):
    '''initialise scale factors and add to self.active_parameters'''
    self.initialise_scale_term(reflection_table)
    self.initialise_decay_term(reflection_table)
    self.initialise_absorption_scales(reflection_table, self.scaling_options['lmax'])

  def initialise_decay_term(self, reflection_table):
    '''calculate the relative, normalised rotation angle. Here this is called
    normalised time to allow a different rotation interval compare to the scale
    correction. A SmoothScaleFactor_1D object is then initialised'''
    rotation_interval = 15.0
    rotation_interval = rotation_interval + 0.001
    reflection_table['normalised_time_values'] = reflection_table['xyz'].parts()[2]/rotation_interval
    'define the highest and lowest gridpoints: go out two further than the max/min int values'
    highest_parameter_value = int((max(reflection_table['normalised_time_values'])//1)+3)#was +2
    lowest_parameter_value = int((min(reflection_table['normalised_time_values'])//1)-2)#was -1
    n_decay_parameters =  highest_parameter_value - lowest_parameter_value + 1
    self.g_decay = SF.SmoothBScaleFactor1D(0.0, n_decay_parameters, reflection_table['d'])
    #self.g_decay.set_normalised_values(reflection_table['normalised_time_values'])
    self.g_parameterisation['g_decay'] = self.g_decay

  def initialise_scale_term(self, reflection_table):
    '''calculate the relative, normalised rotation angle.
    A SmoothScaleFactor_1D object is then initialised'''
    rotation_interval = 15.0
    rotation_interval = rotation_interval + 0.001
    reflection_table['norm_rot_angle'] = reflection_table['xyz'].parts()[2] / rotation_interval
    'define the highest and lowest gridpoints: go out two further than the max/min int values'
    highest_parameter_value = int((max(reflection_table['norm_rot_angle'])//1)+3)#was +2
    lowest_parameter_value = int((min(reflection_table['norm_rot_angle'])//1)-2)#was -1
    n_scale_parameters = highest_parameter_value - lowest_parameter_value + 1
    self.g_scale = SF.SmoothScaleFactor1D(1.0, n_scale_parameters)
    #self.g_scale.set_normalised_values(reflection_table['norm_rot_angle'])
    self.g_parameterisation['g_scale'] = self.g_scale

  def initialise_absorption_scales(self, reflection_table, lmax):
    n_abs_params = 0
    for i in range(lmax):
      n_abs_params += (2*(i+1))+1
    self.g_absorption = SF.SHScaleFactor(0.0, n_abs_params,
      sph_harm_table(reflection_table, lmax))
    self.g_parameterisation['g_absorption'] = self.g_absorption

  def calc_absorption_constraint(self, apm):
    idx = apm.active_parameterisation.index('g_absorption')
    start_idx = apm.cumulative_active_params[idx]
    end_idx = apm.cumulative_active_params[idx+1]
    weight = 1e2
    abs_params = apm.x[start_idx:end_idx]
    residual = (weight * (abs_params)**2)
    gradient = (2 * weight * abs_params)
    #need to make sure gradient is returned is same size as gradient calculated in target fn-
    #would only be triggered if refining absorption as same time as another scale factor.
    gradient_vector = flex.double([])
    for i, param in enumerate(apm.active_parameterisation):
      if param != 'g_absorption':
        gradient_vector.extend(flex.double([0.0]*apm.active_params_list[i]))
      elif param == 'g_absorption':
        gradient_vector.extend(gradient)
    return (residual, gradient_vector)"""


def load_data(filename):
    data_file = open(filename)
    data = pickle.load(data_file)
    data_file.close()
    return data


def run_main(reflections, experiments, params):
    # (reflections, ms) = load_data('test_dataset_mu5.pickle')
    loaded_reflections = test_data_manager(reflections, experiments, params)

    minimised = mf.LBFGS_optimiser(
        loaded_reflections, param_name=["g_scale", "g_decay", "g_absorption"]
    ).return_data_manager()
    # print list(minimised.g_absorption.inverse_scales)
    # print list(minimised.apm.active_parameters)
    minimised.expand_scales_to_all_reflections()

    Rpim, Rmeas = R_pim_meas(minimised)
    print("R_meas is %s" % (Rmeas))
    print("R_pim is %s" % (Rpim))

    plot_smooth_scales(minimised, outputfile="Smooth_scale_factors.png")
    plot_absorption_surface(minimised)
    print("Saved plots of correction factors")

    print(len(reflections))

    minimised.save_reflection_table("synthetic_scaled.pickle")
    print("Saved output to %s" % ("synthetic_scaled.pickle"))


def generate_test_input():
    (reflections, ms) = load_data("test_dataset_mu0p2_smalldetector_P4_rot0.pickle")

    # json.dump(datablock, open(datablock_json, 'w'))

    experiments = ExperimentList()
    exp_dict = {
        "__id__": "crystal",
        "real_space_a": [20.0, 0.0, 0.0],
        "real_space_b": [0.0, 20.0, 0.0],
        "real_space_c": [0.0, 0.0, 15.0],
        "space_group_hall_symbol": " P 4",
    }
    experiments.crystal = Crystal.from_dict(exp_dict)
    experiments.scan = Scan(image_range=[0, 360], oscillation=[0.0, 1.0])
    experiments.beam = Beam(s0=(0.0, 0.0, 1.01))
    experiments.goniometer = Goniometer((1.0, 0.0, 0.0))

    phil_scope = phil.parse(
        """
      include scope dials_scratch.jbe.scaling_code.scaling_options.phil_scope
  """,
        process_includes=True,
    )

    optionparser = OptionParser(phil=phil_scope, check_format=False)
    parameters, _ = optionparser.parse_args(
        args=None, quick_parse=True, show_diff_phil=False
    )
    parameters.__inject__("scaling_method", "aimless")
    parameters.scaling_options.__inject__("multi_mode", False)
    return (reflections, experiments, parameters)


##########
if __name__ == "__main__":
    (reflections, ms) = load_data("test_dataset_mu0p2_smalldetector_P4_rot0.pickle")

    reflections, experiments, params = generate_test_input()

    run_main(reflections, experiments, params)
