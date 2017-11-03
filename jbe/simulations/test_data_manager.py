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
from dials_scratch.jbe.scaling_code.data_quality_assessment import R_meas, R_pim
from dials_scratch.jbe.scaling_code.target_Ih import single_Ih_table
import matplotlib.pyplot as plt
from dials_scratch.jbe.scaling_code.data_manager_functions import aimless_Data_Manager
from dials_scratch.jbe.scaling_code.data_plotter import (plot_data_decay,
plot_data_absorption, plot_data_modulation, plot_smooth_scales, plot_absorption_surface)

class test_data_manager(aimless_Data_Manager):
  def __init__(self, reflections, miller_set):
    self.reflection_table = reflections
    self.miller_set = miller_set
    self.initial_keys = [key for key in self.reflection_table.keys()]
    self.reflection_table['inverse_scale_factor'] = flex.double(
      [1.0] * len(self.reflection_table))
    self.reflection_table['Ih_values'] = flex.double([0.0] * len(self.reflection_table))
    self.reflection_table = self.map_indices_to_asu(self.reflection_table)
    'assign initial weights (will be statistical weights at this point)'
    self.weights_for_scaling = self.update_weights_for_scaling(
      self.reflection_table, weights_filter=False)
    #aimless initialisation
    self.g_absorption = None
    self.g_scale = None
    self.g_decay = None
    self.n_active_params = 0
    self.active_parameters = flex.double([])
    self.scaling_options = {}
    self.scaling_options['lmax'] = 4
    self.scaling_options['Isigma_min'] = 0.0
    self.scaling_options['d_min'] = 0.0
    self.scaling_options['scaling_method'] = 'aimless'
    '''bin reflections, determine outliers, extract reflections and weights for
    scaling and set normalised values.'''
    self.initialise_scale_factors(self.reflection_table)
    (reflections_for_scaling, weights_for_scaling) = (
      self.extract_reflections_for_scaling(self.reflection_table))
    self.Ih_table = single_Ih_table(reflections_for_scaling, weights_for_scaling.get_weights())
    '''refactor the next two operations into extract_reflections?
    reset the normalised values within the scale_factor object to current'''
    self.g_scale.set_normalised_values(reflections_for_scaling[
      'normalised_rotation_angle'])
    self.g_decay.set_normalised_values(reflections_for_scaling[
      'normalised_time_values'])
    self.g_decay.set_d_values(reflections_for_scaling['d'])
    self.g_absorption.set_values(sph_harm_table(reflections_for_scaling,
                                                self.scaling_options['lmax']))

  def map_indices_to_asu(self, reflection_table):
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
    self.active_parameters.extend(self.g_scale.get_scale_factors())
    print "extended params by %s g_scale scale factors" % (len(self.g_scale.get_scale_factors()))
    self.active_parameters.extend(self.g_decay.get_scale_factors())
    print "extended params by %s g_decay scale factors" % (len(self.g_decay.get_scale_factors()))
    self.active_parameters.extend(self.g_absorption.get_scale_factors())
    print "extended params by %s g_absorption scale factors" % (len(self.g_absorption.get_scale_factors()))

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
    self.g_decay = SF.SmoothScaleFactor_1D_Bfactor(0.0, n_decay_parameters, reflection_table['d'])
    #self.g_decay.set_normalised_values(reflection_table['normalised_time_values'])
    self.n_active_params += n_decay_parameters
    self.n_g_decay_params = n_decay_parameters

  def initialise_scale_term(self, reflection_table):
    '''calculate the relative, normalised rotation angle.
    A SmoothScaleFactor_1D object is then initialised'''
    rotation_interval = 15.0
    rotation_interval = rotation_interval + 0.001
    reflection_table['normalised_rotation_angle'] = reflection_table['xyz'].parts()[2] / rotation_interval
    'define the highest and lowest gridpoints: go out two further than the max/min int values'
    highest_parameter_value = int((max(reflection_table['normalised_rotation_angle'])//1)+3)#was +2
    lowest_parameter_value = int((min(reflection_table['normalised_rotation_angle'])//1)-2)#was -1
    n_scale_parameters = highest_parameter_value - lowest_parameter_value + 1
    self.g_scale = SF.SmoothScaleFactor_1D(1.0, n_scale_parameters)
    #self.g_scale.set_normalised_values(reflection_table['normalised_rotation_angle'])
    self.n_active_params += n_scale_parameters
    self.n_g_scale_params = n_scale_parameters

  def initialise_absorption_scales(self, reflection_table, lmax):
    n_abs_params = 0
    for i in range(lmax):
      n_abs_params += (2*(i+1))+1
    self.g_absorption = SF.SphericalAbsorption_ScaleFactor(0.0, n_abs_params,
      sph_harm_table(reflection_table, lmax))
    self.n_active_params += n_abs_params
    self.n_g_abs_params = n_abs_params

  def calc_absorption_constraint(self):
    n_g_scale = self.n_g_scale_params
    n_g_decay = self.n_g_decay_params
    return (0.0 * (self.active_parameters[n_g_scale + n_g_decay:])**2)


def load_data(filename):
  data_file = open(filename)
  data = pickle.load(data_file)
  data_file.close()
  return data

(reflections, ms) = load_data('test_dataset_mu5.pickle')

loaded_reflections = test_data_manager(reflections, ms)

minimised = mf.LBFGS_optimiser(loaded_reflections,
                                        param_name=None).return_data_manager()
#print list(minimised.active_parameters)
minimised.expand_scales_to_all_reflections()

Rmeas = R_meas(minimised)
Rpim = R_pim(minimised)
print "R_meas is %s" % (Rmeas)
print "R_pim is %s" % (Rpim)

plot_smooth_scales(minimised, outputfile='Smooth_scale_factors.png')
plot_absorption_surface(minimised)
print "Saved plots of correction factors"

minimised.save_reflection_table('integrated_scaled.pickle')
print "Saved output to %s" % ('integrated_scaled.pickle')