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
    self.g_parameterisation['g_decay'] = self.g_decay

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
    self.g_parameterisation['g_scale'] = self.g_scale

  def initialise_absorption_scales(self, reflection_table, lmax):
    n_abs_params = 0
    for i in range(lmax):
      n_abs_params += (2*(i+1))+1
    self.g_absorption = SF.SphericalAbsorption_ScaleFactor(0.0, n_abs_params,
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
    return (residual, gradient_vector)

def load_data(filename):
  data_file = open(filename)
  data = pickle.load(data_file)
  data_file.close()
  return data

def run_main(reflections, ms):
  #(reflections, ms) = load_data('test_dataset_mu5.pickle')


  loaded_reflections = test_data_manager(reflections, ms)

  minimised = mf.LBFGS_optimiser(loaded_reflections,
                                 param_name=['g_scale','g_absorption']).return_data_manager()
  print list(minimised.g_absorption.get_scale_factors())
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

##########
if __name__ == "__main__":
  (reflections, ms) = load_data('test_dataset_mu0p2_smalldetector_P4_rot0.pickle')

  
  '''m_idx = [(1,0,0),(0,1,0), (-1,0,0), (0,-1,0)]
  intensities = flex.double([5000.0,10000.0,5000.0,10000.0])
  variances = flex.double([1000.0,1000.0,1000.0,1000.0])
  s2d = flex.vec3_double([(1.0,0.0,0.0),(0.0,1.0,0.0),(-1.0,0.0,0.0),(0.0,-1.0,0.0)])
  ds = [2.0,2.0, 2.0,2.0]
  refls = [(1.0,0.0,0.0),(1.0,0.0,0.0),(1.0,0.0,0.0),(1.0,0.0,0.0)]
  xyzpositions = flex.vec3_double(refls)
  miller_indices = flex.miller_index(m_idx)

  reflections = flex.reflection_table()
  reflections['miller_index'] = miller_indices
  reflections['intensity'] = intensities
  reflections['variance'] = variances
  reflections['d'] = flex.double(ds)
  reflections['s2d'] = s2d
  reflections['xyz'] = xyzpositions

  ms = miller.set(crystal_symmetry=crystal.symmetry(
    space_group_symbol="P4",
        unit_cell=(6,6,6,90,90,90)),
        anomalous_flag=False,
        indices=miller_indices)'''

  run_main(reflections,ms)