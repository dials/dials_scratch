"""
Tests for the scaling model classes.
"""
from dials.array_family import flex
from dials.util.options import OptionParser
from dials_scratch.jbe.scaling_code.model.model import ScalingModelBase,\
  KBScalingModel, PhysicalScalingModel, ArrayScalingModel
from libtbx import phil
from dxtbx.model.experiment_list import ExperimentList
from dxtbx.model import Crystal, Scan, Beam, Goniometer, Detector, Experiment

def generated_single_exp():
  experiments = ExperimentList()
  exp_dict = {"__id__" : "crystal", "real_space_a": [1.0, 0.0, 0.0],
              "real_space_b": [0.0, 1.0, 0.0], "real_space_c": [0.0, 0.0, 2.0],
              "space_group_hall_symbol": " C 2y"}
  crystal = Crystal.from_dict(exp_dict)
  scan = Scan(image_range=[0, 90], oscillation=[0.0, 1.0])
  beam = Beam(s0=(0.0, 0.0, 1.01))
  goniometer = Goniometer((1.0, 0.0, 0.0))
  detector = Detector()
  experiments.append(Experiment(beam=beam, scan=scan, goniometer=goniometer,
    detector=detector, crystal=crystal))
  return experiments

def generated_param():
  phil_scope = phil.parse('''
      include scope dials_scratch.jbe.scaling_code.scaling_options.phil_scope
  ''', process_includes=True)

  optionparser = OptionParser(phil=phil_scope, check_format=False)
  parameters, _ = optionparser.parse_args(args=None, quick_parse=True,
    show_diff_phil=False)
  return parameters

def test_ScalingModelBase():
  """Test for base scaling model class"""

  class SM_base_filler(ScalingModelBase):
    """Fill in abstract methid"""
    def consecutive_refinement_order(self):
      pass

  SM_base = SM_base_filler(configdict={})
  assert not SM_base.is_scaled
  SM_base.set_scaling_model_as_scaled()
  assert SM_base.is_scaled
  SM_base.set_scaling_model_as_unscaled()
  assert not SM_base.is_scaled
  assert SM_base.configdict == {}
  assert not SM_base.components
  _ = SM_base.to_dict()
  SM_base.set_error_model([1.0])
  assert SM_base.configdict['error_model_parameters'] == [1.0]
  _ = SM_base.configure_reflection_table(1.0, 2.0, 3.0) #Check method exists

def test_KBScalingModel():
  """Test for the KB Scaling Model."""

  # Test standard initialisation method.
  configdict = {'corrections':['scale', 'decay']}
  parameters_dict = {
      'scale': {'parameters' : flex.double([1.2]),
                'parameter_esds' : flex.double([0.1])},
      'decay': {'parameters' : flex.double([0.01]),
                'parameter_esds' : flex.double([0.02])}}
  KBmodel = KBScalingModel(parameters_dict, configdict)
  assert KBmodel.id_ == 'KB'
  assert 'scale' in KBmodel.components
  assert 'decay' in KBmodel.components
  assert list(KBmodel.components['scale'].parameters) == [1.2]
  assert list(KBmodel.components['decay'].parameters) == [0.01]
  assert list(KBmodel.components['scale'].parameter_esds) == [0.1]
  assert list(KBmodel.components['decay'].parameter_esds) == [0.02]

  # Test from_dict initialisation method.
  KB_dict = {"__id__": "KB", "is_scaled": True, "scale": {
    "n_parameters": 1, "parameters": [0.5], "est_standard_devs" : [0.05]},
    "configuration_parameters": {"corrections": ["scale"]}}
  KBmodel = KBScalingModel.from_dict(KB_dict)
  assert KBmodel.is_scaled is True
  assert 'scale' in KBmodel.components
  assert 'decay' not in KBmodel.components
  assert list(KBmodel.components['scale'].parameters) == [0.5]
  assert list(KBmodel.components['scale'].parameter_esds) == [0.05]

  new_dict = KBmodel.to_dict()
  assert new_dict == KB_dict

def test_PhysicalScalingModel():
  """Test the PhysicalScalingModel class."""
  configdict = {'corrections':['scale', 'decay', 'absorption'],
    's_norm_fac': 1.0, 'scale_rot_interval': 2.0, 'd_norm_fac': 1.0,
    'decay_rot_interval': 2.0, 'lmax' : 1}

  parameters_dict = {
      'scale': {'parameters' : flex.double([1.2, 1.1]), 'parameter_esds' : None},
      'decay': {'parameters' : flex.double([0.1, 0.2]), 'parameter_esds' : None},
      'absorption': {'parameters' : flex.double([0.01, 0.01, 0.01]),
      'parameter_esds' : None}}

  # Test standard factory initialisation
  physicalmodel = PhysicalScalingModel(parameters_dict, configdict)
  assert physicalmodel.id_ == 'physical'
  comps = physicalmodel.components
  assert 'scale' in comps
  assert 'absorption' in comps
  assert 'decay' in comps
  assert list(comps['scale'].parameters) == [1.2, 1.1]
  assert list(comps['decay'].parameters) == [0.1, 0.2]
  assert list(comps['absorption'].parameters) == [0.01, 0.01, 0.01]

  # Test configure reflection table
  rt = flex.reflection_table()
  rt['xyzobs.px.value'] = flex.vec3_double([(0.1, 0.1, 0.1), (0.1, 0.1, 0.1)])
  rt['s1'] = flex.vec3_double([(0.1, 0.1, 1.1), (0.1, 0.1, 1.1)])
  rt['d'] = flex.double([1.0, 1.0])
  exp = generated_single_exp()[0]
  params = generated_param()
  rt = physicalmodel.configure_reflection_table(rt, exp, params)
  assert physicalmodel.components['scale'].col_name in rt
  assert physicalmodel.components['decay'].col_name in rt

  # Test normalise components.
  physicalmodel.components['scale'].update_reflection_data(rt)
  physicalmodel.components['scale'].calculate_scales_and_derivatives()
  physicalmodel.components['decay'].update_reflection_data(rt)
  physicalmodel.components['decay'].calculate_scales_and_derivatives()
  physicalmodel.normalise_components()
  assert list(physicalmodel.components['scale'].inverse_scales) == [1.0, 1.0]
  assert list(physicalmodel.components['decay'].inverse_scales) == [1.0, 1.0]

  # Test from_dict initialisation method.
  physical_dict = {"__id__": "physical", "is_scaled": True, "scale": {
    "n_parameters": 2, "parameters": [0.5, 1.0], "est_standard_devs" : [0.05, 0.1]},
    "configuration_parameters": {"corrections": ["scale"], "s_norm_fac": 0.1,
        "scale_rot_interval": 10.0}}
  physicalmodel = PhysicalScalingModel.from_dict(physical_dict)
  assert physicalmodel.id_ == 'physical'
  assert 'scale' in physicalmodel.components
  assert 'absorption' not in physicalmodel.components
  assert 'decay' not in physicalmodel.components
  assert list(physicalmodel.components['scale'].parameters) == [0.5, 1.0]
  assert list(physicalmodel.components['scale'].parameter_esds) == [0.05, 0.1]

  new_dict = physicalmodel.to_dict()
  assert new_dict == physical_dict

def test_ArrayScalingModel():
  configdict = {'corrections':['decay', 'absorption'], 'n_res_param': 2,
    'n_time_param': 2, 'resmin' : 1.0, 'res_bin_width' : 1.0,
    'time_norm_fac' : 1.0, 'time_rot_interval' : 1.0, 'n_x_param' : 2,
    'n_y_param' : 2, 'xmin' : 0.0, 'ymin' : 0.0, 'x_bin_width' : 1.0,
    'y_bin_width' : 2.0}
  
  parameters_dict = {
      'decay': {'parameters' : flex.double([1.2, 1.1, 1.0, 0.9]),
        'parameter_esds' : None},
      'absorption': {'parameters' : flex.double([0.1, 0.2, 0.1, 0.2, 0.1, 0.2,
        0.1, 0.2]), 'parameter_esds' : None}}

  # Test standard factory initialisation
  arraymodel = ArrayScalingModel(parameters_dict, configdict)
  assert arraymodel.id_ == 'array'
  assert 'decay' in arraymodel.components
  assert 'absorption' in arraymodel.components
  assert 'modulation' not in arraymodel.components
  assert list(arraymodel.components['decay'].parameters) == [1.2, 1.1, 1.0, 0.9]
  assert list(arraymodel.components['absorption'].parameters) == [
    0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2]

  # Test configure reflection table
  rt = flex.reflection_table()
  rt['xyzobs.px.value'] = flex.vec3_double([(0.1, 0.1, 0.1), (0.1, 0.1, 0.1)])
  rt['s1'] = flex.vec3_double([(0.1, 0.1, 1.1), (0.1, 0.1, 1.1)])
  rt['d'] = flex.double([1.0, 1.0])
  exp = generated_single_exp()[0]
  params = generated_param()
  rt = arraymodel.configure_reflection_table(rt, exp, params)
  for comp in ['decay', 'absorption']:
    for col_name in arraymodel.components[comp].col_names:
      assert col_name in rt

  # Test from_dict initialisation method.
  array_dict = {"__id__": "array", "is_scaled": True, "decay": {
    "n_parameters": 4, "parameters": [0.5, 1.0, 0.4, 1.0],
    "est_standard_devs" : [0.05, 0.1, 0.05, 0.1]},
    "configuration_parameters": {"corrections": ["decay"], 'n_res_param': 2,
    'n_time_param': 2, 'resmin' : 1.0, 'res_bin_width' : 1.0,
    'time_norm_fac' : 1.0, 'time_rot_interval' : 1.0}}
  arraymodel = ArrayScalingModel.from_dict(array_dict)
  assert arraymodel.id_ == 'array'
  comps = arraymodel.components
  assert 'modulation' not in comps
  assert 'absorption' not in comps
  assert 'decay' in comps
  assert list(comps['decay'].parameters) == [0.5, 1.0, 0.4, 1.0]
  assert list(comps['decay'].parameter_esds) == [0.05, 0.1, 0.05, 0.1]

  new_dict = arraymodel.to_dict()
  assert new_dict == array_dict
  