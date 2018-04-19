"""
Module of library functions, to perform core scaling operations on reflection
tables and experiment lists. Some functions, such as create_scaling_model and
merging statistics calculations are called from the main dials.scale script,
whereas others are provided as library functions for calling from custom
scripts. The functions defined here should ideally only require reflection
tables and ExperimentList objects (and sometimes phil_scope objects if
necessary), and return common dials objects such as reflection tables and
ExperimentLists.
"""
import pkg_resources
from libtbx import phil
import iotbx.merging_statistics
from cctbx import miller, crystal
from dials.util.options import OptionParser
from dials_scratch.jbe.scaling_code.scaler_factory import SingleScalerFactory,\
  TargetScalerFactory

def scale_against_target(reflection_table, experiment, target_reflection_table,
  target_experiment, params=None, model='KB'):
  """Determine scale factors for a single dataset, by scaling against a target
  reflection table. Requires a single reflection table for the reflections to
  scale and the target dataset, and an ExperimentList for both datasets. The
  params option can also be specified, if None then the default scaling
  configuration is used. The scaling model can be specified individually.

  Returns the reflection table, with added columns 'inverse_scale_factor' and
  'inverse_scale_factor_variance'."""

  assert model in ['physical', 'array', 'KB']
  if not params:
    phil_scope = phil.parse('''
      include scope dials_scratch.jbe.scaling_code.scaling_options.phil_scope
      include scope dials_scratch.jbe.scaling_code.scaling_refiner.scaling_refinery_phil_scope
    ''', process_includes=True)
    optionparser = OptionParser(phil=phil_scope, check_format=False)
    params, _ = optionparser.parse_args(args=None, quick_parse=True)
    params.__inject__('model', model)

  reflections = [reflection_table, target_reflection_table]
  experiment.append(target_experiment[0])
  experiments = create_scaling_model(params, experiment, reflections)
  scaler = TargetScalerFactory.create(params, experiments, reflections,
    is_scaled_list=[False, True])
  scaler.perform_scaling()
  scaler.expand_scales_to_all_reflections(calc_cov=True)
  return scaler.unscaled_scalers[0].reflection_table


def scale_single_dataset(reflection_table, experiment, params=None,
    model='physical'):
  """Determine scale factors for a single dataset. Requires a reflection table
  and an ExperimentList with a single experiment. A custom params option can be
  specified, if not the default scaling params option will be used, with default
  configuration options. The model can be individually specified.

  Returns the reflection table, with added columns 'inverse_scale_factor' and
  'inverse_scale_factor_variance'."""

  assert model in ['physical', 'array']
  if not params:
    phil_scope = phil.parse('''
      include scope dials_scratch.jbe.scaling_code.scaling_options.phil_scope
      include scope dials_scratch.jbe.scaling_code.scaling_refiner.scaling_refinery_phil_scope
    ''', process_includes=True)
    optionparser = OptionParser(phil=phil_scope, check_format=False)
    params, _ = optionparser.parse_args(args=None, quick_parse=True)
    params.__inject__('model', model)

  experiments = create_scaling_model(params, experiment, [reflection_table])
  scaler = SingleScalerFactory.create(params, experiments[0], reflection_table)
  scaler.perform_scaling()
  scaler.perform_scaling(engine=params.scaling_refinery.full_matrix_engine,
    max_iterations=params.scaling_refinery.full_matrix_max_iterations)
  scaler.expand_scales_to_all_reflections(calc_cov=True)
  return scaler.reflection_table

def create_scaling_model(params, experiments, reflections):
  """Create or load a scaling model for multiple datasets."""
  for i, (exp, refl) in enumerate(zip(experiments, reflections)):
    model = experiments.scaling_models()[i]
    if params.scaling_options.target_intensities and i == len(reflections)-1:
      for entry_point in pkg_resources.iter_entry_points('dxtbx.scaling_model_ext'):
        if entry_point.name == 'KB':
          #finds relevant extension in dials.extensions.scaling_model_ext
          factory = entry_point.load().factory()
          exp.scaling_model = factory.create(params, exp, refl)
          exp.scaling_model.set_scaling_model_as_scaled()
    elif model is not None:
      exp.scaling_model = model
    else:
      for entry_point in pkg_resources.iter_entry_points('dxtbx.scaling_model_ext'):
        if entry_point.name == params.model:
          #finds relevant extension in dials.extensions.scaling_model_ext
          factory = entry_point.load().factory()
          exp.scaling_model = factory.create(params, exp, refl)
  return experiments

def calculate_merging_statistics(reflection_table, experiments):
  """Calculate merging statistics for scaled datasets. Datasets are selected
  from the reflection table based on their id, and a list of dataset statistics
  objects and dataset ids are returned."""
  results = []
  ids = []
  dataset_ids = list(set(reflection_table['id']))
  if len(dataset_ids) == 1:
    results.append(calculate_single_merging_stats(reflection_table, experiments[0]))
    ids.append(dataset_ids[0])
  else:
    for dataset_id in dataset_ids:
      refls = reflection_table.select(reflection_table['id'] == dataset_id)
      results.append(calculate_single_merging_stats(refls, experiments[0]))
      ids.append(dataset_id)
  return results, ids

def calculate_single_merging_stats(reflection_table, experiment):
  """Calculate the merging stats for a single dataset."""
  bad_refl_sel = reflection_table.get_flags(
    reflection_table.flags.bad_for_scaling, all=False)
  r_t = reflection_table.select(~bad_refl_sel)
  u_c = experiment.crystal.get_unit_cell().parameters()
  s_g = experiment.crystal.get_space_group()
  miller_set = miller.set(crystal_symmetry=crystal.symmetry(unit_cell=u_c,
    space_group=s_g), indices=r_t['miller_index'], anomalous_flag=False)
  i_obs = miller.array(miller_set, data=r_t['intensity']/
    r_t['inverse_scale_factor'])
  i_obs.set_observation_type_xray_intensity()
  i_obs.set_sigmas((r_t['variance']**0.5)/r_t['inverse_scale_factor'])
  #dataset_id = list(set(reflection_table['id']))[0]
  result = iotbx.merging_statistics.dataset_statistics(
    i_obs=i_obs, n_bins=20, anomalous=False, sigma_filtering=None,
    use_internal_variance=True, eliminate_sys_absent=False)
  return result
