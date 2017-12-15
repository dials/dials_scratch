from __future__ import division
import dials

def read_experiments(filename):
  from dxtbx.model.experiment_list import ExperimentListFactory
  return ExperimentListFactory.from_json_file(filename)

def read_reflections(filename):
  from dials.array_family import flex
  r = flex.reflection_table.from_pickle(filename)
  return r

def process_reference(reference):
  ''' Load the reference spots. '''
  from dials.array_family import flex
  from time import time
  from libtbx.utils import Sorry
  if reference is None:
    return None, None
  st = time()
  assert("miller_index" in reference)
  assert("id" in reference)
  mask = reference.get_flags(reference.flags.indexed)
  rubbish = reference.select(mask == False)
  if mask.count(False) > 0:
    reference.del_selected(mask == False)
  if len(reference) == 0:
    raise Sorry('''
      Invalid input for reference reflections.
      Expected > %d indexed spots, got %d
    ''' % (0, len(reference)))
  mask = reference.get_flags(reference.flags.centroid_outlier)
  if mask.count(True) > 0:
    rubbish.extend(reference.select(mask))
    reference.del_selected(mask)
  mask = reference['miller_index'] == (0, 0, 0)
  if mask.count(True) > 0:
    rubbish.extend(reference.select(mask))
    reference.del_selected(mask)
  mask = reference['id'] < 0
  if mask.count(True) > 0:
    raise Sorry('''
      Invalid input for reference reflections.
      %d reference spots have an invalid experiment id
    ''' % mask.count(True))
  if (reference['panel'] == reference['shoebox'].panels()).count(False) > 0:
    raise RuntimeError('reflection table "panel" column does not match "shoebox" panel')
  return reference

def compute_reference(experiments, reflections, params):
  from dials.algorithms.integration.parallel_integrator import ReferenceCalculatorManager

  
  reference_manager = ReferenceCalculatorManager(
    experiments, 
    reflections, 
    params)
  print reference_manager.summary()

  for task in reference_manager.tasks():
    result = task()
    reference_manager.accumulate(result)
  reference_manager.finalize()

  return reference_manager.result(), reference_manager.reference


if __name__ == '__main__':
  from os.path import join
  import sys
  from dials.util import log
  from dials.array_family import flex


  log.config()

  experiments_filename = sys.argv[1]
  reflections_filename = sys.argv[2]

  experiments = read_experiments(experiments_filename)
  reflections = read_reflections(reflections_filename)
  reflections = process_reference(reflections) 

  from dials.command_line.integrate import phil_scope

  params = phil_scope.extract()
  
  params.integration.block.size = 20
  params.integration.block.units = "frames"
  params.integration.mp.nproc = 8
  #params.integration.mp.nproc = 1
  params.profile.gaussian_rs.fitting.grid_method="circular_grid"
  params.integration.use_dynamic_mask=True


  predicted = flex.reflection_table.from_predictions_multi(
    experiments,
    dmin=params.prediction.d_min,
    dmax=params.prediction.d_max,
    margin=params.prediction.margin,
    force_static=params.prediction.force_static,
    padding=params.prediction.padding)


  matched, reference, unmatched = predicted.match_with_reference(reflections)
    
  predicted.compute_bbox(experiments)

  print "Dynamic Mask: ", experiments[0].imageset.has_dynamic_mask()

  print "Read %d reflections" % len(reflections)


  from time import time
  st = time()
  
  reflections, reference = compute_reference(experiments, predicted, params)
  print "Num used in modelling", reflections.get_flags(reflections.flags.used_in_modelling).count(True)
  print "Time taken: ", time() - st

  import pickle
  pickle.dump(reference, open("reference_profiles.pickle", "w"))
