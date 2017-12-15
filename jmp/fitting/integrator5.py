from __future__ import division
import dials

def read_experiments(filename):
  from dxtbx.model.experiment_list import ExperimentListFactory
  return ExperimentListFactory.from_json_file(filename)

def read_reference(filename):
  import cPickle as pickle
  return pickle.load(open(filename))

def read_reflections(filename):
  from dials.array_family import flex
  r = flex.reflection_table.from_pickle(filename)
  del r['shoebox']
  return r

def construct_reference(experiments, reference):
    from dials.algorithms.integration.parallel_integrator import GaussianRSReferenceProfileData
    from dials.algorithms.integration.parallel_integrator import GaussianRSMultiCrystalReferenceProfileData
    from dials.algorithms.integration.parallel_integrator import ReferenceProfileData
    from dials.algorithms.profile_model.modeller import CircleSampler
    from dials.array_family import flex
    from dials.algorithms.profile_model.gaussian_rs.transform import TransformSpec

    assert len(reference) % 9 == 0
    num_scan_points = len(reference) // 9
    
    data_spec = GaussianRSMultiCrystalReferenceProfileData()
    for e in experiments:

      sampler = CircleSampler(
        e.detector[0].get_image_size(),
        e.scan.get_array_range(),
        num_scan_points)


      spec = TransformSpec(
        e.beam,
        e.detector,
        e.goniometer,
        e.scan,
        e.profile.sigma_b(deg=False),
        e.profile.sigma_m(deg=False),
        e.profile.n_sigma() * 1.5,
        grid_size)

      temp = reference

      reference = ReferenceProfileData()
      for d, m in temp:
        reference.append(d, m)

      spec = GaussianRSReferenceProfileData(reference, sampler, spec)

      data_spec.append(spec)
    return data_spec


def integrate(experiments, reflections, reference, params):
  from dials.algorithms.integration.parallel_integrator import IntegrationManager

  
  integrator_manager = IntegrationManager(
    experiments, 
    reflections, 
    reference,
    params)
  print integrator_manager.summary()

  for task in integrator_manager.tasks():
    result = task()
    integrator_manager.accumulate(result)
  integrator_manager.finalize()

  return integrator_manager.result()


if __name__ == '__main__':
  from os.path import join
  import sys
  from dials.util import log


  log.config()

  experiments_filename = sys.argv[1]
  reference_filename = sys.argv[2]
  grid_size = int(sys.argv[3])
  if sys.argv[4] == "True":
    detector_space = True
  elif sys.argv[4] == "False":
    detector_space = False
  else:
    raise RuntimeError(sys.argv[4])
  reflections_filename = sys.argv[5]

  experiments = read_experiments(experiments_filename)
  reflections = read_reflections(reflections_filename)
  reference = read_reference(reference_filename)

  reference = construct_reference(experiments, reference[0])

  print "Dynamic Mask: ", experiments[0].imageset.has_dynamic_mask()

  print "Read %d reflections" % len(reflections)

  detector_space = True
  deconvolution = False

  from dials.command_line.integrate import phil_scope

  params = phil_scope.extract()

  params.integration.mp.nproc = 8
  
  from time import time
  st = time()

  from dials.array_family import flex 
  reflections["intensity.prf_old.value"] = reflections["intensity.prf.value"]
  reflections["intensity.prf_old.variance"] = reflections["intensity.prf.variance"]
  reflections["intensity.prf.value"] = flex.double(len(reflections))
  reflections["intensity.prf.variance"] = flex.double(len(reflections))
  
  reflections = integrate(experiments, reflections, reference, params)
  print "Num profile fitted", reflections.get_flags(reflections.flags.integrated_prf).count(True)
  print "Time taken: ", time() - st

  print list(reflections.keys())

  reflections.as_pickle("integrated.pickle")
