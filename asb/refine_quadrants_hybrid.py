#!/usr/bin/env dials.python
from __future__ import division
import sys

from libtbx.phil import command_line, parse
from dxtbx.serialize import load as load_dxtbx
from dxtbx.model.experiment.experiment_list import ExperimentList, Experiment

from dials.model.serialize import load as load_dials
from dials.array_family import flex
from dials.algorithms.refinement import RefinerFactory

"""This script is intended for quadrant refinement using still shot data
collected on a CSPAD detector. This version of the script uses a 'hybrid
minimiser'. Rather than a single joint refinement job of all crystals and the
detector, only the detector parameters are refined at first (using all data)
then each crystal is refined individually. This forms one macrocycle."""

def load_input(exp_path, ref_path):

  refs = load_dials.reflections(ref_path)
  exp = load_dxtbx.experiment_list(exp_path , check_format=False)[0]
  return refs, exp

class ExperimentFromCrystal(object):

  def __init__(self, reference_beam, reference_detector):

    self.reference_beam = reference_beam
    self.reference_detector = reference_detector
    return

  def __call__(self, crystal):

    return Experiment(beam=self.reference_beam,
                      detector=self.reference_detector,
                      crystal=crystal)

class DetectorRefiner(object):

  user_phil=parse("""
  refinement{
    parameterisation {
      beam.fix=all
      crystal.fix=all
      detector.hierarchy_level=1
      sparse=True
    }
    target.jacobian_max_nref=100000
    reflections.do_outlier_rejection=True
    reflections.minimum_number_of_reflections=1
  }""")
  from dials.data.refinement import phil_scope as refinement_phil
  working_phil = refinement_phil.fetch(sources=[user_phil])

  def __call__(self, experiments, reflections):

    self.working_phil.show()
    params = self.working_phil.extract()
    refiner = RefinerFactory.from_parameters_data_experiments(
        params, reflections, experiments, verbosity=1)
    refiner.run()
    experiments = refiner.get_experiments()

    return experiments

class CrystalRefiners(object):

  user_phil=parse("""
  refinement{
    parameterisation {
      beam.fix=all
      detector.fix=all
    }
    reflections.do_outlier_rejection=True
    reflections.minimum_number_of_reflections=1
  }""")
  from dials.data.refinement import phil_scope as refinement_phil
  working_phil = refinement_phil.fetch(sources=[user_phil])

  def __call__(self, experiments, reflections):

    self.working_phil.show()
    params = self.working_phil.extract()

    for iexp, exp in enumerate(experiments):

      print "Refining crystal", iexp
      # reflection subset for a single experiment
      refs = reflections.select(reflections['id'] == iexp)
      refs['id'] = flex.int(len(refs),0)
      # experiment list for a single experiment
      exps=ExperimentList()
      exps.append(exp)
      refiner = RefinerFactory.from_parameters_data_experiments(
        params, refs, exps, verbosity=1)
      # do refinement
      refiner.run()
      refined_exps = refiner.get_experiments()
      # replace this experiment with the refined one
      experiments[iexp] = refined_exps[0]

    return experiments

if __name__ =="__main__":

  if len(sys.argv) != 2: exit("please pass the path to a phil file")
  phil = sys.argv[1]

  master_phil = parse("""
    input
      .multiple = True
    {
      experiments = None
        .type = path
      reflections = None
        .type = path
    }
    n_macrocycles = 1
      .type = int(value_min=1)
    """)

  cmd_line = command_line.argument_interpreter(master_params=master_phil)
  working_phil = cmd_line.process_and_fetch(args=(phil,))
  working_params = working_phil.extract()

  for input in working_params.input:
    print input.experiments, input.reflections

  assert len(working_params.input) > 1
  print len(working_params.input), "datasets specified as input"

  e = enumerate(working_params.input)
  i, line = e.next()
  reflections, exp = load_input(line.experiments, line.reflections)
  assert reflections['id'].all_eq(0)
  experiment_from_crystal=ExperimentFromCrystal(exp.beam, exp.detector)

  experiments=ExperimentList()
  experiments.append(experiment_from_crystal(exp.crystal))

  for i, line in e:
    refs, exp = load_input(line.experiments, line.reflections)
    refs['id'] = flex.int(len(refs),i)
    reflections.extend(refs)
    experiments.append(experiment_from_crystal(exp.crystal))

  dr = DetectorRefiner()
  cr = CrystalRefiners()

  for cycle in range(working_params.n_macrocycles):

    print "MACROCYCLE %02d" % (cycle + 1)
    print "=============\n"
    # first run: multi experiment joint refinement of detector with fixed beam and
    # crystals
    experiments = dr(experiments, reflections)

    # second run
    experiments = cr(experiments, reflections)

  # save the refined experiments
  from dxtbx.model.experiment.experiment_list import ExperimentListDumper
  dump = ExperimentListDumper(experiments)
  experiments_filename = "refined_experiments.json"
  dump.as_json(experiments_filename)
  print "refined geometry written to {0}".format(experiments_filename)
