#!/usr/bin/env cctbx.python

#
#  Copyright (C) (2015) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

"""
Development testing of restraints_parameterisation and associated classes

"""

# Python and cctbx imports
from __future__ import division
import os
from libtbx.phil import parse
import libtbx.load_env # required for libtbx.env.find_in_repositories
from dxtbx.model.experiment.experiment_list import ExperimentListFactory
from dials.algorithms.refinement import RefinerFactory
from dials.array_family import flex

from dials.algorithms.refinement.restraints import RestraintsParameterisation

# The phil scope
from dials.algorithms.refinement.refiner import phil_scope
user_phil = parse('''
refinement
{
  parameterisation
  {
    crystal
    {
      unit_cell
      {
        restraints
        {
          tie_to_target
          {
            values=10,10,10,90,90,90
            sigmas=1,1,1,1,1,1
            id=0
          }
        }
      }
    }
  }
}
''')

working_phil = phil_scope.fetch(source=user_phil)
working_params = working_phil.extract()

def test1(params):

  dials_regression = libtbx.env.find_in_repositories(
    relative_path="dials_regression",
    test=os.path.isdir)

  # use the multi stills test data
  data_dir = os.path.join(dials_regression, "refinement_test_data", "multi_stills")
  experiments_path = os.path.join(data_dir, "combined_experiments.json")
  pickle_path = os.path.join(data_dir, "combined_reflections.pickle")

  experiments = ExperimentListFactory.from_json_file(experiments_path,
                check_format=False)
  reflections = flex.reflection_table.from_pickle(pickle_path)

  refiner = RefinerFactory.from_parameters_data_experiments(params,
        reflections, experiments)

  # hack to extract the restraints parameterisation from the Refiner
  rp = refiner._target._restraints_parameterisation

  print "do something"
  return

def test2():
  #from libtbx.test_utils import approx_equal
  #from scitbx.array_family import flex
  from math import pi
  from random import gauss
  from dials.test.algorithms.refinement.setup_geometry import Extract

  from dxtbx.model.experiment.experiment_list import ExperimentList, Experiment
  #from dials.algorithms.refinement.prediction import ScansRayPredictor, \
  #  ExperimentsPredictor

  #### Import model parameterisations

  from dials.algorithms.refinement.parameterisation.prediction_parameters import \
      XYPhiPredictionParameterisation
  from dials.algorithms.refinement.parameterisation.detector_parameters import \
      DetectorParameterisationSinglePanel
  from dials.algorithms.refinement.parameterisation.beam_parameters import \
      BeamParameterisation
  from dials.algorithms.refinement.parameterisation.crystal_parameters import \
      CrystalOrientationParameterisation, \
      CrystalUnitCellParameterisation

  overrides = """geometry.parameters.crystal.a.length.range = 10 50
  geometry.parameters.crystal.b.length.range = 10 50
  geometry.parameters.crystal.c.length.range = 10 50"""

  master_phil = parse("""
      include scope dials.test.algorithms.refinement.geometry_phil
      """, process_includes=True)

  models = Extract(master_phil, overrides)

  mydetector = models.detector
  mygonio = models.goniometer
  mycrystal = models.crystal
  mybeam = models.beam

  # Build a mock scan for a 72 degree sweep
  sweep_range = (0., pi/5.)
  from dxtbx.model.scan import scan_factory
  sf = scan_factory()
  myscan = sf.make_scan(image_range = (1,720),
                        exposure_times = 0.1,
                        oscillation = (0, 0.1),
                        epochs = range(720),
                        deg = True)

  # Create parameterisations of these models

  det_param = DetectorParameterisationSinglePanel(mydetector)
  s0_param = BeamParameterisation(mybeam, mygonio)
  xlo_param = CrystalOrientationParameterisation(mycrystal)
  xluc_param = CrystalUnitCellParameterisation(mycrystal)

  # Create an ExperimentList
  experiments = ExperimentList()
  experiments.append(Experiment(
        beam=mybeam, detector=mydetector, goniometer=mygonio, scan=myscan,
        crystal=mycrystal, imageset=None))

  # Build a prediction parameterisation
  pred_param = XYPhiPredictionParameterisation(experiments,
                 detector_parameterisations = [det_param],
                 beam_parameterisations = [s0_param],
                 xl_orientation_parameterisations = [xlo_param],
                 xl_unit_cell_parameterisations = [xluc_param])

  # Build a restraints parameterisation
  rp = RestraintsParameterisation(detector_parameterisations = [det_param],
               beam_parameterisations = [s0_param],
               xl_orientation_parameterisations = [xlo_param],
               xl_unit_cell_parameterisations = [xluc_param])

  # make a unit cell target
  sigma = 1.
  uc = mycrystal.get_unit_cell().parameters()
  print uc
  target_uc = [gauss(e, sigma) for e in uc]
  print target_uc

  rp.add_restraints_to_target_xl_unit_cell(experiment_id=0, values=target_uc,
                                           sigma=[sigma]*6)

  # get analytical values and gradients
  vals, grads = rp.get_values_and_gradients()

  # for debugging, print finite difference gradients of the unit cell parameters
  # versus

  # get finite difference gradients
  p_vals = pred_param.get_param_vals()
  deltas = [1.e-7] * len(p_vals)

  fd_grad=[]
  for i in range(len(deltas)):

    val = p_vals[i]

    p_vals[i] -= deltas[i] / 2.
    pred_param.set_param_vals(p_vals)

    rev_state, _ = rp.get_values_and_gradients()
    rev_state = flex.double(rev_state)

    p_vals[i] += deltas[i]
    pred_param.set_param_vals(p_vals)

    fwd_state, _ = rp.get_values_and_gradients()
    fwd_state = flex.double(fwd_state)

    p_vals[i] = val

    fd = (fwd_state - rev_state) / deltas[i]
    fd_grad.append(fd)

  # for comparison, fd_grad is a list of flex.doubles, each of which corresponds
  # to a column of the sparse matrix grads.
  for i, fd in enumerate(fd_grad):
    # extract dense column from the sparse matrix
    an = grads.col(i).as_dense_vector()

    print list(an.round(6))
    print list(fd.round(6))
    print


  # enter interactive console
  from dials.util.command_line import interactive_console; interactive_console()

if __name__ == '__main__':

  test1(working_params)
  test2()
