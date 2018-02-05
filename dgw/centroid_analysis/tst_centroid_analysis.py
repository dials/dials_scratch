#!/usr/bin/env cctbx.python

#
#  Copyright (C) (2017) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

"""
Tests the CentroidAnalyser class used in refinement

"""

from __future__ import absolute_import, division
import os
import libtbx.load_env # required for libtbx.env.find_in_repositories
from libtbx.test_utils import approx_equal
from dials.array_family import flex
from dials.algorithms.refinement.analysis.centroid_analysis import \
  CentroidAnalyser

def test1():
  """Test centroid analysis on an indexed.pickle"""

  if not libtbx.env.has_module("dials_regression"):
    print "Skipping test1 in " + __file__ + " as dials_regression not present"
    return

  dials_regression = libtbx.env.find_in_repositories(
    relative_path="dials_regression",
    test=os.path.isdir)

  # Use dials_regression/indexing_test_data/i04_weak_data/indexed.pickle for
  # this test. This file is somewhat malformed, in that the observed and
  # calculated phi centroids are identical. We should not produce a phi
  # periodogram in that case
  data_dir = os.path.join(dials_regression, "indexing_test_data", "i04_weak_data")
  pickle_path = os.path.join(data_dir, "indexed.pickle")

  # load the reflections and do analysis
  rt = flex.reflection_table.from_pickle(pickle_path)
  results = CentroidAnalyser(rt)()

  #### Check that results are as we expect ####

  # Only one experiment id
  assert len(results) == 1
  results = results[0]

  # Residuals per block
  assert results['nblocks'] == 80
  assert approx_equal(results['block_size'], 1.010557248)
  expected = [0.335558, 0.335558, 0.329863, 0.331119, 0.330501, 0.331638,
    0.342241, 0.326378, 0.333448, 0.336584, 0.326387, 0.331208, 0.324562,
    0.333608, 0.331141, 0.325593, 0.336456, 0.333159, 0.323751, 0.340046,
    0.327231, 0.321859, 0.331268, 0.331076, 0.320894, 0.323385, 0.326196,
    0.320985, 0.327399, 0.324754, 0.324914, 0.329531, 0.328205, 0.330763,
    0.324561, 0.32797, 0.331537, 0.32979, 0.332977, 0.334632, 0.324989,
    0.32758, 0.333216, 0.339781, 0.311982, 0.327621, 0.34168, 0.322333,
    0.321525, 0.3229, 0.324373, 0.329953, 0.324727, 0.320513, 0.332656,
    0.319292, 0.322898, 0.316097, 0.319094, 0.33926, 0.316256, 0.335599,
    0.319271, 0.321485, 0.332974, 0.31678, 0.333504, 0.323779, 0.326404,
    0.334334, 0.334769, 0.329136, 0.327918, 0.334443, 0.333598, 0.316668,
    0.317806, 0.319651, 0.326678, 0.326678]
  assert approx_equal(results['av_x_resid_per_block'], expected, eps=5.e-6)
  expected = [0.140451, 0.140451, 0.134973, 0.133594, 0.138062, 0.13792,
    0.142743, 0.135466, 0.139534, 0.140333, 0.136268, 0.138599, 0.142424,
    0.137319, 0.139688, 0.14409, 0.140079, 0.146143, 0.135093, 0.134561,
    0.141375, 0.142994, 0.143815, 0.129358, 0.151952, 0.136534, 0.137577,
    0.14513, 0.136317, 0.142892, 0.138412, 0.142925, 0.131937, 0.142954,
    0.144306, 0.134972, 0.143383, 0.129731, 0.139279, 0.13215, 0.137761,
    0.142486, 0.135102, 0.122229, 0.155227, 0.134371, 0.129951, 0.151823,
    0.130035, 0.143184, 0.143669, 0.139766, 0.143137, 0.141744, 0.134653,
    0.144242, 0.14377, 0.152644, 0.140361, 0.13361, 0.137553, 0.133924,
    0.149607, 0.155575, 0.149489, 0.146792, 0.14263, 0.148653, 0.139416,
    0.132103, 0.131738, 0.13657, 0.144084, 0.140233, 0.137387, 0.133706,
    0.138388, 0.14073, 0.138223, 0.138223]
  assert approx_equal(results['av_y_resid_per_block'], expected, eps=5.e-6)

  # For this data no interval widths have been suggested
  assert results['x_interval'] is None
  assert results['y_interval'] is None
  assert results['phi_interval'] is None

  # It isn't necessary to test the content of the periodograms here, as this is
  # tested in scitbx/math. However, do check that the phi periodogram has been
  # skipped whilst the others exist
  assert results['x_periodogram'] is not None
  assert results['y_periodogram'] is not None
  assert results['phi_periodogram'] is None

  print "OK"
  return


def test2():
  """Test on simulated data"""

  # Get models for reflection prediction
  import dials.test.algorithms.refinement.setup_geometry as setup_geometry

  from libtbx.phil import parse
  overrides = """geometry.parameters.crystal.a.length.value = 77
  geometry.parameters.crystal.b.length.value = 77
  geometry.parameters.crystal.c.length.value = 37"""

  master_phil = parse("""
      include scope dials.test.algorithms.refinement.geometry_phil
      """, process_includes=True)

  from dxtbx.model import Crystal
  models = setup_geometry.Extract(master_phil)
  crystal = Crystal(
    real_space_a=(2.62783398111729, -63.387215823567125, -45.751375737456975),
    real_space_b=(15.246640559660356, -44.48254330406616, 62.50501032727026),
    real_space_c=(-76.67246874451074, -11.01804131886244, 10.861322446352226),
    space_group_symbol="I 2 3")
  detector = models.detector
  goniometer = models.goniometer
  beam = models.beam

  # Build a mock scan for a 180 degree sweep
  from dxtbx.model import ScanFactory
  sf = ScanFactory()
  scan = sf.make_scan(image_range = (1,1800),
                        exposure_times = 0.1,
                        oscillation = (0, 0.1),
                        epochs = range(1800),
                        deg = True)

  # Build an experiment list
  from dxtbx.model.experiment_list import ExperimentList, Experiment
  experiments = ExperimentList()
  experiments.append(Experiment(
        beam=beam, detector=detector, goniometer=goniometer,
        scan=scan, crystal=crystal, imageset=None))

  # Generate all indices in a 1.5 Angstrom sphere
  from dials.algorithms.spot_prediction import IndexGenerator
  from cctbx.sgtbx import space_group, space_group_symbols
  resolution = 1.5
  index_generator = IndexGenerator(crystal.get_unit_cell(),
                  space_group(space_group_symbols(1).hall()).type(), resolution)
  indices = index_generator.to_array()

  # Predict rays within the sweep range
  from dials.algorithms.refinement.prediction import ScansRayPredictor
  sweep_range = scan.get_oscillation_range(deg=False)
  ray_predictor = ScansRayPredictor(experiments, sweep_range)
  obs_refs = ray_predictor(indices)

  # Take only those rays that intersect the detector
  from dials.algorithms.spot_prediction import ray_intersection
  intersects = ray_intersection(detector, obs_refs)
  obs_refs = obs_refs.select(intersects)

  # Make a reflection predictor and re-predict for all these reflections. The
  # result is the same, but we gain also the flags and xyzcal.px columns
  from dials.algorithms.refinement.prediction import ExperimentsPredictor
  ref_predictor = ExperimentsPredictor(experiments)
  obs_refs['id'] = flex.int(len(obs_refs), 0)
  obs_refs = ref_predictor(obs_refs)

  # Copy 'observed' centroids from the predicted ones, applying sinusoidal
  # offsets
  obs_x, obs_y, obs_z = obs_refs['xyzcal.mm'].parts()

  # obs_z is in range (0, pi). Calculate offsets for phi at twice that
  # frequency
  im_width = scan.get_oscillation(deg=False)[1]
  z_off = flex.sin(2 * obs_z) * im_width
  obs_z += z_off

  # Calculate offsets for x
  pixel_size = detector[0].get_pixel_size()
  x_off = flex.sin(20 * obs_z) * pixel_size[0]

  # Calculate offsets for y with a phase-shifted sine wave
  from math import pi
  y_off = flex.sin(4 * obs_z + pi/6) * pixel_size[1]

  # Incorporate the offsets into the 'observed' centroids
  obs_z += z_off
  obs_x += x_off
  obs_y += y_off
  obs_refs['xyzobs.mm.value'] = flex.vec3_double(obs_x, obs_y, obs_z)

  # Now do centroid analysis of the residuals
  results = CentroidAnalyser(obs_refs, debug=True)()

  # FIXME this test shows that the suggested interval width heuristic is not
  # yet robust. This simulation function seems a useful direction to proceed
  # in though
  raise RuntimeError('test2 failed')

  print "OK"
  return

if __name__ == '__main__':

  # simple test on an indexed.pickle
  test1()

  # A test using simulated data with sinusodial residuals
  test2()

