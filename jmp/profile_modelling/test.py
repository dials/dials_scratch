from __future__ import division
from dials_scratch.jmp.profile_modelling import ReciprocalLatticePointSpread
from dials_scratch.jmp.profile_modelling import WavelengthSpread
from dials_scratch.jmp.profile_modelling import MosaicBlockAngularSpread
from dials_scratch.jmp.profile_modelling import CovarianceMatrix
from dials_scratch.jmp.profile_modelling import Model3D
from dials.array_family import flex
from scitbx import matrix

def tst_reciprocal_lattice_point_spread():
  A = matrix.sqr((
    1.0, 2.0, 3.0,
    3.0, 1.0, 2.0,
    2.0, 3.0, 1.0))

  parameters = (
    0.1,
    0.001, 0.2,
    0.002, 0.003, 0.3)

  spread = ReciprocalLatticePointSpread(A.elems, parameters)

  cov = spread.covariance

  spread.covariance = cov

  assert all(abs(p1-p2) < 1e-7 for p1, p2 in zip(spread.parameters, parameters))

  print 'OK'

def tst_mosaic_block_angular_spread():

  r = (0.1, 0.2, -0.3)
  spread = MosaicBlockAngularSpread(r, 0.001)

  cov = spread.covariance

  spread.covariance = cov

  assert abs(spread.parameter - 0.001) < 1e-7

  print 'OK'

def tst_wavelength_spread():

  s0 = (0.2, -0.1, 1)
  r = (0.1, 0.2, -0.3)
  spread = WavelengthSpread(s0, r, 0.01)

  cov = spread.covariance

  spread.covariance = cov

  assert abs(spread.parameter - 0.01) < 1e-7

  print 'OK'

def tst_covariance_matrix():

  A = matrix.sqr((
    1.0, 2.0, 3.0,
    3.0, 1.0, 2.0,
    2.0, 3.0, 1.0))

  s0 = (0.2, -0.1, 1)
  r = (0.1, 0.2, -0.3)

  cov = CovarianceMatrix(A, s0, r, use_wavelength_spread=False)

  assert cov.use_mosaic_block_angular_spread == True
  assert cov.use_wavelength_spread == False
  assert cov.num_parameters() == 7

  parameters = (
    0.1,
    0.001, 0.2,
    0.002, 0.003, 0.3)

  cov.compose(
    reciprocal_lattice_point_spread=parameters,
    mosaic_block_angular_spread=0.01)

  assert all(abs(p1-p2) < 1e-7 for p1, p2 in
             zip(cov.reciprocal_lattice_point_spread().parameters, parameters))
  assert abs(cov.mosaic_block_angular_spread().parameter - 0.01) < 1e-7

  cov = CovarianceMatrix(A, s0, r)

  cov.parameters = flex.double((0.1, 0.001, 0.2, 0.002, 0.004, 0.002, 0.002, 0.0003))

  sigma = matrix.sqr(cov.covariance)

  from numpy.linalg import eigvals
  eigen_values = eigvals(sigma.as_list_of_lists())
  assert all(e > 0 for e in eigen_values)

  print 'OK'

def tst_model_3d():

  from dials.array_family import flex
  from dxtbx.model.experiment_list import ExperimentListFactory
  experiments_filename = '/home/upc86896/Data/bag_training/processed_profile/profile_model/experiments.json'
  reflections_filename = '/home/upc86896/Data/bag_training/processed_profile/profile_model/reflections.pickle'

  experiments = ExperimentListFactory.from_json_file(experiments_filename)
  reflections = flex.reflection_table.from_pickle(reflections_filename)

  selection = reflections.get_flags(reflections.flags.used_in_refinement)
  reflections = reflections.select(selection)
  reflections.sort("intensity.sum.value", reverse=True)
  reflections = reflections[0:10]


  for refl in reflections:


    beam = experiments[0].beam
    panel = experiments[0].detector[0]
    goniometer = experiments[0].goniometer
    scan = experiments[0].scan
    A = experiments[0].crystal.get_A()
    s0 = experiments[0].beam.get_s0()

    h = matrix.col(refl['miller_index'])
    r0 = matrix.sqr(A) * h
    x, y, z = refl['xyzcal.px']

    cov = CovarianceMatrix(A, s0, r0,
                           use_wavelength_spread=False,
                           use_mosaic_block_angular_spread=False)
    cov.compose(reciprocal_lattice_point_spread=(
      0.001,
      0.000, 0.001,
      0.000, 0.000, 0.001))
    sigma = cov.covariance


    model = Model3D(beam, panel, goniometer, scan, sigma, r0)

    r1 = model.coord(x, y, z)
    f = model.f(x, y, z)

    d = matrix.col(r1) - r0
    assert d.length() < 1e-3
    assert f > 0

  print 'OK'



def run():
  tst_reciprocal_lattice_point_spread()
  tst_mosaic_block_angular_spread()
  tst_wavelength_spread()
  tst_covariance_matrix()
  tst_model_3d()

if __name__ == '__main__':
  run()
