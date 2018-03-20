#!/usr/bin/env python
#
# dials.potato.py
#
#  Copyright (C) 2018 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import absolute_import, division
from libtbx.phil import parse
from dials.algorithms.profile_model.gaussian_rs.calculator import ComputeEsdBeamDivergence
from dials.algorithms.profile_model.gaussian_rs import MaskCalculator
from dials.algorithms.profile_model.gaussian_rs import BBoxCalculator
from dials.algorithms.profile_model.gaussian_rs import CoordinateSystem2d
from dials_scratch.jmp.potato.profile_refiner import ProfileRefinerData
from dials_scratch.jmp.potato.profile_refiner import ProfileRefiner
from dials_scratch.jmp.potato.profile_refiner import print_eigen_values_and_vectors
from dials_scratch.jmp.potato.crystal_refiner import CrystalRefiner
from dials_scratch.jmp.potato.parameterisation import SimpleMosaicityParameterisation
from dials_scratch.jmp.potato.model import compute_change_of_basis_operation
from dials_scratch.jmp.potato.model import SimpleMosaicityModel
from dials_scratch.jmp.potato import chisq_quantile
from dials_scratch.jmp.potato import Predictor
from dials_scratch.jmp.potato import BBoxCalculator as BBoxCalculatorNew
from dials_scratch.jmp.potato import MaskCalculator as MaskCalculatorNew
from dials.algorithms.spot_prediction import IndexGenerator
from scitbx.linalg import eigensystem
from dials.array_family import flex
from scitbx import matrix
from math import pi, sqrt, floor, ceil, exp
from dials.algorithms.shoebox import MaskCode
import logging

logger = logging.getLogger(__name__)

# Parameters
phil_scope = parse('''

  refinement {

    n_macro_cycles=3
      .type = int

  }

  prediction {
    d_min = None
      .type = float

    probability = 0.997
      .type = float
  }

  debug {
    output {
      shoeboxes = True
        .type = bool

      plots = False
        .type = bool

      print_shoeboxes = False
        .type = bool
    }
  }


''')


class InitialIntegrator(object):
  '''
  A class to do an initial integration of strong spots

  '''

  def __init__(self, params, experiments, reflections):
    '''
    Do the initial integration

    '''

    # Save the experiments and reflections
    self.experiments = experiments
    self.reflections = reflections

    # Do the processing
    self._update_observed_reflection_predictions()
    self._compute_sigma_d()
    self._compute_bbox()
    self._allocate_shoebox()
    self._extract_shoebox()
    self._compute_mask()
    self._compute_background()
    self._compute_intensity()

    # Print shoeboxes
    if params.debug.output.print_shoeboxes:
      self._print_shoeboxes()

  def _update_observed_reflection_predictions(self):
    '''
    Make sure we have the correct reciprocal lattice vector

    '''
    print "Updating predictions for %d reflections" % len(self.reflections)

    # Get stuff from experiment
    A = matrix.sqr(self.experiments[0].crystal.get_A())
    s0 = matrix.col(self.experiments[0].beam.get_s0())

    # Compute the vector to the reciprocal lattice point
    # since this is not on the ewald sphere, lets call it s2
    h = self.reflections['miller_index']
    s1 = flex.vec3_double(len(h))
    s2 = flex.vec3_double(len(h))
    for i in range(len(self.reflections)):
      r = A*matrix.col(h[i])
      s2[i] = s0 + r
      s1[i] = matrix.col(s2[i]).normalize()*s0.length()
    self.reflections['s1'] = s1
    self.reflections['s2'] = s2

    # Compute the ray intersections
    xyzpx = flex.vec3_double()
    xyzmm = flex.vec3_double()
    for i in range(len(s2)):
      ss = s1[i]
      mm = self.experiments[0].detector[0].get_ray_intersection(ss)
      px = self.experiments[0].detector[0].millimeter_to_pixel(mm)
      xyzpx.append(px + (0,))
      xyzmm.append(mm + (0,))
    self.reflections['xyzcal.mm'] = xyzmm
    self.reflections['xyzcal.px'] = xyzpx

  def _compute_sigma_d(self):
    '''
    Compute and initial spot size estimate

    '''
    print "Computing initial sigma d estimate for %d reflections" % len(self.reflections)
    compute_sigma_d = ComputeEsdBeamDivergence(
      self.experiments[0].detector,
      self.reflections)
    self.sigma_d = compute_sigma_d.sigma()
    print "Sigma D: %.5f degrees" % (self.sigma_d * 180 / pi)
    print ""

  def _compute_bbox(self):
    '''
    Compute the bounding box

    '''

    print "Computing the bounding box for %d reflections" % len(self.reflections)

    # Initialise the bounding box calculator
    compute_bbox = BBoxCalculator(
      self.experiments[0].crystal,
      self.experiments[0].beam,
      self.experiments[0].detector,
      self.experiments[0].goniometer,
      self.experiments[0].scan,
      self.sigma_d * 6,
      0)

    # Compute the bounding box
    bbox = compute_bbox(
      self.reflections['s1'],
      self.reflections['xyzcal.px'].parts()[2],
      self.reflections['panel'])

    # Set in the reflection table
    self.reflections['bbox'] = bbox

  def _allocate_shoebox(self):
    '''
    Allocate the shoebox

    '''
    self.reflections['shoebox'] = flex.shoebox(
      self.reflections['panel'],
      self.reflections['bbox'],
      allocate=True)

  def _compute_mask(self):
    '''
    Compute the spot mask

    '''
    print "Creating the foreground mask for %d reflections" % len(self.reflections)

    # Initialise the mask calculator
    mask_foreground = MaskCalculator(
      self.experiments[0].crystal,
      self.experiments[0].beam,
      self.experiments[0].detector,
      self.experiments[0].goniometer,
      self.experiments[0].scan,
      self.sigma_d * 3,
      0)

    # Compute the reflection mask
    mask_foreground(
      self.reflections['shoebox'],
      self.reflections['s1'],
      self.reflections['xyzcal.px'].parts()[2],
      self.reflections['panel'])

  def _extract_shoebox(self):
    '''
    Extract the shoebox

    '''
    print "Extracting shoebox from image for %d reflections" % len(self.reflections)
    self.reflections.extract_shoeboxes(self.experiments[0].imageset)

  def _compute_background(self):
    '''
    Compute the reflection background

    '''
    print "Computing background for %d reflections" % len(self.reflections)
    self.reflections.compute_background(self.experiments)

  def _compute_intensity(self):
    '''
    Compute the reflection intensity

    '''
    print "Computing intensity for %d reflections" % len(self.reflections)
    self.reflections.compute_summed_intensity()
    print "%d reflections integrated" % self.reflections.get_flags(
      self.reflections.flags.integrated_sum).count(True)

  def _print_shoeboxes(self):
    '''
    Print the shoeboxes

    '''
    sbox = self.reflections['shoebox']
    for r in range(len(sbox)):
      data = sbox[r].data
      mask = sbox[r].mask
      print mask.as_numpy_array()
      print data.as_numpy_array()


class Refiner(object):
  '''
  A class to do the refinement

  '''

  def __init__(self, experiments, reflections, sigma_d, params):
    '''
    Initialise the refiner

    '''

    # Save the experiments and reflections
    self.params = params
    self.experiments = experiments
    self.reflections = reflections
    self.sigma_d = sigma_d
    self._profile_parameters = None

    # Preprocess the reflections
    self._preprocess()

    # Do the macro cycles of refinement between refining the profile parameters
    # and refining the crystal orientation and unit cell
    for cycle in range(self.params.refinement.n_macro_cycles):
      print ""
      print "Macro cycle %d" % (cycle+1)
      self._refine_profile()
      self._refine_crystal()

    # Create the profile model
    self.profile_model = SimpleMosaicityParameterisation(self._profile_parameters)


  def _preprocess(self):
    '''
    Preprocess the reflections

    '''

    # Filter based on centroid distance
    self._filter_reflections_based_on_centroid_distance()

    # Make some plots
    if self.params.debug.output.plots:
      self._plot_distance_from_ewald_sphere()

    # Construct the profile refiner data
    self._refiner_data = ProfileRefinerData.from_reflections(
      self.experiments[0],
      self.reflections)

  def _filter_reflections_based_on_centroid_distance(self):
    '''
    Filter reflections too far from predicted position

    '''
    Xobs, Yobs, _ = self.reflections['xyzobs.px.value'].parts()
    Xcal, Ycal, _ = self.reflections['xyzcal.px'].parts()
    D = flex.sqrt((Xobs-Xcal)**2 + (Yobs-Ycal)**2)
    selection = D < 2
    self.reflections = self.reflections.select(selection)
    print "Selected %d reflections with centroid-prediction distance < 2px" % len(self.reflections)

  def _refine_profile(self):
    '''
    Do the profile refinement

    '''
    print ""
    print "Refining profile parmameters"

    # If we have no profile parameters to from sigma d
    if self._profile_parameters is None:
      self._profile_parameters = matrix.col((
        self.sigma_d,
        0,
        self.sigma_d,
        0,
        0,
        self.sigma_d))

    # Create the parameterisation
    parameterisation = SimpleMosaicityParameterisation(self._profile_parameters)

    # Create the refiner and refine
    refiner = ProfileRefiner(
      parameterisation,
      self._refiner_data)
    refiner.refine()

    # Set the profile parameters
    self._profile_parameters = refiner.parameters

  def _refine_crystal(self):
    '''
    Do the crystal parameter refinement

    '''
    print ""
    print "Refining crystal unit cell and orientation parameters"

    # Create the parameterisation
    model = SimpleMosaicityParameterisation(self._profile_parameters)

    # Create the refiner and refine
    refiner = CrystalRefiner(
      self.experiments[0],
      self.reflections,
      model)

    # Set the experiments
    self.experiments[0] = refiner.experiment

  def _plot_distance_from_ewald_sphere(self):
    '''
    Plot distance from Ewald sphere

    '''
    from matplotlib import pylab
    s0 = matrix.col(self.experiments[0].beam.get_s0())
    s2 = self.reflections['s2']
    D = flex.double(s0.length() - matrix.col(s).length() for s in s2)
    Dmean = flex.sum(D) / len(D)
    Dvar = flex.sum(flex.double([(d - Dmean)**2 for d in D])) / len(D)
    print "Mean D: ", Dmean
    print "Variance: ", Dvar
    pylab.hist(D, bins=20)
    pylab.show()


class FinalIntegrator(object):
  '''
  Do the final refinement

  '''

  def __init__(self, experiments, reflections, profile_model):
    '''
    Initialise the refiner

    '''

    # Save some stuff
    self.experiments = experiments
    self.reflections = reflections
    self.profile_model = profile_model

    # Do the processing
    self._compute_bbox()
    self._allocate_shoebox()
    self._extract_shoebox()
    self._compute_mask()
    self._compute_background()
    self._compute_intensity()
    self._compute_partiality()

  def _compute_bbox(self):
    '''
    Compute the bounding box

    '''

    # Get the sigma
    sigma = self.profile_model.sigma()

    # Compute the bounding boxes
    calculator = BBoxCalculatorNew(self.experiments[0], sigma)
    calculator.compute(self.reflections)

    # Select reflections within detector
    x0, x1, y0, y1, _, _ = self.reflections["bbox"].parts()
    xsize, ysize = self.experiments[0].detector[0].get_image_size()
    selection = (x1 > 0) & (y1 > 0) & (x0 < xsize) & (y0 < ysize)
    self.reflections = self.reflections.select(selection)
    print "Filtered reflecions with bbox outside image range"
    print "Kept %d reflections" % len(self.reflections)

  def _allocate_shoebox(self):
    '''
    Allocate the shoebox

    '''
    self.reflections['shoebox'] = flex.shoebox(
      self.reflections['panel'],
      self.reflections['bbox'],
      allocate=True)

  def _compute_mask(self):
    '''
    Compute the reflection mask

    '''
    sigma = self.profile_model.sigma()
    calculator = MaskCalculatorNew(self.experiments[0], sigma)
    calculator.compute(self.reflections)

  def _extract_shoebox(self):
    '''
    Extract the shoebox

    '''
    print "Extracting shoebox from image for %d reflections" % len(self.reflections)
    self.reflections.extract_shoeboxes(self.experiments[0].imageset)

  def _compute_background(self):
    '''
    Compute the reflection background

    '''
    print "Computing background for %d reflections" % len(self.reflections)
    self.reflections.compute_background(self.experiments)

  def _compute_intensity(self):
    '''
    Compute the reflection intensity

    '''
    print "Computing intensity for %d reflections" % len(self.reflections)
    self.reflections.compute_summed_intensity()
    print "%d reflections integrated" % self.reflections.get_flags(
      self.reflections.flags.integrated_sum).count(True)

  def _compute_partiality(self):
    '''
    Compute the partiality

    '''
    s0 = matrix.col(self.experiments[0].beam.get_s0())
    partiality = flex.double(len(self.reflections))
    for k in range(len(self.reflections)):
      s1 = matrix.col(self.reflections[k]['s1'])
      s2 = matrix.col(self.reflections[k]['s2'])
      sbox = self.reflections[k]['shoebox']

      sigma = self.profile_model.sigma()
      R = compute_change_of_basis_operation(s0, s2)
      S = R*sigma*R.transpose()
      mu = R*s2
      assert(abs(1-mu.normalize().dot(matrix.col((0,0,1)))) < 1e-7)

      S11 = matrix.sqr((
        S[0], S[1],
        S[3], S[4]))
      S12 = matrix.col((S[2], S[5]))
      S21 = matrix.col((S[6], S[7])).transpose()
      S22 = S[8]

      mu1 = matrix.col((mu[0], mu[1]))
      mu2 = mu[2]
      partiality[k] = exp(-0.5*(s0.length()-mu2) * (1/S22) * (s0.length()-mu2))
    self.reflections['partiality'] = partiality



class Integrator(object):
  '''
  Class to perform integration of stills in the following way:

  1. Do an initial integration of strong reflections
  2. Refine profile and crystal parameters
  3. Do a final integration of all reflections

  '''

  def __init__(self,
               experiments,
               reflections,
               params=None):
    '''
    Initialise the integrator

    '''

    # Only use single experiment at the moment
    if len(experiments) > 1:
      raise RuntimeError('Only 1 experiment can be processed')

    # Set the parameters
    if params is not None:
      self.params = params
    else:
      self.params = phil_scope.extract(parse(""))

    # Save some stuff
    self.experiments = experiments
    self.reflections = reflections
    self.profile_model = None
    self.sigma_d = None

  def initial_integration(self):
    '''
    Do an initial integration of the strong spots

    '''
    integrator = InitialIntegrator(
      self.params,
      self.experiments,
      self.reflections)
    self.reflections = integrator.reflections
    self.sigma_d = integrator.sigma_d

  def refine(self):
    '''
    Do the refinement of profile and crystal parameters

    '''
    refiner = Refiner(
      self.experiments,
      self.reflections,
      self.sigma_d,
      self.params)
    self.experiments = refiner.experiments
    self.reflections = refiner.reflections
    self.profile_model = refiner.profile_model

  def predict(self):
    '''
    Predict the reflections

    '''
    print ""
    print "Predicting reflections"

    # Set a resolution range
    if self.params.prediction.d_min is None:
      s0 = self.experiments[0].beam.get_s0()
      d_min = self.experiments[0].detector.get_max_resolution(s0)
    else:
      d_min = self.params.predictions.d_min

    # Create the index generator
    index_generator = IndexGenerator(
      self.experiments[0].crystal.get_unit_cell(),
      self.experiments[0].crystal.get_space_group().type(),
      d_min)

    # Get an array of miller indices
    miller_indices_to_test = index_generator.to_array()
    print "Generated %d miller indices" % len(miller_indices_to_test)

    # Get the covariance matrix
    sigma = self.profile_model.sigma()

    # Create the predictor
    predictor = Predictor(
      self.experiments[0],
      sigma,
      self.params.prediction.probability)

    # Do the prediction
    self.reflections = predictor.predict(miller_indices_to_test)
    print "Predicted %d reflections" % len(self.reflections)

  def integrate(self):
    '''
    Do an final integration of the predicted

    '''
    integrator = FinalIntegrator(
      self.experiments,
      self.reflections,
      self.profile_model)
    self.reflections = integrator.reflections

    # Delete shoeboxes if necessary
    if not self.params.debug.output.shoeboxes:
      del self.reflections['shoebox']
