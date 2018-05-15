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
from dials.algorithms.refinement.refinement_helpers import corrgram
from dials_scratch.jmp.potato.refiner import RefinerData
from dials_scratch.jmp.potato.refiner import Refiner as ProfileRefiner
from dials_scratch.jmp.potato.refiner import print_eigen_values_and_vectors
from dials_scratch.jmp.potato.parameterisation import ModelState
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
import matplotlib

# Set matplotlib backend
matplotlib.use("agg", warn=False)

logger = logging.getLogger("dials." + __name__)

# Parameters
phil_scope = parse('''

  profile
  {

    rlp_mosaicity {

      anisotropic = True
        .type = bool

    }

    wavelength_spread {

      model = *delta gaussian
        .type = choice

    }

    unit_cell {

      fixed = False
        .type = bool

    }

    orientation {

      fixed = False
        .type = bool

    }

  }

  refinement {

    max_centroid_distance = 2
      .type = float

    n_macro_cycles = 1
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


class Indexer(object):
  '''
  A class to reindex the strong spot list

  '''

  def __init__(self, params, experiments, reflections):
    '''
    Do the indexing

    '''

    # Save some state
    self.params = params
    self.experiments = experiments
    self.reflections = reflections

    # Do the processing
    self._index()
    self._predict()
    self._filter_reflections_based_on_centroid_distance()

  def _index(self):
    '''
    Index the strong spots

    '''

    # Get some stuff from experiment
    A = matrix.sqr(self.experiments[0].crystal.get_A())
    s0 = matrix.col(self.experiments[0].beam.get_s0())
    detector = self.experiments[0].detector

    # Index all the reflections
    xyz_list = self.reflections['xyzobs.px.value']
    miller_index = self.reflections['miller_index']
    selection = flex.size_t()
    for i in range(len(self.reflections)):

      # Get the observed pixel coordinate
      x, y, _ = xyz_list[i]

      # Get the lab coord
      s1 = matrix.col(detector[0].get_pixel_lab_coord((x,y))).normalize()*s0.length()

      # Get the reciprocal lattice vector
      r = s1 - s0

      # Compute the fractional miller index
      hf = A.inverse() * r

      # Compute the integer miller index
      h = matrix.col((
        int(floor(hf[0] + 0.5)),
        int(floor(hf[1] + 0.5)),
        int(floor(hf[2] + 0.5))))

      # If its not indexed as 0, 0, 0 then append
      if h != matrix.col((0, 0, 0)):
        miller_index[i] = h
        selection.append(i)

    # Print some info
    logger.info("Reindexed %d/%d input reflections" % (
      len(selection),
      len(self.reflections)))

    # Select all the indexed reflections
    self.reflections = self.reflections.select(selection)

  def _predict(self):
    '''
    Predict the position of the spots

    '''

    # Get some stuff from experiment
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
    logger.info("Do prediction for %d reflections" % len(self.reflections))

  def _filter_reflections_based_on_centroid_distance(self):
    '''
    Filter reflections too far from predicted position

    '''
    Xobs, Yobs, _ = self.reflections['xyzobs.px.value'].parts()
    Xcal, Ycal, _ = self.reflections['xyzcal.px'].parts()
    D = flex.sqrt((Xobs-Xcal)**2 + (Yobs-Ycal)**2)
    selection = D < self.params.refinement.max_centroid_distance
    self.reflections = self.reflections.select(selection)
    logger.info("Selected %d reflections with centroid-prediction distance < %d pixels" % (
      len(self.reflections),
      self.params.refinement.max_centroid_distance))



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

  def _compute_sigma_d(self):
    '''
    Compute and initial spot size estimate

    '''
    logger.info("Computing initial sigma d estimate for %d reflections" % len(self.reflections))
    compute_sigma_d = ComputeEsdBeamDivergence(
      self.experiments[0].detector,
      self.reflections)
    self.sigma_d = compute_sigma_d.sigma()
    logger.info("Sigma D: %.5f degrees" % (self.sigma_d * 180 / pi))
    logger.info("")

  def _compute_bbox(self):
    '''
    Compute the bounding box

    '''

    logger.info("Computing the bounding box for %d reflections" % len(self.reflections))

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
    logger.info("Creating the foreground mask for %d reflections" % len(self.reflections))

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
    logger.info("Extracting shoebox from image for %d reflections" % len(self.reflections))
    self.reflections.extract_shoeboxes(self.experiments[0].imageset)

  def _compute_background(self):
    '''
    Compute the reflection background

    '''
    logger.info("Computing background for %d reflections" % len(self.reflections))
    self.reflections.compute_background(self.experiments)

  def _compute_intensity(self):
    '''
    Compute the reflection intensity

    '''
    logger.info("Computing intensity for %d reflections" % len(self.reflections))
    self.reflections.compute_summed_intensity()
    logger.info("%d reflections integrated" % self.reflections.get_flags(
      self.reflections.flags.integrated_sum).count(True))

  def _print_shoeboxes(self):
    '''
    Print the shoeboxes

    '''
    sbox = self.reflections['shoebox']
    for r in range(len(sbox)):
      data = sbox[r].data
      mask = sbox[r].mask
      logger.info(mask.as_numpy_array())
      logger.info(data.as_numpy_array())


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

    # Set the M params
    self.M_params = flex.double((sigma_d, 0, sigma_d, 0, 0, sigma_d))

    # Set the L params
    if self.params.profile.wavelength_spread.model == "gaussian":
      self.L_params = flex.double((sigma_d,))
    elif self.params.profile.wavelength_spread.model == "delta":
      self.L_params = flex.double([0])
    else:
      raise RuntimeError("Unknown wavelength spread model %s" %
                         self.params.profile.wavelength_spread.model)

    # Preprocess the reflections
    self._preprocess()

    # Do the macro cycles of refinement between refining the profile parameters
    # and refining the crystal orientation and unit cell
    for cycle in range(self.params.refinement.n_macro_cycles):
      logger.info("")
      logger.info("Macro cycle %d" % (cycle+1))
      self._refine_profile()

  def _preprocess(self):
    '''
    Preprocess the reflections

    '''

    # Make some plots
    if self.params.debug.output.plots:
      self._plot_distance_from_ewald_sphere()

    # Construct the profile refiner data
    self._refiner_data = RefinerData.from_reflections(
      self.experiments[0],
      self.reflections)

  def _refine_profile(self):
    '''
    Do the profile refinement

    '''
    logger.info("")
    logger.info("Refining profile parmameters")

    # Create the parameterisation
    state = ModelState(
      self.experiments[0],
      fix_orientation       = self.params.profile.orientation.fixed,
      fix_unit_cell         = self.params.profile.unit_cell.fixed,
      fix_wavelength_spread = self.params.profile.wavelength_spread.model == "delta")

    # Set the parameters
    state.set_M_params(self.M_params)
    state.set_L_params(self.L_params)

    # Create the refiner and refine
    refiner = ProfileRefiner(state, self._refiner_data)
    refiner.refine()

    # Set the profile parameters
    self.M_params = state.get_M_params()
    self.L_params = state.get_L_params()

    # Set the mosaicity
    self.experiments[0].crystal.mosaicity = state.get_M()

    # Plot the corrgram
    if self.params.debug.output.plots:
      self._plot_corrgram(
        refiner.correlation(),
        refiner.labels())

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
    fig, ax1 = pylab.subplots(figsize=(10,8))
    ax1.hist(D, bins=max(5, min(int(0.2*len(s2)), 20)))
    ax1.set_xlabel("Distance from Ewald sphere (epsilon)")
    ax1.axvline(x=0, color='black')
    ax1.set_title("Mean(epsilon) = %.2e, Variance(epsilon) = %.2e" % (
      Dmean,
      Dvar))
    fig.savefig("epsilon_distribution.png", dpi=300)
    fig.clf()

  def _plot_corrgram(self, corrmat, labels):
    '''
    Plot a corrgram of correlations between parameters

    '''
    plt = corrgram(corrmat, labels)
    plt.savefig("corrgram.png", dpi=300)
    plt.clf()


class FinalIntegrator(object):
  '''
  Do the final refinement

  '''

  def __init__(self, params, experiments, reflections):
    '''
    Initialise the refiner

    '''

    # Save some stuff
    self.experiments = experiments
    self.reflections = reflections

    # Do the processing
    self._compute_bbox()
    self._allocate_shoebox()
    self._extract_shoebox()
    self._compute_mask()
    self._compute_background()
    self._compute_intensity()
    self._compute_partiality()

    # Plot the partialities
    if params.debug.output.plots:
      self._plot_partiality()

  def _compute_bbox(self):
    '''
    Compute the bounding box

    '''

    # Get the sigma
    sigma = self.experiments[0].crystal.mosaicity

    # Compute the bounding boxes
    calculator = BBoxCalculatorNew(self.experiments[0], sigma)
    calculator.compute(self.reflections)

    # Select reflections within detector
    x0, x1, y0, y1, _, _ = self.reflections["bbox"].parts()
    xsize, ysize = self.experiments[0].detector[0].get_image_size()
    selection = (x1 > 0) & (y1 > 0) & (x0 < xsize) & (y0 < ysize)
    self.reflections = self.reflections.select(selection)
    logger.info("Filtered reflections with bbox outside image range")
    logger.info("Kept %d reflections" % len(self.reflections))

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
    sigma = self.experiments[0].crystal.mosaicity
    calculator = MaskCalculatorNew(self.experiments[0], sigma)
    calculator.compute(self.reflections)

  def _extract_shoebox(self):
    '''
    Extract the shoebox

    '''
    logger.info("Extracting shoebox from image for %d reflections" % len(self.reflections))
    self.reflections.extract_shoeboxes(self.experiments[0].imageset)

  def _compute_background(self):
    '''
    Compute the reflection background

    '''
    logger.info("Computing background for %d reflections" % len(self.reflections))
    self.reflections.compute_background(self.experiments)

  def _compute_intensity(self):
    '''
    Compute the reflection intensity

    '''
    logger.info("Computing intensity for %d reflections" % len(self.reflections))
    self.reflections.compute_summed_intensity()
    logger.info("%d reflections integrated" % self.reflections.get_flags(
      self.reflections.flags.integrated_sum).count(True))

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

      sigma = self.experiments[0].crystal.mosaicity
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

  def _plot_partiality(self):
    '''
    Plot the partiality

    '''
    from matplotlib import pylab
    P = self.reflections['partiality']
    fig, ax1 = pylab.subplots(figsize=(10,8))
    ax1.hist(P, bins=max(5, int(min(0.2*len(P)), 20)))
    ax1.set_xlabel("Scale factor")
    fig.savefig("partiality.png", dpi=300)
    fig.clf()


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
    self.sigma_d = None

  def reindex_strong_spots(self):
    '''
    Reindex the strong spot list

    '''
    indexer = Indexer(
      self.params,
      self.experiments,
      self.reflections)
    self.reflections = indexer.reflections

  def integrate_strong_spots(self):
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

  def predict(self):
    '''
    Predict the reflections

    '''
    logger.info("")
    logger.info("Predicting reflections")

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
    logger.info("Generated %d miller indices" % len(miller_indices_to_test))

    # Get the covariance matrix
    sigma = self.experiments[0].crystal.mosaicity

    # Create the predictor
    predictor = Predictor(
      self.experiments[0],
      sigma,
      self.params.prediction.probability)

    # Do the prediction
    self.reflections = predictor.predict(miller_indices_to_test)
    logger.info("Predicted %d reflections" % len(self.reflections))

  def integrate(self):
    '''
    Do an final integration of the predicted

    '''
    integrator = FinalIntegrator(
      self.params,
      self.experiments,
      self.reflections)
    self.reflections = integrator.reflections

    # Delete shoeboxes if necessary
    if not self.params.debug.output.shoeboxes:
      del self.reflections['shoebox']
