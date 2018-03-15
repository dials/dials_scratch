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
from dials.algorithms.spot_prediction import IndexGenerator
from scitbx.linalg import eigensystem
from dials.array_family import flex
from scitbx import matrix
from math import pi, sqrt, floor, ceil, exp
from dials.algorithms.shoebox import MaskCode
import logging

logger = logging.getLogger(__name__)


class Predictor(object):

  def __init__(self, experiment, parameters, dmin=None):

    print ""
    print "Predicting reflections"

    # Set a resolution range
    if dmin is None:
      s0 = experiment.beam.get_s0()
      dmin = experiment.detector.get_max_resolution(s0)

    # Create the index generator
    index_generator = IndexGenerator(
      experiment.crystal.get_unit_cell(),
      experiment.crystal.get_space_group().type(),
      dmin)

    # Get an array of miller indices
    miller_indices_to_test = index_generator.to_array()
    print "Generated %d miller indices" % len(miller_indices_to_test)

    # Get the covariance matrix
    sigma = SimpleMosaicityParameterisation(parameters).sigma()
    sigma_inv = sigma.inverse()

    # Compute quantile
    quantile = chisq_quantile(3, 0.997)

    # Get stuff from experiment
    A = matrix.sqr(experiment.crystal.get_A())
    s0 = matrix.col(experiment.beam.get_s0())

    # Loop through miller indices and check each is in range
    print "Checking reflections against MVN quantile %f with Mahalabonis distance %f" % (
      0.997, sqrt(quantile))
    panel = experiment.detector[0]
    miller_indices = flex.miller_index()
    entering = flex.bool()
    s1_list = flex.vec3_double()
    s2_list = flex.vec3_double()
    xyzcalpx = flex.vec3_double()
    xyzcalmm = flex.vec3_double()
    panel_list = flex.size_t()
    for h in miller_indices_to_test:

      r = A * h
      s2 = s0 + r
      s3 = s2.normalize()*s0.length()

      d = ((s3-s2).transpose()*sigma_inv*(s3-s2))[0]

      if d < quantile:
        e = s2.length() < s0.length()

        R = compute_change_of_basis_operation(s0, s2)

        S =R*sigma*R.transpose()
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

        mubar = mu1 + S12*(1/S22)*(s0.length()-mu2)
        v = matrix.col((mubar[0], mubar[1], s0.length())).normalize()*s0.length()
        s1 = R.transpose()*v

        try:
          xymm = panel.get_ray_intersection(s1)
        except Exception:
          continue

        xypx = panel.millimeter_to_pixel(xymm)
        miller_indices.append(h)
        entering.append(e)
        panel_list.append(0)
        s1_list.append(s1)
        s2_list.append(s2)
        xyzcalpx.append((xypx[0], xypx[1], 0))
        xyzcalmm.append((xymm[0], xymm[1], 0))

    self._reflections = flex.reflection_table()
    self._reflections['miller_index'] = miller_indices
    self._reflections['entering'] = entering
    self._reflections['s1'] = s1_list
    self._reflections['s2'] = s2_list
    self._reflections['xyzcal.px'] = xyzcalpx
    self._reflections['xyzcal.mm'] = xyzcalmm
    self._reflections['panel'] = panel_list
    self._reflections['id'] = flex.int(len(self._reflections), 0)
    print "Predicted %d reflections" % len(self._reflections)


  def reflections(self):
    return self._reflections


class BBoxCalculatorNew(object):

  def __init__(self, experiment, parameters):
    self.experiment = experiment
    self.parameters = parameters

  def compute(self, reflections):

    print "Computing bbox for %d reflections" % len(reflections)

    # Compute quantile
    quantile = chisq_quantile(2, 0.997)
    D = sqrt(quantile) * 2
    bbox = flex.int6()
    print "ML: %f" % D
    for i in range(len(reflections)):
      s1 = matrix.col(reflections[i]['s1'])
      s2 = matrix.col(reflections[i]['s2'])


      s0 = matrix.col(self.experiment.beam.get_s0())

      # Ensure our values are ok
      assert s1.length() > 0

      sigma = SimpleMosaicityParameterisation(self.parameters).sigma()
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

      mubar = mu1 + S12*(1/S22)*(s0.length()-mu2)
      Sbar = S11 - S12*(1/S22)*S21

      eigen_decomposition = eigensystem.real_symmetric(Sbar.as_flex_double_matrix())
      Q = matrix.sqr(eigen_decomposition.vectors())
      L = matrix.diag(eigen_decomposition.values())
      max_L = max(L)

      delta = sqrt(max_L) * D


      p1 = mubar + matrix.col((-delta, -delta))
      p2 = mubar + matrix.col((-delta, +delta))
      p3 = mubar + matrix.col((+delta, -delta))
      p4 = mubar + matrix.col((+delta, +delta))

      p1 = matrix.col((p1[0], p1[1], s0.length())).normalize()
      p2 = matrix.col((p2[0], p2[1], s0.length())).normalize()
      p3 = matrix.col((p3[0], p3[1], s0.length())).normalize()
      p4 = matrix.col((p4[0], p4[1], s0.length())).normalize()

      x1 = R.transpose()*p1
      x2 = R.transpose()*p2
      x3 = R.transpose()*p3
      x4 = R.transpose()*p4

      xy1 = self.experiment.detector[0].get_ray_intersection_px(x1)
      xy2 = self.experiment.detector[0].get_ray_intersection_px(x2)
      xy3 = self.experiment.detector[0].get_ray_intersection_px(x3)
      xy4 = self.experiment.detector[0].get_ray_intersection_px(x4)

      xx = (xy1[0], xy2[0], xy3[0], xy4[0])
      yy = (xy1[1], xy2[1], xy3[1], xy4[1])
      x0, x1 = int(floor(min(xx)))-1, int(ceil(max(xx)))+1
      y0, y1 = int(floor(min(yy)))-1, int(ceil(max(yy)))+1
      assert x1 > x0
      assert y1 > y0
      bbox.append((x0, x1, y0, y1, 0, 1))

    reflections['bbox'] = bbox

    x0, x1, y0, y1, _, _ = bbox.parts()
    xsize, ysize = self.experiment.detector[0].get_image_size()
    selection = (x1 > 0) & (y1 > 0) & (x0 < xsize) & (y0 < ysize)

    reflections = reflections.select(selection)
    print "Filtered reflecions with bbox outside image range"
    print "Kept %d reflections" % len(reflections)
    return reflections

class MaskCalculatorNew(object):

  def __init__(self, experiment, parameters):
    self.experiment = experiment
    self.parameters = parameters

  def compute(self, reflections):
    print "Computing mask for %d reflections" % len(reflections)

    # Compute quantile
    quantile = chisq_quantile(2, 0.997)
    D = quantile
    panel = self.experiment.detector[0]
    print "ML: %f" % sqrt(D)
    for k in range(len(reflections)):
      s1 = matrix.col(reflections[k]['s1'])
      s2 = matrix.col(reflections[k]['s2'])
      sbox = reflections[k]['shoebox']


      s0 = matrix.col(self.experiment.beam.get_s0())
      s0_length = s0.length()

      # Ensure our values are ok
      assert s1.length() > 0

      sigma = SimpleMosaicityParameterisation(self.parameters).sigma()
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

      mubar = mu1 + S12*(1/S22)*(s0.length()-mu2)
      Sbar = S11 - S12*(1/S22)*S21
      Sbar_inv = Sbar.inverse()

      mask = sbox.mask
      x0, x1, y0, y1, _, _ = sbox.bbox
      cs = CoordinateSystem2d(s0, s2)
      for j in range(mask.all()[1]):
        for i in range(mask.all()[2]):
          ii = i + x0
          jj = j + y0
          s1 = panel.get_pixel_lab_coord((ii,jj))
          s2 = panel.get_pixel_lab_coord((ii+1,jj))
          s3 = panel.get_pixel_lab_coord((ii,jj+1))
          s4 = panel.get_pixel_lab_coord((ii+1,jj+1))
          s1 = matrix.col(s1).normalize() * s0_length
          s2 = matrix.col(s2).normalize() * s0_length
          s3 = matrix.col(s3).normalize() * s0_length
          s4 = matrix.col(s4).normalize() * s0_length
          x1 = matrix.col(cs.from_beam_vector(s1))
          x2 = matrix.col(cs.from_beam_vector(s2))
          x3 = matrix.col(cs.from_beam_vector(s3))
          x4 = matrix.col(cs.from_beam_vector(s4))
          d1 = ((x1 - mubar).transpose()*Sbar_inv*(x1 - mubar))[0]
          d2 = ((x2 - mubar).transpose()*Sbar_inv*(x2 - mubar))[0]
          d3 = ((x3 - mubar).transpose()*Sbar_inv*(x3 - mubar))[0]
          d4 = ((x4 - mubar).transpose()*Sbar_inv*(x4 - mubar))[0]
          if min([d1, d2, d3, d4]) < D:
            mask[0,j,i] = mask[0,j,i] | MaskCode.Foreground
          else:
            mask[0,j,i] = mask[0,j,i] | MaskCode.Background


class Integrator(object):
  '''
  Class to perform integration of stills in the following way:

  1. Refine profile model
  2. Refine crystal orientation
  3. Integrate reflections

  '''

  def __init__(self,
               experiments,
               reflections):

    # Only use single experiment at the moment
    if len(experiments) > 1:
      raise RuntimeError('Only 1 experiment can be processed')

    # Save some stuff
    self.experiments = experiments
    self.reflections = reflections
    self.n_macro_cycles = 3
    self._profile_parameters = None

  def initial_integration(self):
    '''
    Do an initial integration of the strong spots

    '''
    self._compute_sigma_d()
    self._compute_bbox()
    self._allocate_shoebox()
    self._extract_shoebox()
    self._compute_mask()
    self._compute_background()
    self._compute_intensity()

  def refine(self):
    '''
    Do the refinement of profile and crystal parameters

    '''

    # Preprocess the reflections
    self._preprocess()

    # Do the macro cycles of refinement between refining the profile parameters
    # and refining the crystal orientation and unit cell
    for cycle in range(self.n_macro_cycles):
      print ""
      print "Macro cycle %d" % (cycle+1)
      self._refine_profile()
      self._refine_crystal()

  def predict(self):
    predictor = Predictor(self.experiments[0], self._profile_parameters)
    self.reflections = predictor.reflections()

  def integrate(self):

    self._compute_bbox_new()
    self._allocate_shoebox()
    self._extract_shoebox()
    self._compute_mask_new()
    self._compute_background()
    self._compute_intensity()
    self._compute_partiality()

  def _compute_bbox_new(self):
    calculator = BBoxCalculatorNew(self.experiments[0], self._profile_parameters)
    self.reflections = calculator.compute(self.reflections)

  def _compute_mask_new(self):
    calculator = MaskCalculatorNew(self.experiments[0], self._profile_parameters)
    calculator.compute(self.reflections)

  def _compute_partiality(self):
    s0 = matrix.col(self.experiments[0].beam.get_s0())
    partiality = flex.double(len(self.reflections))
    for k in range(len(self.reflections)):
      s1 = matrix.col(self.reflections[k]['s1'])
      s2 = matrix.col(self.reflections[k]['s2'])
      sbox = self.reflections[k]['shoebox']

      sigma = SimpleMosaicityParameterisation(self._profile_parameters).sigma()
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
    compute_bbox = BBoxCalculator(
      self.experiments[0].crystal,
      self.experiments[0].beam,
      self.experiments[0].detector,
      self.experiments[0].goniometer,
      self.experiments[0].scan,
      self.sigma_d * 6,
      0)

    bbox = compute_bbox(
      self.reflections['s1'],
      self.reflections['xyzcal.px'].parts()[2],
      self.reflections['panel'])
    self.reflections['bbox_old'] = self.reflections['bbox']
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
    mask_foreground = MaskCalculator(
      self.experiments[0].crystal,
      self.experiments[0].beam,
      self.experiments[0].detector,
      self.experiments[0].goniometer,
      self.experiments[0].scan,
      self.sigma_d * 3,
      0)
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
    if False:
      for r in self.reflections:
        print r['shoebox'].data.as_numpy_array()

  def _compute_background(self):
    print "Computing background for %d reflections" % len(self.reflections)
    self.reflections.compute_background(self.experiments)
    if False:
      for r in self.reflections:
        d = r['shoebox'].data
        b = r['shoebox'].background
        m = r['shoebox'].mask
        diff = (d-b)
        mask = (diff > 0).as_1d().as_int()
        mask.reshape(diff.accessor())
        diff = diff*mask.as_double().as_float()
        print (flex.floor(diff)).as_numpy_array()

  def _compute_intensity(self):
    print "Computing intensity for %d reflections" % len(self.reflections)
    self.reflections.compute_summed_intensity()
    print "%d reflections integrated" % self.reflections.get_flags(
      self.reflections.flags.integrated_sum).count(True)

  def _preprocess(self):
    '''
    Preprocess the reflections

    '''

    # Don't trust the predictions in the reflection file.
    self._update_observed_reflection_predictions()

    # Construct the profile refiner data
    self._refiner_data = ProfileRefinerData.from_reflections(
      self.experiments[0],
      self.reflections)

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
    s2 = flex.vec3_double(len(h))
    for i in range(len(self.reflections)):
      r = A*matrix.col(h[i])
      s2[i] = s0 + r
    self.reflections['s2'] = s2

  def _refine_profile(self):
    '''
    Do the profile refinement

    '''
    print ""
    print "Refining profile parmameters"
    if self._profile_parameters is None:
      self._profile_parameters = matrix.col((
        self.sigma_d,
        0,
        self.sigma_d,
        0,
        0,
        self.sigma_d))
    parameterisation = SimpleMosaicityParameterisation(self._profile_parameters)
    refiner = ProfileRefiner(
      parameterisation,
      self._refiner_data)
    refiner.refine()
    self._profile_parameters = refiner.parameters

  def _refine_crystal(self):
    '''
    Do the crystal parameter refinement

    '''
    print ""
    print "Refining crystal unit cell and orientation parameters"
    model = SimpleMosaicityParameterisation(self._profile_parameters)
    refiner = CrystalRefiner(
      self.experiments[0],
      self.reflections,
      model)
    self.experiments[0] = refiner.experiment
