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
from dials_scratch.jmp.stills.potato.profile_refiner import ProfileRefinerData
from dials_scratch.jmp.stills.potato.profile_refiner import ProfileRefiner
from dials_scratch.jmp.stills.potato.profile_refiner import print_eigen_values_and_vectors
from dials_scratch.jmp.stills.potato.crystal_refiner import CrystalRefiner
from dials.array_family import flex
from scitbx import matrix
import logging

logger = logging.getLogger(__name__)

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
    self.n_macro_cycles = 5
    self._profile_parameters = None

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
    pass

  def integrate(self):
    pass

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
    refiner = ProfileRefiner(
      self._refiner_data,
      self._profile_parameters)
    refiner.refine()
    self._profile_parameters = refiner.parameters

  def _refine_crystal(self):
    '''
    Do the crystal parameter refinement

    '''
    print ""
    print "Refining crystal unit cell and orientation parameters"
    refiner = CrystalRefiner(
      self.experiments[0],
      self.reflections,
      self._profile_parameters)
    self.experiments[0] = refiner.experiment
