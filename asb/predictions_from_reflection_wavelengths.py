from __future__ import division
from dials.array_family import flex
import math

"""
Experimental functions for alternate means of predicting reflection central impacts using a mosaic model and a bandpass.
"""

def estimate_bandpass(reflections):
  """
  Estimate a minimum and maximum wavelength from the column reflection_wavelength_from_pixels
  At present, this is just mean +/- sigma
  """
  stats = flex.mean_and_variance(reflections['reflection_wavelength_from_pixels'])
  sigma = stats.unweighted_sample_standard_deviation()
  wavelength_min = stats.mean() - sigma
  wavelength_max = stats.mean() + sigma
  return wavelength_min, wavelength_max

def refine_wavelengths(experiments, reflections, initial_mosaic_parameters, refine_bandpass=False):
  """ Simplex minimizer to refine per pixel wavelengths by refining the mosaic model and bandpass for
      each experiment """

  from scitbx.simplex import simplex_opt
  class simplex_minimizer(object):
    """Class for refining mosaic parameters """
    def __init__(self, experiments, reflections, initial_mosaic_parameters, refine_bandpass=False):
      """ Initialize the minimizer and perform the minimization
      @param experiments ExperimentList
      @param reflections flex.reflection_table
      @param initial_mosaic_parameters Tuple of domain size (angstroms) and half mosaic angle (degrees)
      @param refine_bandpass If True, refine band pass for each experiment.
      """
      self.experiments = experiments
      self.reflections = reflections
      self.refine_bandpass = refine_bandpass
      self.x = flex.double(initial_mosaic_parameters)
      self.n = 2
      if refine_bandpass:
        for expt_id, expt in enumerate(experiments):
          refls = reflections.select(reflections['id'] == expt_id)
          wavelength_min, wavelength_max = estimate_bandpass(refls)
          expt.crystal.bandpass = wavelength_min, wavelength_max
          # actually refine the midpoint and the width.  That ensures the two values don't cross over.
          self.x.append(wavelength_max - wavelength_min)
          self.x.append((wavelength_min+wavelength_max)/2)
          self.n += 2

      self.starting_simplex = []
      for i in xrange(self.n+1):
        self.starting_simplex.append((0.9+(flex.random_double(self.n)/10))*self.x)

      self.optimizer = simplex_opt( dimension = self.n,
                                    matrix    = self.starting_simplex,
                                    evaluator = self,
                                    tolerance = 1e-1)
      self.x = self.optimizer.get_solution()

    def target(self, vector):
      """ Compute the functional """
      print "Starting target", list(vector[0:2])
      if (vector < 0).count(True) > 0:
        return 1e6
      if self.refine_bandpass:
        for i in xrange(len(self.experiments)):
          width, midpoint = vector[(i*2)+2], vector[(i*2)+3]
          self.experiments[i].crystal.bandpass = midpoint - width, midpoint + width
          if i == 0:
            print "Experiment 0 bandpass:", self.experiments[i].crystal.bandpass

      # recompute the wavelengths
      self.reflections = tophat_vector_wavelengths(self.experiments, self.reflections, vector[0:2])
      stats = flex.mean_and_variance(12398.4/self.reflections['reflection_wavelength_from_mosaicity_and_bandpass'])
      print "Mean energy: %.1f +/- %.1f"%(stats.mean(), stats.unweighted_sample_standard_deviation()), vector[0], vector[1]
      self.reflections = predictions_from_per_reflection_energies(experiments, self.reflections, 'reflection_wavelength_from_mosaicity_and_bandpass', 'mosbandp')

      # Miminimze the deltaXY RMSD of observations-predictions
      rmsd = math.sqrt(flex.sum((self.reflections['xyzobs.mm.value']-self.reflections['xyzcal.mm.mosbandp']).norms()**2)/len(self.reflections))
      print "Dataset RMSD (microns)", rmsd * 1000
      return rmsd
      # Alternative formulation of target. Minimize summed squared differences between computed wavelengths for each reflection.
      #f = flex.sum((self.reflections['reflection_wavelength_from_mosaicity_and_bandpass']-self.reflections['reflection_wavelength_from_pixels'])**2)
      #print "Functional", f
      #return f

  minimizer = simplex_minimizer(experiments, reflections, initial_mosaic_parameters, refine_bandpass=refine_bandpass)
  final_mosaic_parameters = minimizer.x[0:2]
  print "Final mosaic parameters", list(final_mosaic_parameters), experiments[0].crystal.bandpass
  if refine_bandpass:
    for i in xrange(len(experiments)):
      width, midpoint = minimizer.x[(i*2)+2], minimizer.x[(i*2)+3]
      experiments[i].crystal.bandpass = midpoint - width, midpoint + width

  reflections = predictions_from_per_reflection_energies(experiments, minimizer.reflections, 'reflection_wavelength_from_mosaicity_and_bandpass', 'mosbandp')
  return reflections

def tophat_vector_wavelengths(experiments, reflections, mosaic_parameters):
  """
  Given a set of mosaic parameters, use vectors to estimate a wavelength for each reflection

  Details:
  For a given reflection, the reciprocal lattice point vector q = Ah, where A is the reciprocal
  A matrix and h is the reflection's miller index.  Construct a vector e1 orthagonal to s0 and q.
  Construct 4 more vectors of length equal to the magnitude of q, and that lie in the plane that
  is normal to e1:
  q_mos_inner and q_mos_outer: q vectors rotated into or out of the Ewald sphere by an angle equal
  to the half mosaic angle + the angle inscribed by adding an arc length equal to 2/domain size
  (angstroms). Call this angle the combined mosaic angle approximation.
  q_wavelength_min and q_wavelength_max: q vectors rotated on an ewald sphere with radius
  1/bandpass minimum or 1/bandpass maximum.

  Consider now the pairs of vectors wavelength min/max and q_mos inner/outer.  If neither q_mos
  vectors lies between the two wavelength vectors, then assign a refletion's wavelength to
  either wavelength min or max, depending on which is closest.  Otherwise, find two vectors
  that lie between wavelength min/max. For example if q_mos inner lies between them, but outer
  does not, the two vectors will be q_mos inner and q wavelength min.  If both q_mos inner and
  outer lie between wavelength min/max, then the two vetors will be q_mos inner and outer.

  Once the two vectors have been identified, define a new q vector which is the average of the
  two vectors.  Determine the wavelength that would allow this q vector to be in the diffracting
  condition.  Report that wavelength in the column reflection_wavelength_from_mosaicity_and_bandpass.

  Because these determinations involve hard cutoffs instead of gaussians, the wavelengths are
  determined in a manner similar to the overlap of tophat functions.

  @param experiments ExperimentList. If crystal.band_pass is set, use it. Otherwise, estimate the
  band pass using estimate_bandpass
  @param reflections flex.reflection_table Needs to contain the column
  reflection_wavelength_from_pixels
  @param mosaic_parameters Tuple of domain size (angstroms) and half mosaic angle (degrees)

  """

  if 'reflection_wavelength_from_pixels' not in reflections:
    return reflections

  domain_size_ang, half_mosaicity_deg = mosaic_parameters
  print "Computing per-reflection wavelengths from domain size", domain_size_ang, "(ang), half mosaic angle", half_mosaicity_deg, "(deg), and bandpass derived from each image"

  table = flex.reflection_table()

  new_wavelengths = flex.double()

  # Keep track of the various cases
  case_0 = case_1 = case_2 = case_3 = case_4 = case_5 = 0

  for expt_id, expt in enumerate(experiments):
    refls = reflections.select(reflections['id'] == expt_id)
    table.extend(refls)

    # Determine how to obtain the bandpass
    if hasattr(expt.crystal, 'bandpass'):
      wavelength_min, wavelength_max = expt.crystal.bandpass
    else:
      wavelength_min, wavelength_max = estimate_bandpass(refls)
      expt.crystal.bandpass = wavelength_min, wavelength_max

    unit_s0 = flex.vec3_double(len(refls), expt.beam.get_s0()).each_normalize()
    wavelength = expt.beam.get_wavelength()

    q = flex.mat3_double(len(refls), expt.crystal.get_A()) * refls['miller_index'].as_vec3_double()
    e1 = q.cross(unit_s0).each_normalize()

    # length of an arc l = 2pir * angle/2pi = r*angle. So angle = l/r
    combined_mosaic_angle_approximation = ((2/domain_size_ang)/q.norms()) + (half_mosaicity_deg*math.pi/180)
    q_mos_inner = q.rotate_around_origin(e1, -combined_mosaic_angle_approximation)
    q_mos_outer = q.rotate_around_origin(e1,  combined_mosaic_angle_approximation)
    # Z: angle between q and s0
    z_mos_inner = (math.pi-q_mos_inner.angle(unit_s0))
    z_mos_outer = (math.pi-q_mos_outer.angle(unit_s0))
    lmbda_bigger  = flex.cos(z_mos_inner)/(q.norms()/2)
    lmbda_smaller = flex.cos(z_mos_outer)/(q.norms()/2)

    z_wavelength_min = flex.acos(wavelength_min*q.norms()/2)
    z_wavelength_max = flex.acos(wavelength_max*q.norms()/2)
    q_wavelength_min = unit_s0.rotate_around_origin(e1, math.pi+z_wavelength_min)*q.norms()
    q_wavelength_max = unit_s0.rotate_around_origin(e1, math.pi+z_wavelength_max)*q.norms()

    assert (lmbda_smaller<lmbda_bigger).count(False) == 0
    sel = flex.bool(len(refls), True)
    image_wavelengths = flex.double(len(refls), 0)

    q_inner = flex.vec3_double()
    q_outer = flex.vec3_double()

    # Iterate through the reflections and sort each into one of the 6 cases
    for i in xrange(len(refls)):
      # Cases 0 and 1: both q_mos inner and outer are outside of the wavelength min/max vectors
      if lmbda_smaller[i] > wavelength_max:
        image_wavelengths[i] = wavelength_max
        sel[i] = False
        case_0 += 1
        continue
      elif lmbda_bigger[i] < wavelength_min:
        image_wavelengths[i] = wavelength_min
        sel[i] = False
        case_1 += 1
        continue

      sel[i] = True

      # Case 2: q_mos outer is between wavelengths min and max so use q_mos outer
      # Case 3: q_mos outer is outside of wavelength min so use wavelength min
      if lmbda_smaller[i] >= wavelength_min:
        q_outer.append(q_mos_outer[i])
        case_2 += 1
      else:
        q_outer.append(q_wavelength_min[i])
        case_3 += 1

      # Case 4: q_mos inner is between wavelengths min and max so use q_mos inner
      # Case 5: q_mos inner is outside of wavelength max so use wavelength max
      if lmbda_bigger[i] <= wavelength_max:
        q_inner.append(q_mos_inner[i])
        case_4 += 1
      else:
        q_inner.append(q_wavelength_max[i])
        case_5 += 1

    # Compute new reflection wavelengths
    new_q = (q_inner + q_outer)*0.5
    z = math.pi-new_q.angle(unit_s0.select(sel))
    lmbda = flex.cos(z)/(new_q.norms()/2)
    image_wavelengths.set_selected(sel, lmbda)
    new_wavelengths.extend(image_wavelengths)

  assert (new_wavelengths <= 0).count(True) == 0
  reflections['reflection_wavelength_from_mosaicity_and_bandpass'] = new_wavelengths
  print "CASES", case_0, case_1 ,case_2 ,case_3, case_4, case_5

  return reflections

def predictions_from_per_reflection_energies(experiments, reflections, tag, dest):
  """
  Give a new set of energies for each reflection, recompute reflection locations and
  delta psi values
  @param experiments ExperimentList
  @param reflections flex.reflection_table
  @param tag Column name for the per reflection wavelength
  @param dest Predictions will be in put in the columns delpsical.rad.dest,
  xyzcal.mm.dest and xyzcal.px.dest
  """
  if tag not in reflections:
    return reflections

  table = flex.reflection_table()
  all_new_delpsi = flex.double()
  all_new_pred_mm = flex.vec3_double()
  all_new_pred_px = flex.vec3_double()

  for expt_id, expt in enumerate(experiments):
    refls = reflections.select(reflections['id'] == expt_id)
    table.extend(refls)

    unit_s0 = flex.vec3_double(len(refls), expt.beam.get_s0()).each_normalize()
    s0 = unit_s0 * (1/refls[tag])

    # Calculate the reciprocal space vector and required unit vectors
    q = flex.mat3_double(len(refls), expt.crystal.get_A()) * refls['miller_index'].as_vec3_double()
    assert (q.norms() > 0).count(False) == 0
    e1 = q.cross(unit_s0).each_normalize()
    c0 = unit_s0.cross(e1).each_normalize()

    # Calculate the vector rotated to the Ewald sphere
    qq = q.norms()**2
    wavelength = refls[tag]
    a = 0.5 * qq * wavelength;
    tmp = qq - a*a;
    assert (tmp > 0.0).count(False) == 0
    b = flex.sqrt(tmp)
    r = -1.0 * a * unit_s0 + b * c0;

    # Calculate delpsi value
    q0 = q.each_normalize()
    q1 = q0.cross(e1).each_normalize()
    delpsi = -1.0 * flex.atan2(r.dot(q1), r.dot(q0))
    all_new_delpsi.extend(delpsi)

    # Calculate the Ray
    s1 = (s0 + r).each_normalize() * s0.norms()
    pred_mm = flex.vec3_double()
    pred_px = flex.vec3_double()
    for i in xrange(len(refls)):
      panel = expt.detector[refls['panel'][i]]
      pred = panel.get_ray_intersection(s1[i])
      pred_mm.append((pred[0], pred[1], 0))
      pred = panel.millimeter_to_pixel(pred)
      pred_px.append((pred[0], pred[1], 0))
    all_new_pred_mm.extend(pred_mm)
    all_new_pred_px.extend(pred_px)

  table['delpsical.rad.%s'%dest] = all_new_delpsi
  table['xyzcal.mm.%s'%dest] = all_new_pred_mm
  table['xyzcal.px.%s'%dest] = all_new_pred_px
  return table
