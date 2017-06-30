from __future__ import division
from libtbx.phil import parse


phil_scope = parse('''

  spot_finding {

    min_spot_size = 2
      .type = int(value_min=1)
      .help = "The minimum number of pixels in a spot"

    max_spot_size = 100
      .type = int(value_min=1)
      .help = "The maximum number of pixels in a spot"

  }

  indexing {

    tolerance = 0.3
      .type = float(value_min=0.01)
      .help = "The tolerance when indexing strong spots"

  }

  refinement {

    n_macro_cycles = 3
      .type = int(value_min=1)
      .help = "The number of refinement macro cycles"

    sample = 100
      .type = int(value_min=10)
      .help = "The sample size for refinement"

    profile {

      num_samples = 5
        .type = int(value_min=1)
        .help = "The number of pixel samples to take in the integral"

      tolerance = 1e-4
        .type = float(value_min=0)
        .help = "The tolerance to stop the optimization"

      target = *maximum_likelihood least_squares
        .type = choice
        .help = "The target function to use"

    }

  }

  integration {

    prediction {

      min_foreground_pixels = 4
        .type = int(value_min=1)
        .help = "The minimum number of foreground pixels in a prediction"

      min_background_pixels = 10
        .type = int(value_min=1)
        .help = "The minimum number of background pixels"

    }

  }

''')


def golden_section_search(f, a, b, tol=1e-3, callback=None):
  from math import sqrt
  gr = (sqrt(5) + 1) / 2
  c = b - (b - a) / gr
  d = a + (b - a) / gr
  while abs(c - d) > tol:

    fc = f(c)
    fd = f(d)

    callback(c, d, fc, fd)

    if fc < fd:
      b = d
    else:
      a = c

    c = b - (b - a) / gr
    d = a + (b - a) / gr

  return (b + a) / 2



def generate_start(values, offset):
  assert len(values) == len(offset)
  start = [values]
  for j, o in enumerate(offset):
    next = values.deep_copy()
    next[j] += o
    start.append(next)
  return start

class simple_simplex(object):
  def __init__(self, values, offset, evaluator, max_iter):
    from scitbx import simplex
    self.n = len(values)
    self.x = values
    self.starting_simplex = generate_start(values, offset)
    self.fcount = 0

    optimizer = simplex.simplex_opt(dimension=self.n,
                                    matrix=self.starting_simplex,
                                    evaluator=evaluator,
                                    tolerance=1e-7,
                                    max_iter=max_iter)

    self.x = optimizer.get_solution()
    return

  def get_solution(self):
    return self.x

  def target(self, vector):
    print "SCORE"
    score = scorify(vector)
    return score





class CrystalRefiner(object):
  '''
  A class to perform refinement of crystal parameters

  '''

  def __init__(self,
               experiment,
               reflections,
               mosaicity,
               params):
    '''
    Perform the refinement

    :param experiment: The experiment
    :param reflections: The reflections
    :param mosaicity: The mosaicity
    :param params: The parameters

    '''
    from dials.algorithms.refinement.parameterisation.crystal_parameters \
      import CrystalUnitCellParameterisation
    from dials.algorithms.refinement.parameterisation.crystal_parameters \
      import CrystalOrientationParameterisation
    from dials_scratch.jmp.stills import Model
    from dials.array_family import flex
    from scitbx import simplex
    from math import sqrt

    # Store the input
    self.experiment = experiment
    self.reflections = reflections
    self.mosaicity = mosaicity
    self.params = params

    # Get the data and image mask
    data = self.experiment.imageset.get_raw_data(0)[0]
    mask = self.experiment.imageset.get_mask(0)[0]

    # Initialise the model
    self.model = Model(
      beam             = self.experiment.beam,
      detector         = self.experiment.detector,
      crystal          = self.experiment.crystal,
      reflections      = self.reflections,
      image_data       = data,
      image_mask       = mask,
      mosaicity        = self.mosaicity,
      bandpass         = 0.0,
      foreground_limit = 0.3,
      background_limit = 0.5,
      num_samples      = self.params.refinement.profile.num_samples)

    # Get the crystal model and the parameterisation
    self.crystal = self.experiment.crystal
    self.cucp = CrystalUnitCellParameterisation(self.crystal)
    self.cop = CrystalOrientationParameterisation(self.crystal)

    # Get the current values and generate some offsets
    values = flex.double(
      self.cucp.get_param_vals() +
      self.cop.get_param_vals())
    offset = flex.double(
      [0.01 * v for v in self.cucp.get_param_vals()] +
      [0.1, 0.1, 0.1])

    # The optimization history
    self.history = []

    # Get the initial cell and initial score
    initial_cell = self.crystal.get_unit_cell()
    initial_score = self.target(values)

    # Perform the optimization
    optimizer = simple_simplex(values, offset, self, 2000)
    result = optimizer.get_solution()
    print 'Initial cell:', initial_cell
    print 'Final cell:  ', self.crystal.get_unit_cell()

    # Compute RMSD
    xcal, ycal, _ = self.model.observed().parts()
    xobs, yobs, _ = self.model.predicted().parts()
    rmsd_x = sqrt(flex.sum((xcal-xobs)**2) / len(xcal))
    rmsd_y = sqrt(flex.sum((ycal-yobs)**2) / len(ycal))
    print 'RMSD X, Y (px): %f, %f' % (rmsd_x, rmsd_y)

  def target(self, vector):
    '''
    The target function

    '''
    from dials.array_family import flex
    from math import sqrt

    # Get the cell and orientation parameters
    cell_parms = self.cucp.get_param_vals()
    orientation_parms = self.cop.get_param_vals()
    assert len(vector) == len(cell_parms) + len(orientation_parms)

    # Update the cell and orientation parameters
    tst_cell = vector[:len(cell_parms)]
    tst_orientation = vector[len(cell_parms):len(cell_parms) + len(orientation_parms)]
    self.cucp.set_param_vals(tst_cell)
    self.cop.set_param_vals(tst_orientation)

    # Update the model
    self.model.crystal = self.crystal
    self.model.update(pixel_lookup=False)

    # Get the observed and predicted position
    T = self.model.success()
    O = self.model.observed()
    C = self.model.predicted()

    # Select only those successes
    selection = T == True
    Xobs, Yobs, Zobs = O.select(selection).parts()
    Xcal, Ycal, Zcal = C.select(selection).parts()

    # Compute the rmsd between observed and calculated
    score = flex.sum((Xobs-Xcal)**2 + (Yobs-Ycal)**2 + (Zobs-Zcal)**2)

    # Append to the history
    self.history.append(tst_cell, tst_orientation, score)

    # Print some info
    print 'Cell: %.3f %.3f %.3f %.3f %.3f %.3f; Phi: %.3f %.3f %.3f; RMSD: %.3f' % (
      tuple(self.crystal.get_unit_cell().parameters()) +
      tuple(tst_orientation) +
      tuple((sqrt(score / len(Xobs)),)))
    return score


class ProfileRefiner(object):
  '''
  A class to perform the profile refinement

  '''

  def __init__(self,
               experiment,
               reflections,
               params):
    '''
    Do the profile refinement

    :param experiment: The experiment
    :param reflections: The reflection list
    :param params: The parameters

    '''
    from dials_scratch.jmp.stills import Model

    # Save the experiments and reflections
    self.experiment = experiment
    self.reflections = reflections
    self.params = params

    # The optimization history
    self.history = []

    # Get the data and image mask
    data = self.experiment.imageset.get_raw_data(0)[0]
    mask = self.experiment.imageset.get_mask(0)[0]

    # Initialise the model
    self.model = Model(
      beam             = self.experiment.beam,
      detector         = self.experiment.detector,
      crystal          = self.experiment.crystal,
      reflections      = self.reflections,
      image_data       = data,
      image_mask       = mask,
      mosaicity        = 0.1,
      bandpass         = 0.0,
      foreground_limit = 0.3,
      background_limit = 0.5,
      num_samples      = self.params.refinement.profile.num_samples)

    def callback(a, b, fa, fb):
      self.history.append((a, b, fa, fb))
      print "  A = %f, B = %f, Fa = %f, Fb = %f" % (a, b, fa, fb)

    # Compute min and max by taking as fraction of hkl then converting using A
    # matrix
    min_mosaicity = 1e-10
    max_mosaicity = 0.01

    # Perform a golden section search for the profile parameter
    print "Refining profile parameters"
    print "  using '%s' target function" % self.params.refinement.profile.target
    self.profile = golden_section_search(
      self.target,
      min_mosaicity,
      max_mosaicity,
      self.params.refinement.profile.tolerance,
      callback)
    assert self.profile < max_mosaicity
    print "Mosaicity = %f" % (self.profile)

  def target(self, mosaicity):
    '''
    The profile model target function

    '''

    # Set the mosaicity
    self.model.mosaicity = mosaicity

    # Update the model
    self.model.update()

    # Choose the least squares or maximum likelihood scorea
    target_function = self.params.refinement.profile.target
    if target_function == 'least_squares':
      score = self.model.least_squares_score()
    elif target_function == 'maximum_likelihood':
      score = -self.model.maximum_likelihood_score()
    else:
      raise RuntimeError("Unknown target: %s" % target_function)

    # Return the score
    return score


class Integrator(object):
  '''
  Top level integrator class for a single still image

  Given a still experiment do the following:

  1) Find strong spots on the image
  2) Index those spots
  3) Use those spots to refine the profile model
  4) Then refine the crystal model parameters
  5) Predict the full set of reflections
  6) Integrate the reflections

  '''

  def __init__(self, experiment, params):
    '''
    Initialise the class

    :param experiment: The experiment to process
    :param params: The configuration parameters

    '''
    assert len(experiment.detector) == 1
    assert len(experiment.imageset) == 1
    self.experiment = experiment
    self.params = params
    self.reflections = None

  def process(self):
    '''
    Process the image

    '''

    # Find the strong spots
    self.find_spots()

    # Index the strong spots
    self.index_spots()

    # Select a sample of the strong spots
    self.sample_reflections()

    # Iteratively refine profile and crystal models
    for i in range(self.params.refinement.n_macro_cycles):
      self.refine_profile()
      self.refine_crystal()

    # Do a final refinement of profile parameters
    self.refine_profile()

    # Integrate the reflections
    self.integrate()

  def find_spots(self, min_spot_size=2, max_spot_size=100):
    '''
    Find the strong spots on the image

    '''
    from dials.algorithms.spot_finding.threshold import XDSThresholdStrategy
    from dials.model.data import PixelList
    from dials.model.data import PixelListLabeller
    from dials.array_family import flex

    print ""
    print "-" * 80
    print " Finding strong spots"
    print "-" * 80
    print ""

    # Instantiate the threshold function
    threshold = XDSThresholdStrategy()

    # Get the raw data and image mask
    image = self.experiment.imageset.get_raw_data(0)[0]
    mask = self.experiment.imageset.get_mask(0)[0]

    # Threshold the image and create the pixel labeller
    threshold_mask = threshold(image, mask)
    pixel_labeller = PixelListLabeller()
    pixel_labeller.add(PixelList(0, image, threshold_mask))

    # Create the shoebox list from the pixel list
    creator = flex.PixelListShoeboxCreator(
      pixel_labeller,
      0,                                       # Panel number
      0,                                       # Z start
      True,                                    # 2D
      self.params.spot_finding.min_spot_size,  # Min Pixels
      self.params.spot_finding.max_spot_size,  # Max Pixels
      False)                                   # Find hot pixels
    shoeboxes = creator.result()

    # Filter the list to remove large and small spots
    size = creator.spot_size()
    large = size > self.params.spot_finding.max_spot_size
    small = size < self.params.spot_finding.min_spot_size
    bad = large | small
    shoeboxes = shoeboxes.select(~bad)
    print "Discarding %d spots with < %d pixels" % (
      small.count(True),
      self.params.spot_finding.min_spot_size)
    print "Discarding %d spots with > %d pixels" % (
      large.count(True),
      self.params.spot_finding.max_spot_size)

    # Extract the strong spot information
    centroid = shoeboxes.centroid_valid()
    intensity = shoeboxes.summed_intensity()
    observed = flex.observation(shoeboxes.panels(), centroid, intensity)

    # Create the reflection list
    self.reflections = flex.reflection_table(observed, shoeboxes)
    print "Using %d strong spots" % len(self.reflections)

  def index_spots(self):
    '''
    Assign miller indices to the strong spots

    '''
    from dials.array_family import flex
    from scitbx import matrix
    from math import sqrt

    print ""
    print "-" * 80
    print " Indexing strong spots"
    print "-" * 80
    print ""

    # Initialise the arrays
    miller_index = flex.miller_index()
    distance = flex.double()

    # Get the models
    beam = self.experiment.beam
    panel = self.experiment.detector[0]
    crystal = self.experiment.crystal

    # Get some stuff we need for calculations
    s0 = matrix.col(beam.get_s0())
    A = matrix.sqr(crystal.get_A())
    Ainv = A.inverse()

    # For each reflection get the index and distance
    for index, refl in enumerate(self.reflections):
      x, y, z = refl['xyzobs.px.value']
      s1 = matrix.col(panel.get_pixel_lab_coord((x, y)))
      r = s1.normalize() * s0.length() - s0
      h = Ainv * r
      h0 = tuple(int(round(a)) for a in h)
      d = sqrt(sum(map(lambda a,b: (a-b)**2, h, h0)))
      miller_index.append(h0)
      distance.append(d)

    # Assign the miller indices
    self.reflections['miller_index'] = miller_index

    # Select only those spots within the tolerance
    selection = distance < self.params.indexing.tolerance
    self.reflections = self.reflections.select(selection)
    print "Indexed %d strong spots within tolerance %f" % (
      len(self.reflections),
      self.params.indexing.tolerance)

  def sample_reflections(self):
    '''
    Select a sample of reflections

    '''
    from dials.array_family import flex
    from random import sample
    if self.params.refinement.sample:
      print "Selecting %d/%d reflections for refinement" % (
        self.params.refinement.sample,
        len(self.reflections))
      selection = flex.size_t(
        sorted(sample(
          range(len(self.reflections)),
          self.params.refinement.sample)))
      self.reflections = self.reflections.select(selection)

  def refine_profile(self):
    '''
    Refine the profile parameters

    '''
    print ""
    print "-" * 80
    print " Refining profile parameters"
    print "-" * 80
    print ""

    refiner = ProfileRefiner(
      self.experiment,
      self.reflections,
      self.params)
    self.experiment = refiner.experiment
    self.reflections = refiner.reflections
    self.mosaicity = refiner.profile

    # from matplotlib import pylab
    # X1, X2, Y1, Y2 = zip(*refiner.history)
    # X = X1 + X2
    # Y = Y1 + Y2
    # pylab.scatter(X, Y)
    # pylab.show()

  def refine_crystal(self):
    '''
    Refine the crystal parameters

    '''
    print ""
    print "-" * 80
    print " Refining crystal parameters"
    print "-" * 80
    print ""

    refiner = CrystalRefiner(
      self.experiment,
      self.reflections,
      self.mosaicity,
      self.params)
    self.experiment = refiner.experiment
    self.reflections = refiner.reflections

  def integrate(self):
    '''
    Predict and integrate the reflections

    '''
    from dials_scratch.jmp.stills import Model
    from dials.array_family import flex

    print ""
    print "-" * 80
    print " Integrating reflections"
    print "-" * 80
    print ""

    # Get the data and image mask
    data = self.experiment.imageset.get_raw_data(0)[0]
    mask = self.experiment.imageset.get_mask(0)[0]

    # Initialise the model
    model = Model(
      beam             = self.experiment.beam,
      detector         = self.experiment.detector,
      crystal          = self.experiment.crystal,
      reflections      = flex.reflection_table(),
      image_data       = data,
      image_mask       = mask,
      mosaicity        = self.mosaicity,
      bandpass         = 0.0,
      foreground_limit = 0.3,
      background_limit = 0.5,
      num_samples      = self.params.refinement.profile.num_samples,
      predict_all      = True)

    # Update the model
    model.update()

    # Get the predicted image
    self.image_pred = model.image_pred

    # Get the predicted reflections
    self.reflections = model.reflections

    # Add the columns of data
    self.reflections['id'] = flex.size_t(len(self.reflections), 0)
    self.reflections['partial_id'] = flex.size_t(range(len(self.reflections)))
    self.reflections['panel'] = flex.size_t(len(self.reflections), 0)
    self.reflections['intensity.sum.value'] = model.intensity()
    self.reflections['intensity.sum.variance'] = model.variance()
    self.reflections['xyzcal.px'] = model.predicted()
    self.reflections['xyzobs.px'] = model.observed()
    self.reflections['partiality'] = model.scale()
    self.reflections['num_pixels.background'] = model.num_background()
    self.reflections['num_pixels.foreground'] = model.num_foreground()
    self.reflections['shoebox'] = model.shoebox()

    # Set the integrated flag
    self.reflections.set_flags(model.success(), self.reflections.flags.integrated_sum)

    # Select only those reflections actually predicted
    # i.e. num_fg and num_bg > 0
    min_foreground_pixels = self.params.integration.prediction.min_foreground_pixels
    min_background_pixels = self.params.integration.prediction.min_background_pixels
    selection1 = self.reflections['num_pixels.background'] > min_background_pixels
    selection2 = self.reflections['num_pixels.foreground'] > min_foreground_pixels
    selection = selection1 & selection2
    self.reflections = self.reflections.select(selection)
    print "Integrated %d reflections" % len(self.reflections)


if __name__ == '__main__':

  from dxtbx.model.experiment_list import ExperimentListFactory
  from dxtbx.model.experiment_list import ExperimentListDumper
  from dials.array_family import flex
  import sys

  experiments_filename = sys.argv[1]

  # Read the experiments
  experiments = ExperimentListFactory.from_json_file(experiments_filename)

  params = phil_scope.fetch(parse("")).extract()

  integrator = Integrator(experiments[0], params)
  integrator.process()

  reflections = integrator.reflections
  experiments[0] = integrator.experiment

  selection = reflections.get_flags(reflections.flags.integrated_sum)
  partiality = reflections['partiality'].select(selection)
  min_partiality = flex.min(partiality)
  max_partiality = flex.max(partiality)

  print ""
  print "Mosaicity: %f" % integrator.mosaicity
  print "Min partiality: %f, Max partiality: %f" % (
    min_partiality, max_partiality)
  print ""

  from matplotlib import pylab
  pylab.hist(partiality, bins=100)
  pylab.show()

  reflections.as_pickle("integrated.pickle")
  ExperimentListDumper(experiments).as_json("integrated_experiments.json")
