#!/usr/bin/env dials.python

from libtbx.phil import parse
from libtbx.utils import Sorry
from dials.array_family import flex
from math import pi, sqrt, floor

phil_scope = parse('''

  random_seed = 42
    .type = int

  centroid_sigmas
  {
    sigX_px = 0.25
      .help = "sigma for random normal error to be added to the observed"
              "centroid in pixels along the fast axis direction"
      .type = float(value_min=0)

    sigY_px = 0.25
      .help = "sigma for random normal error to be added to the observed"
              "centroid in pixels along the slow axis direction"
      .type = float(value_min=0)

    sigZ_px = 0.25
      .help = "sigma for random normal error to be added to the observed"
              "centroid in images along the rotation axis direction"
      .type = float(value_min=0)
  }
''')

help_message = '''

Simulated indexed observations for refinement testing. Predictions will be
made for all reflections in the input reflections .pickle file, for each of the
experimental geometries described by the input experiments .json files.

Examples::

  dials.python create_indexed.py observed.pickle experiments1.json experiments2.json

'''

# constants
TWO_PI = 2.0 * pi
RAD_TO_DEG = 180. / pi

class Script(object):
  '''Class to run script.'''

  def __init__(self):
    '''Setup the script.'''
    from dials.util.options import OptionParser
    import libtbx.load_env

    usage = "usage: %s [options] experiments.json" % libtbx.env.dispatcher_name
    self.parser = OptionParser(
      usage=usage,
      phil=phil_scope,
      read_experiments=True,
      read_reflections=True,
      check_format=False,
      epilog=help_message)

  def run(self):

    self.params, options = self.parser.parse_args()
    if len(self.params.input.experiments) == 0:
      self.parser.print_help()
      raise Sorry("No experiments found in the input")
    if len(self.params.input.reflections) == 0:
      self.parser.print_help()
      raise Sorry("No reflection data found in the input")

    import random
    random.seed(self.params.random_seed)
    flex.set_random_seed(self.params.random_seed)

    from dials.util.options import flatten_reflections
    self.original_reflections = flatten_reflections(self.params.input.reflections)[0]

    print "Number of reflections loaded: {0}".format(len(self.original_reflections))

    # make errors to add to observed centroids
    self.set_random_error()

    for wrapper in self.params.input.experiments:
      self.create_indexed(experiments=wrapper.data, filename=wrapper.filename)

  def set_random_error(self):

    import random
    nref = len(self.original_reflections)

    self.shiftX_px = flex.double([random.gauss(0,
        self.params.centroid_sigmas.sigX_px) for _ in xrange(nref)])
    self.shiftY_px = flex.double([random.gauss(0,
        self.params.centroid_sigmas.sigY_px) for _ in xrange(nref)])
    self.shiftZ_px = flex.double([random.gauss(0,
        self.params.centroid_sigmas.sigZ_px) for _ in xrange(nref)])
    return

  def create_indexed(self, experiments, filename):

    print "Simulating indexed observations for {0}".format(filename)

    # Extract the experiment
    exp=experiments[0]

    #print "Experiment scan range is {0},{1}".format(*exp.scan.get_oscillation_range(deg=True))

    # Copy essential columns of the original reflections into a new table
    # (imageset_id is required for reciprocal_lattice_viewer)
    obs_refs = flex.reflection_table()
    cols = ['id', 'imageset_id', 'miller_index', 'panel', 's1', 'flags', 'entering',
            'xyzobs.mm.value', 'xyzobs.mm.variance', 'xyzcal.mm',
            'xyzobs.px.value', 'xyzobs.px.variance', 'xyzcal.px']
    for k in cols:
      if k in self.original_reflections.keys():
        obs_refs[k] = self.original_reflections[k]

    x_obs, y_obs, phi_obs = obs_refs['xyzobs.mm.value'].parts()

    #print "Original observation scan range is {0},{1}".format(
    #    flex.min(phi_obs) * RAD_TO_DEG,
    #    flex.max(phi_obs) *  RAD_TO_DEG)

    # Reset panel number to zero prior to prediction
    obs_refs['panel'] = obs_refs['panel'] * 0

    # Predict new centroid positions
    from dials.algorithms.refinement.prediction import ExperimentsPredictor
    predictor=ExperimentsPredictor(experiments)
    predictor(obs_refs)

    # Get vectors to add as errors to the predicted centroid positions to form
    # new 'observations'
    shift_px = flex.vec3_double(self.shiftX_px,
                                self.shiftY_px,
                                self.shiftZ_px)
    px_size_mm = exp.detector[0].get_pixel_size()
    image_width_rad = exp.scan.get_oscillation(deg=False)[1]
    shift = flex.vec3_double(px_size_mm[0] * self.shiftX_px,
                             px_size_mm[1] * self.shiftY_px,
                             image_width_rad * self.shiftZ_px)
    obs_refs['xyzobs.px.value'] = obs_refs['xyzcal.px'] + shift_px
    obs_refs['xyzobs.mm.value'] = obs_refs['xyzcal.mm'] + shift

    x, y, z = obs_refs['xyzobs.mm.value'].parts()

    #print "Simulated observation scan range is {0},{1}".format(
    #    flex.min(z) * RAD_TO_DEG,
    #    flex.max(z) * RAD_TO_DEG)

    # Keep original variance estimates from spot-finding for centroids in
    # pixels/images, but rescale for centroids in mm/rad for this experiment
    var_x_px, var_y_px, var_z_px = obs_refs['xyzobs.px.variance'].parts()
    var_x_mm = var_x_px * px_size_mm[0]**2
    var_y_mm = var_y_px * px_size_mm[1]**2
    var_z_rd = var_z_px * image_width_rad**2
    obs_refs['xyzobs.mm.variance'] = flex.vec3_double(var_x_mm, var_y_mm, var_z_rd)

    # Reset observed s1 vectors
    from dials.algorithms.refinement.refinement_helpers import set_obs_s1
    obs_refs['s1'] = flex.vec3_double(len(obs_refs))
    set_obs_s1(obs_refs, experiments)

    import os
    outname = os.path.splitext(filename)[0] + '.pickle'
    print "Saving reflections to {0}".format(outname)
    obs_refs.as_pickle(outname)
    print

if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
