from __future__ import absolute_import, division

import logging
logger = logging.getLogger('indigo')

#try:
#  # try importing scipy.linalg before any cctbx modules, otherwise we
#  # sometimes get a segmentation fault/core dump if it is imported after
#  # scipy.linalg is a dependency of sklearn.cluster.DBSCAN
#  import scipy.linalg # import dependency
#except ImportError:
#  pass

import copy

import libtbx
from libtbx.utils import Sorry
import iotbx.phil
from dials.util.options import OptionParser
from dials.util.options import flatten_reflections
from dials.util.options import flatten_datablocks
from dials.util.options import flatten_experiments
import cctbx.miller
from dials.array_family import flex
import operator

help_message = '''

Experimental indexing algorithm for electron diffraction still shots. Requires a
known unit cell, low order diffraction spots and good geometry calibration

Example::

  dials.python indigo.py datablock.json strong.pickle unit_cell=79,79,37,90,90,90 space_group=P43212

'''

phil_scope = iotbx.phil.parse("""\
include scope dials.algorithms.indexing.indexer.master_phil_scope

indexing{
  low_res_spot_match
    .expert_level = 1
  {
    candidate_spots
    {
      d_min = 16.0
        .type = float(value_min=0)

      dstar_tolerance = 5.0
        .help = "Number of sigmas from the centroid position for which to"
                "calculate d* bands"
        .type = float
    }
  }
}

output {
  experiments = experiments.json
    .type = path
  split_experiments = False
    .type = bool
  reflections = indexed.pickle
    .type = path
  unindexed_reflections = None
    .type = path
  log = indigo.log
    .type = str
  debug_log = indigo.debug.log
    .type = str
}

verbosity = 1
  .type = int(value_min=0)
  .help = "The verbosity level"
""", process_includes=True)

# local overrides for refiner.phil_scope
phil_overrides = iotbx.phil.parse('''
refinement
{
  verbosity = 1
}
''')

working_phil = phil_scope.fetch(sources=[phil_overrides])
master_params = working_phil.extract()


from dials.algorithms.indexing.indexer import indexer_base
class indexer_low_res_spot_match(indexer_base):

  def __init__(self, reflections, imagesets, params):
    super(indexer_low_res_spot_match, self).__init__(reflections, imagesets, params)

  @staticmethod
  def from_parameters(reflections, imagesets,
                      known_crystal_models=None, params=None):

    if params is None:
      params = master_params

    if known_crystal_models is not None:
      from dials.algorithms.indexing.known_orientation \
           import indexer_known_orientation
      if params.indexing.known_symmetry.space_group is None:
        params.indexing.known_symmetry.space_group \
          = known_crystal_models[0].get_space_group().info()
      idxr = indexer_known_orientation(
        reflections, imagesets, params, known_crystal_models)
    else:
      has_stills = False
      has_sweeps = False
      for imageset in imagesets:
        if imageset.get_goniometer() is None or imageset.get_scan() is None or \
            imageset.get_scan().get_oscillation()[1] == 0:
          if has_sweeps:
            raise Sorry("Please provide only stills or only sweeps, not both")
          has_stills = True
        else:
          if has_stills:
            raise Sorry("Please provide only stills or only sweeps, not both")
          has_sweeps = True
      assert not (has_stills and has_sweeps)
      use_stills_indexer = has_stills

      if not (params.indexing.stills.indexer is libtbx.Auto or params.indexing.stills.indexer.lower() == 'auto'):
        if params.indexing.stills.indexer == 'stills':
          use_stills_indexer = True
        elif params.indexing.stills.indexer == 'sweeps':
          use_stills_indexer = False
        else:
          assert False

      if params.indexing.basis_vector_combinations.max_refine is libtbx.Auto:
        if use_stills_indexer:
          params.indexing.basis_vector_combinations.max_refine = 5
        else:
          params.indexing.basis_vector_combinations.max_refine = 50

      if use_stills_indexer:
        # Ensure the indexer and downstream applications treat this as set of stills
        from dxtbx.imageset import ImageSet#, MemImageSet
        reset_sets = []
        for i in xrange(len(imagesets)):
          imagesweep = imagesets.pop(0)
          imageset = ImageSet(imagesweep.data(), imagesweep.indices())
          # if isinstance(imageset, MemImageSet):
          #   imageset = MemImageSet(imagesweep._images, imagesweep.indices())
          # else:
          #   imageset = ImageSet(imagesweep.reader(), imagesweep.indices())
          #   imageset._models = imagesweep._models
          imageset.set_scan(None)
          imageset.set_goniometer(None)
          reset_sets.append(imageset)
        imagesets.extend(reset_sets)

      # Parameters are set differently depending on use_stills_indexer, but
      # either way you get an indexer_low_res_spot_match
      idxr = indexer_low_res_spot_match(reflections, imagesets, params=params)

    return idxr

  def find_lattices(self):

    try:
      assert self.target_symmetry_primitive is not None
      assert self.target_symmetry_primitive.unit_cell() is not None
    except AssertionError:
      raise Sorry("indigo requires a known unit_cell=a,b,c,aa,bb,cc")

    self._low_res_spot_match()
    crystal_models = self.candidate_crystal_models
    experiments = ExperimentList()
    for cm in crystal_models:
      for imageset in self.imagesets:
        experiments.append(Experiment(imageset=imageset,
                                      beam=imageset.get_beam(),
                                      detector=imageset.get_detector(),
                                      goniometer=imageset.get_goniometer(),
                                      scan=imageset.get_scan(),
                                      crystal=cm))
    return experiments

  def _low_res_spot_match(self):

    # Construct a library of candidate low res indices with their d* values
    self.candidate_hkls = self._calc_candidate_hkls(
        self.target_symmetry_primitive,
        self.params.low_res_spot_match.candidate_spots.d_min)

    # Take a subset of the observations at the same resolution
    spot_dstar = self.reflections['rlp'].norms()
    dstar_max = 1./self.params.low_res_spot_match.candidate_spots.d_min
    sel = spot_dstar <= dstar_max
    self.spots = self.reflections.select(sel)
    self.spots['dstar'] = spot_dstar.select(sel)

    # Calculate d* bands encompassing each spot
    self._calc_dstar_bands()

    # First search: match each observation with candidate indices within the
    # resolution band
    seeds = self._seed_spots()

    self.candidate_crystal_models = None
    raise NotImplementedError("Nothing to see here")

  @staticmethod
  def _calc_candidate_hkls(symmetry, d_min):
    hkl_list = cctbx.miller.build_set(symmetry, anomalous_flag=False,
      d_min=d_min)
    rt = flex.reflection_table()
    rt['miller_index'] = hkl_list.indices()
    rt['dstar'] = 1. / hkl_list.d_spacings().data()
    return rt

  def _calc_dstar_bands(self):
    # Convert sigma_x and sigma_y in pixels for each reflection into sigma_r
    # in radial direction from the beam centre

    # XXX In what circumstance might there be more than one imageset?
    detector = self.imagesets[0].get_detector()
    beam = self.imagesets[0].get_beam()

    beam_centre = flex.vec2_double(len(self.spots))
    pnl_ids = set(self.spots['panel'])
    for pnl in pnl_ids:
      sel = self.spots['panel'] == pnl
      panel = detector[pnl]
      bc = panel.get_ray_intersection_px(beam.get_s0())
      beam_centre.set_selected(sel, bc)

    # Lab coordinate of the beam centre
    panel = detector[self.spots[0]['panel']]
    bc = panel.get_ray_intersection(beam.get_s0())
    bc_lab = panel.get_lab_coord(bc)

    # Lab coordinate of each spot
    spot_lab = flex.vec3_double(len(self.spots))
    pnl_ids = set(self.spots['panel'])
    for pnl in pnl_ids:
      sel = self.spots['panel'] == pnl
      panel = detector[pnl]
      obs = self.spots['xyzobs.mm.value'].select(sel)
      x_mm, y_mm, _ = obs.parts()
      spot_lab.set_selected(sel, panel.get_lab_coord(flex.vec2_double(x_mm, y_mm)))

    # Radius and radial direction vectors for each spot
    radius = spot_lab - bc_lab
    unit_r = radius.each_normalize()

    # Calculate radius and direction vector for each spot from the beam centre
    x, y, _ = radius.parts()
    x2, y2 = flex.pow2(x), flex.pow2(y)
    r2 = x2 + y2

    # Calc error along unit_radius direction with simple error propagation
    # that assumes no covariance between x and y centroid errors
    sig_x2, sig_y2, _ = self.spots['xyzobs.mm.variance'].parts()
    var_r = (x2 / r2) * sig_x2 + (y2 / r2) * sig_y2
    sig_r = flex.sqrt(var_r)

    # Pixel coordinates at limits of the band
    tol = self.params.low_res_spot_match.candidate_spots.dstar_tolerance
    outer_spot_lab = spot_lab + unit_r * (tol * sig_r)
    inner_spot_lab = spot_lab - unit_r * (tol * sig_r)

    # Set d* at band limits
    inv_lambda = 1./beam.get_wavelength()
    s1_outer = outer_spot_lab.each_normalize() * inv_lambda
    s1_inner = inner_spot_lab.each_normalize() * inv_lambda
    self.spots['dstar_outer'] = (s1_outer - beam.get_s0()).norms()
    self.spots['dstar_inner'] = (s1_inner - beam.get_s0()).norms()

    return

  def _seed_spots(self):
    # As the first stage of search, determine a list of seed spots for further
    # stages. Order these by distance of observation d* from the candidate
    # reflection d*

    result = []
    for i, spot in enumerate(self.spots):
      sel = ((self.candidate_hkls['dstar'] <= spot['dstar_outer']) &
             (self.candidate_hkls['dstar'] >= spot['dstar_inner']))
      cands = self.candidate_hkls.select(sel)
      for c in cands:
        dist = abs(c['dstar'] - spot['dstar'])
        result.append({'spot_id':i,
                       'miller_index':c['miller_index'],
                       'dist':dist})

    result.sort(key=operator.itemgetter('dist'))
    return result


def run(args):
  import libtbx.load_env
  from libtbx.utils import Sorry
  from dials.util import log
  usage = "%s [options] datablock.json strong.pickle" %libtbx.env.dispatcher_name

  parser = OptionParser(
    usage=usage,
    phil=working_phil,
    read_reflections=True,
    read_datablocks=True,
    read_experiments=True,
    check_format=False,
    epilog=help_message)

  params, options = parser.parse_args(show_diff_phil=False)

  # Configure the logging
  log.config(
    params.verbosity,
    info=params.output.log,
    debug=params.output.debug_log)

  from dials.util.version import dials_version
  logger.info(dials_version())

  # Log the diff phil
  diff_phil = parser.diff_phil.as_str()
  if diff_phil is not '':
    logger.info('The following parameters have been modified:\n')
    logger.info(diff_phil)

  datablocks = flatten_datablocks(params.input.datablock)
  experiments = flatten_experiments(params.input.experiments)
  reflections = flatten_reflections(params.input.reflections)

  if len(datablocks) == 0:
    if len(experiments) > 0:
      imagesets = experiments.imagesets()
    else:
      parser.print_help()
      return
  else:
    imagesets = []
    for datablock in datablocks:
      imagesets.extend(datablock.extract_imagesets())
  if len(experiments):
    known_crystal_models = experiments.crystals()
  else:
    known_crystal_models = None

  if len(reflections) == 0:
    raise Sorry("No reflection lists found in input")
  if len(reflections) > 1:
    assert len(reflections) == len(imagesets)
    from scitbx.array_family import flex
    for i in range(len(reflections)):
      reflections[i]['imageset_id'] = flex.int(len(reflections[i]), i)
      if i > 0:
        reflections[0].extend(reflections[i])

  reflections = reflections[0]

  for imageset in imagesets:
    if (imageset.get_goniometer() is not None and
        imageset.get_scan() is not None and
        imageset.get_scan().get_oscillation()[1] == 0):
      imageset.set_goniometer(None)
      imageset.set_scan(None)

  idxr = indexer_low_res_spot_match.from_parameters(
    reflections, imagesets,
    known_crystal_models=known_crystal_models,
    params=params)
  idxr.index()
  refined_experiments = idxr.refined_experiments
  reflections = copy.deepcopy(idxr.refined_reflections)
  reflections.extend(idxr.unindexed_reflections)
  if len(refined_experiments):
    if params.output.split_experiments:
      logger.info("Splitting experiments before output")
      from dxtbx.model.experiment_list import ExperimentList
      refined_experiments = ExperimentList(
        [copy.deepcopy(re) for re in refined_experiments])
    logger.info("Saving refined experiments to %s" % params.output.experiments)
    idxr.export_as_json(refined_experiments,
                        file_name=params.output.experiments)
    logger.info("Saving refined reflections to %s" % params.output.reflections)
    idxr.export_reflections(
      reflections, file_name=params.output.reflections)

    if params.output.unindexed_reflections is not None:
      logger.info("Saving unindexed reflections to %s"
           %params.output.unindexed_reflections)
      idxr.export_reflections(idxr.unindexed_reflections,
                              file_name=params.output.unindexed_reflections)

  return

if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
