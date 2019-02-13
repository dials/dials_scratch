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
from math import pi
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
from scitbx import matrix
from scitbx.math import superpose
from dxtbx.model import Crystal
from dxtbx.model.experiment_list import Experiment, ExperimentList

TWO_PI = 2.0 * pi
FIVE_DEG = TWO_PI * 5./360.

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
      d_min = 15.0
        .type = float(value_min=0)

      dstar_tolerance = 4.0
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

  #DEBUG - remove after algorithm development
  def debug(self):
    # load up known indices
    from dials.array_family import flex
    indexed = flex.reflection_table.from_pickle("indexed.pickle")
    idx_dstar = indexed['rlp'].norms()
    dstar_max = 1./self.params.low_res_spot_match.candidate_spots.d_min
    sel = idx_dstar <= dstar_max
    indexed = indexed.select(sel)
    return indexed

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

    # Set reciprocal space orthogonalisation matrix
    uc = self.target_symmetry_primitive.unit_cell()
    self.Bmat = matrix.sqr(uc.fractionalization_matrix()).transpose()

    self._low_res_spot_match()
    crystal_model, n_indexed = self.choose_best_orientation_matrix(
      self.candidate_crystal_models)
    if crystal_model is not None:
      crystal_models = [crystal_model]
    else:
      crystal_models = []
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
    self._calc_candidate_hkls()

    # Take a subset of the observations at the same resolution and calculate
    # some values that will be needed for the search
    self._calc_obs_data()

    # First search: match each observation with candidate indices within the
    # resolution band
    self._calc_seeds_and_stems()

    # Further search will try to index 3 spots within resolution bands and
    # tolerated reciprocal space distance from one another
    triplets = []
    for seed in self.seeds:
      stems = self._pairs_with_seed(seed)

      for stem in stems:
        branches = self._triplets_with_seed_and_stem(seed, stem)

        for branch in branches:
          triplets.append((seed, stem, branch))

    # DEBUG - check through these results to see which are correct
    idx = self.debug()
    triplets.sort(key=lambda x: x[2]['residual_rlp_dist_total'])

    correct = {0:set((4, 3, 0)), 1:set((1,1,1)), 2:set((1,1,1)), 3:set((4,3,0))}
    print("soln    seed    stem    branch")
    for i, t in enumerate(triplets):
      res = [(e['spot_id'], e['miller_index']) for e in t]
      res.sort(key=operator.itemgetter(0))
      nwrong=0
      for r in res:
        spot = r[0]
        pos_hkl = set([abs(e) for e in r[1]])
        if pos_hkl != correct[spot]: nwrong +=1
      if nwrong > 0: continue
      miller1 = "{}: ({: 2d}, {: 2d}, {: 2d})".format(res[0][0], *res[0][1])
      miller2 = "{}: ({: 2d}, {: 2d}, {: 2d})".format(res[1][0], *res[1][1])
      miller3 = "{}: ({: 2d}, {: 2d}, {: 2d})".format(res[2][0], *res[2][1])
      print("{0:03d} {1} {2} {3}".format(i, miller1, miller2, miller3))

    candidate_crystal_models = []
    for triplet in triplets:
      model = self._fit_crystal_model(triplet)
      if model:
        candidate_crystal_models.append(model)
      if len(candidate_crystal_models) == self.params.basis_vector_combinations.max_refine:
        break

    # At this point, either extract basis vectors from the candidate_crystal_models
    # and pass to indexer_base.find_candidate_orientation_matrices, or just
    # use the candidate_crystal_models as it is. Currently, try the latter
    self.candidate_crystal_models = candidate_crystal_models

  def _calc_candidate_hkls(self):
    # 1 ASU
    hkl_list = cctbx.miller.build_set(self.target_symmetry_primitive,
        anomalous_flag=False,
        d_min=self.params.low_res_spot_match.candidate_spots.d_min)
    rt = flex.reflection_table()
    rt['miller_index'] = hkl_list.indices()
    rt['dstar'] = 1. / hkl_list.d_spacings().data()
    self.candidate_hkls = rt

    # P1 indices with separate Friedel pairs
    hkl_list = cctbx.miller.build_set(self.target_symmetry_primitive,
        anomalous_flag=True,
        d_min=self.params.low_res_spot_match.candidate_spots.d_min)
    hkl_list_p1 = hkl_list.expand_to_p1()
    rt = flex.reflection_table()
    rt['miller_index'] = hkl_list_p1.indices()
    rt['dstar'] = 1. / hkl_list_p1.d_spacings().data()
    self.candidate_hkls_p1 = rt
    return

  def _calc_obs_data(self):
    """Calculates a set of low resolution observations to try to match to
    indices. Each observation will record its d* value as well as
    tolerated d* bands and a 'clock angle'"""

    # First select low resolution spots only
    spot_dstar = self.reflections['rlp'].norms()
    dstar_max = 1./self.params.low_res_spot_match.candidate_spots.d_min
    sel = spot_dstar <= dstar_max
    self.spots = self.reflections.select(sel)
    self.spots['dstar'] = spot_dstar.select(sel)

    # XXX In what circumstance might there be more than one imageset?
    detector = self.imagesets[0].get_detector()
    beam = self.imagesets[0].get_beam()

    # Lab coordinate of the beam centre, using the first spot's panel
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

    # Radius vectors for each spot
    radius = spot_lab - bc_lab

    # Usually the radius vectors would all be in a single plane, but this might
    # not be the case if the spots are on different panels. To put them on the
    # same plane, project onto fast/slow of the panel used to get the beam
    # centre
    df = flex.vec3_double(len(self.spots), detector[0].get_fast_axis())
    ds = flex.vec3_double(len(self.spots), detector[0].get_slow_axis())
    clock_dirs = (radius.dot(df) * df + radius.dot(ds) * ds).each_normalize()

    # From this, find positive angles of each vector around a clock, using the
    # fast axis as 12 o'clock
    angs = clock_dirs.angle(detector[0].get_fast_axis())
    dots = clock_dirs.dot(detector[0].get_slow_axis())
    sel = dots < 0 # select directions in the second half of the clock face
    angs.set_selected(sel, (TWO_PI - angs.select(sel)))
    self.spots['clock_angle'] = angs

    # Project radius vectors onto fast/slow of the relevant panels
    df = flex.vec3_double(len(self.spots))
    ds = flex.vec3_double(len(self.spots))
    for pnl in pnl_ids:
      sel = self.spots['panel'] == pnl
      panel = detector[pnl]
      df.set_selected(sel, panel.get_fast_axis())
      ds.set_selected(sel, panel.get_slow_axis())
    panel_dirs = (radius.dot(df) * df + radius.dot(ds) * ds).each_normalize()

    # Calc error along each panel direction with simple error propagation
    # that assumes no covariance between x and y centroid errors.
    x = panel_dirs.dot(df)
    y = panel_dirs.dot(ds)
    x2, y2 = flex.pow2(x), flex.pow2(y)
    r2 = x2 + y2
    sig_x2, sig_y2, _ = self.spots['xyzobs.mm.variance'].parts()
    var_r = (x2 / r2) * sig_x2 + (y2 / r2) * sig_y2
    sig_r = flex.sqrt(var_r)

    # Pixel coordinates at limits of the band
    tol = self.params.low_res_spot_match.candidate_spots.dstar_tolerance
    outer_spot_lab = spot_lab + panel_dirs * (tol * sig_r)
    inner_spot_lab = spot_lab - panel_dirs * (tol * sig_r)

    # Set d* at band limits
    inv_lambda = 1./beam.get_wavelength()
    s1_outer = outer_spot_lab.each_normalize() * inv_lambda
    s1_inner = inner_spot_lab.each_normalize() * inv_lambda
    self.spots['dstar_outer'] = (s1_outer - beam.get_s0()).norms()
    self.spots['dstar_inner'] = (s1_inner - beam.get_s0()).norms()

    return

  def _calc_seeds_and_stems(self):
    # As the first stage of search, determine a list of seed spots for further
    # stages. Order these by distance of observation d* from the candidate
    # reflection d*

    # First the 'seeds' (in 1 ASU)
    result = []
    for i, spot in enumerate(self.spots):
      sel = ((self.candidate_hkls['dstar'] <= spot['dstar_outer']) &
             (self.candidate_hkls['dstar'] >= spot['dstar_inner']))
      cands = self.candidate_hkls.select(sel)
      for c in cands:
        r_dst = abs(c['dstar'] - spot['dstar'])
        result.append({'spot_id':i,
                       'miller_index':c['miller_index'],
                       'residual_dstar':r_dst,
                       'clock_angle':spot['clock_angle']})

    result.sort(key=operator.itemgetter('residual_dstar'))
    self.seeds = result

    # Now the 'stems' to use in second search level, using all indices
    result = []
    for i, spot in enumerate(self.spots):
      sel = ((self.candidate_hkls_p1['dstar'] <= spot['dstar_outer']) &
             (self.candidate_hkls_p1['dstar'] >= spot['dstar_inner']))
      cands = self.candidate_hkls_p1.select(sel)
      for c in cands:
        r_dst = abs(c['dstar'] - spot['dstar'])
        result.append({'spot_id':i,
                       'miller_index':c['miller_index'],
                       'residual_dstar':r_dst,
                       'clock_angle':spot['clock_angle']})

    result.sort(key=operator.itemgetter('residual_dstar'))
    self.stems = result
    return

  def _pairs_with_seed(self, seed):

    seed_rlp = matrix.col(self.spots[seed['spot_id']]['rlp'])

    result = []
    for cand in self.stems:
      # Don't check the seed spot itself
      if cand['spot_id'] == seed['spot_id']:
        continue

      # Skip spots at a very similar clock angle, which probably belong to the
      # same line of indices from the origin
      angle_diff = cand['clock_angle'] - seed['clock_angle']
      angle_diff = abs(((angle_diff + pi) % TWO_PI) - pi)
      if angle_diff < FIVE_DEG:
        continue

      # Calculate the plane normal for the plane containing the seed and stem.
      # Skip pairs of Miller indices that belong to the same line
      seed_vec = self.Bmat * seed['miller_index']
      cand_vec = self.Bmat * cand['miller_index']
      try:
        seed_vec.cross(cand_vec).normalize()
      except ZeroDivisionError:
        continue

      # Compare expected reciprocal space distance with observed distance
      cand_rlp = matrix.col(self.spots[cand['spot_id']]['rlp'])
      obs_dist = (cand_rlp - seed_rlp).length()
      exp_dist = (seed_vec - cand_vec).length()
      r_dist = abs(obs_dist - exp_dist)

      # If the distance difference is larger than the sum of the tolerated
      # d* bands then reject the candidate
      band1 = (self.spots[seed['spot_id']]['dstar_outer'] -
               self.spots[seed['spot_id']]['dstar_inner'])
      band2 = (self.spots[cand['spot_id']]['dstar_outer'] -
               self.spots[cand['spot_id']]['dstar_inner'])
      assert band1 > 0 # FIXME can remove after algorithm finished
      assert band2 > 0 # FIXME can remove after algorithm finished
      if r_dist > band1 + band2:
        continue

      # copy cand to a new dictionary, include the reciprocal space residual
      # distance and plane normal, then add to result.
      stem = cand.copy()
      stem['residual_rlp_dist'] = r_dist
      stem['plane_normal'] = seed_vec.cross(cand_vec).normalize()
      result.append(stem)

    result.sort(key=operator.itemgetter('residual_rlp_dist'))
    return result

  def _triplets_with_seed_and_stem(self, seed, stem):

    seed_rlp = matrix.col(self.spots[seed['spot_id']]['rlp'])
    stem_rlp = matrix.col(self.spots[stem['spot_id']]['rlp'])

    result = []
    for cand in self.stems:
      # Don't check the seed or first stem spot themselves
      if cand['spot_id'] == seed['spot_id'] or cand['spot_id'] == stem['spot_id']:
        continue

      # Compare expected reciprocal space distances with observed distances
      cand_rlp = matrix.col(self.spots[cand['spot_id']]['rlp'])
      obs_dist1 = (cand_rlp - seed_rlp).length()
      seed_vec = self.Bmat * seed['miller_index']
      cand_vec = self.Bmat * cand['miller_index']
      exp_dist1 = (seed_vec - cand_vec).length()
      r_dist1 = abs(obs_dist1 - exp_dist1)

      obs_dist2 = (cand_rlp - stem_rlp).length()
      stem_vec = self.Bmat * stem['miller_index']
      exp_dist2 = (stem_vec - cand_vec).length()
      r_dist2 = abs(obs_dist2 - exp_dist2)

      # If either of the distance differences is larger than the sum of the
      # tolerated d* bands then reject the candidate
      band1 = (self.spots[seed['spot_id']]['dstar_outer'] -
               self.spots[seed['spot_id']]['dstar_inner'])
      band2 = (self.spots[stem['spot_id']]['dstar_outer'] -
               self.spots[stem['spot_id']]['dstar_inner'])
      band3 = (self.spots[cand['spot_id']]['dstar_outer'] -
               self.spots[cand['spot_id']]['dstar_inner'])
      assert band1 > 0 # FIXME can remove after algorithm finished
      assert band2 > 0 # FIXME can remove after algorithm finished
      assert band3 > 0 # FIXME can remove after algorithm finished
      if r_dist1 > band1 + band3 or r_dist2 > band2 + band3:
        continue

      # Calculate co-planarity of the candidate. If plane_score is too high
      # (corresponding to the relp being more than two degrees off the plane)
      # then reject the candidate
      plane_score = abs(cand_vec.normalize().dot(stem['plane_normal']))
      if plane_score > 0.035:
        continue

      # copy cand to a new dictionary, include the total reciprocal space
      # residual distance as a measure of the how well this candidate matches
      # and add to result. Keep plane_score as well in case this is useful for
      # ranking potential solutions later.
      branch = cand.copy()
      branch['residual_rlp_dist_total'] = (r_dist1 + r_dist2 +
          stem['residual_rlp_dist'])
      branch['plane_score'] = plane_score

      result.append(branch)

    #result.sort(key=operator.itemgetter('residual_rlp_dist_total'))
    return result

  def _fit_crystal_model(self, triplet):

    # Reciprocal lattice points of the observations
    sel=flex.size_t([e['spot_id'] for e in triplet])
    reference = self.spots['rlp'].select(sel)

    # Ideal relps from the known cell
    other = flex.vec3_double([self.Bmat * e['miller_index'] for e in triplet])

    # Find U matrix that takes ideal relps to the reference
    fit = superpose.least_squares_fit(reference, other)

    # Construct a crystal model
    UB = fit.r * self.Bmat
    xl = Crystal(A=UB, space_group_symbol="P1")

    return xl

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
