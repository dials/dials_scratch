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
      limit_resolution_by = *n_spots d_min
        .type = choice

      d_min = 15.0
        .type = float(value_min=0)

      n_spots = 10
        .type = int

      dstar_tolerance = 4.0
        .help = "Number of sigmas from the centroid position for which to"
                "calculate d* bands"
        .type = float
    }
    debug_reflections = None
      .type = path

    use_P1_indices_as_seeds = False
      .type = bool

    search_depth = *triplets quads
      .type = choice

    bootstrap_crystal = False
      .type = bool
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

class CompleteGraph(object):

  def __init__(self, seed_vertex):

    self.vertices = [seed_vertex]
    self.weight = [{0:0.0}]
    self.total_weight = 0.0

  def factory_add_vertex(self, vertex, weights_to_other):
    # Return a new graph as a copy of this with an extra vertex added. This
    # is a factory rather than a change in-place because CompleteGraph ought
    # to be immutable to implement __hash__
    g = copy.deepcopy(self)

    current_len = len(g.vertices)
    assert len(weights_to_other) == current_len
    g.vertices.append(vertex)
    node = current_len

    # Update distances from other nodes to the new one
    for i, w in enumerate(weights_to_other):
      g.weight[i][node] = w

    # Add distances to other nodes from this one
    weights_to_other.append(0.0)
    to_other = {}
    for i, w in enumerate(weights_to_other):
      to_other[i] = w
    g.weight.append(to_other)

    # Update the total weight
    g.total_weight += sum(weights_to_other)

    # Sort the vertices and weights by spot_id
    l = zip(g.vertices, g.weight)
    l.sort(key=lambda v_w: v_w[0]['spot_id'])
    v, w = zip(*l)
    g.vertices = [e for e in v]
    g.weight = [e for e in w]

    return g

  def __hash__(self):
    h = tuple((e['spot_id'], e['miller_index']) for e in self.vertices)
    return hash(h)

  def __eq__(self, other):
    for a, b in zip(self.vertices, other.vertices):
      if a['spot_id'] != b['spot_id']:
        return False
      if a['miller_index'] != b['miller_index']:
        return False
    return True

  def __ne__(self, other):
    return not self == other

from dials.algorithms.indexing.indexer import indexer_base
class indexer_low_res_spot_match(indexer_base):

  def __init__(self, reflections, imagesets, params):
    super(indexer_low_res_spot_match, self).__init__(reflections, imagesets, params)

  def debug(self):
    # load up known indices
    from dials.array_family import flex
    indexed = flex.reflection_table.from_pickle(
        self.params.low_res_spot_match.debug_reflections)
    idx_dstar = indexed['rlp'].norms()
    dstar_max = 1./(self.params.low_res_spot_match.candidate_spots.d_min - 0.5)
    sel = idx_dstar <= dstar_max
    indexed = indexed.select(sel)
    return indexed

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

    # Take a subset of the observations at the same resolution and calculate
    # some values that will be needed for the search
    self._calc_obs_data()

    # Construct a library of candidate low res indices with their d* values
    self._calc_candidate_hkls()

    # First search: match each observation with candidate indices within the
    # acceptable resolution band
    self._calc_seeds_and_stems()

    # Second search: match seed spots with another spot from a different
    # reciprocal lattice row, such that the observed reciprocal space distances
    # are within tolerances
    stems = []
    if self.params.low_res_spot_match.use_P1_indices_as_seeds:
      seeds = self.stems
    else:
      seeds = self.seeds
    for seed in seeds:
      stems.extend(self._pairs_with_seed(seed))
    stems = list(set(stems)) # filter duplicates

    # Further search iterations: extend to more spots within tolerated distances
    triplets = []
    for stem in stems:
      triplets.extend(self._extend_by_candidates(stem))
    triplets = list(set(triplets)) # filter duplicates

    if self.params.low_res_spot_match.debug_reflections:
      idx = self.debug()
      # add known indices to the correct spots
      hkl = flex.miller_index(len(self.spots))
      for ref in idx:
        for ispot, spot in enumerate(self.spots):
          if ref['xyzobs.mm.value'] == spot['xyzobs.mm.value']:
            hkl[ispot] = ref['miller_index']
      self.spots['known_hkl'] = hkl

    branches = triplets
    if self.params.low_res_spot_match.search_depth == 'quads':
      quads = []
      for triplet in triplets:
        #quads.extend(self._extend_by_candidates(triplet, debug=True))
        quads.extend(self._extend_by_candidates(triplet))
      quads = list(set(quads)) # filter duplicates
      branches = quads

    # Sort branches by total deviation of observed distances from expected
    branches.sort(key=operator.attrgetter('total_weight'))

    candidate_crystal_models = []
    for branch in branches:
      model = self._fit_crystal_model(branch)
      if model:
        candidate_crystal_models.append(model)
      if len(candidate_crystal_models) == self.params.basis_vector_combinations.max_refine:
        break

    # sort models by rmsd
    #candidate_crystal_models.sort(key=operator.attrgetter('rms'))

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

    spot_dstar = self.reflections['rlp'].norms()
    if self.params.low_res_spot_match.candidate_spots.limit_resolution_by == "n_spots":
      n_spots = self.params.low_res_spot_match.candidate_spots.n_spots
      n_spots = min(n_spots, len(self.reflections) - 1)
      dstar_max = flex.sorted(spot_dstar)[n_spots - 1]
      self.params.low_res_spot_match.candidate_spots.d_min = 1./dstar_max

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
    self.spots['dstar_band'] = self.spots['dstar_outer'] - self.spots['dstar_inner']

    return

  def _calc_seeds_and_stems(self):
    # As the first stage of search, determine a list of seed spots for further
    # stages. Order these by distance of observed d* from the candidate
    # reflection's canonical d*

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
                       'rlp_datum': self.Bmat * c['miller_index'],
                       'residual_dstar':r_dst,
                       'clock_angle':spot['clock_angle']})

    result.sort(key=operator.itemgetter('residual_dstar'))
    self.seeds = result

    # Now the 'stems' to use in second search level, using all indices in P 1
    result = []
    for i, spot in enumerate(self.spots):
      sel = ((self.candidate_hkls_p1['dstar'] <= spot['dstar_outer']) &
             (self.candidate_hkls_p1['dstar'] >= spot['dstar_inner']))
      cands = self.candidate_hkls_p1.select(sel)
      for c in cands:
        r_dst = abs(c['dstar'] - spot['dstar'])
        result.append({'spot_id':i,
                       'miller_index':c['miller_index'],
                       'rlp_datum': self.Bmat * c['miller_index'],
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
      seed_vec = seed['rlp_datum']
      cand_vec = cand['rlp_datum']
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
      band1 = self.spots[seed['spot_id']]['dstar_band']
      band2 = self.spots[cand['spot_id']]['dstar_band']
      if r_dist > band1 + band2:
        continue

      # Store the seed-stem match as a 2-node graph
      g = CompleteGraph({'spot_id':seed['spot_id'],
                         'miller_index':seed['miller_index'],
                         'rlp_datum':seed['rlp_datum']})
      g = g.factory_add_vertex({'spot_id':cand['spot_id'],
                                'miller_index':cand['miller_index'],
                                'rlp_datum':cand['rlp_datum']},
                                weights_to_other=[r_dist])

      # FIXME at the moment, monkey patching the graph object to contain extra
      # data about the plane between relps. Later move to a co-planarity check
      # within _extend_by_candidates, in which planarity of all vertices
      # in the graph is checked
      g.plane_normal = seed_vec.cross(cand_vec).normalize()

      result.append(g)

    result.sort(key=operator.attrgetter('total_weight'))
    return result

  def _extend_by_candidates(self, graph):

    existing_ids = [e['spot_id'] for e in graph.vertices]
    obs_relps = [matrix.col(self.spots[e]['rlp']) for e in existing_ids]
    exp_relps = [e['rlp_datum'] for e in graph.vertices]

    result = []

    for cand in self.stems:
      # Don't check spots already matched
      if cand['spot_id'] in existing_ids:
        continue

      # Compare expected reciprocal space distances with observed distances
      cand_rlp = matrix.col(self.spots[cand['spot_id']]['rlp'])
      cand_vec = cand['rlp_datum']

      obs_dists = [(cand_rlp - rlp).length() for rlp in obs_relps]
      exp_dists = [(vec - cand_vec).length() for vec in exp_relps]

      residual_dist = [abs(a - b) for (a, b) in zip(obs_dists, exp_dists)]

      # If any of the distance differences is larger than the sum of the
      # tolerated d* bands then reject the candidate
      candidate_band = self.spots[cand['spot_id']]['dstar_band']
      bad_candidate = False
      for r_dist, spot_id in zip(residual_dist, existing_ids):
        relp_band = self.spots[spot_id]['dstar_band']
        if r_dist > relp_band + candidate_band:
          bad_candidate = True
          break
      if bad_candidate:
        continue

      # Calculate co-planarity of the candidate. If plane_score is too high
      # (corresponding to the relp being more than two degrees off the plane)
      # then reject the candidate.
      # FIXME implement better plane score
      plane_score = abs(cand_vec.normalize().dot(graph.plane_normal))
      if plane_score > 0.035:
        continue

      # Construct a graph including the candidate node
      g = graph.factory_add_vertex({'spot_id':cand['spot_id'],
                                    'miller_index':cand['miller_index'],
                                    'rlp_datum':cand['rlp_datum']},
                                    weights_to_other=residual_dist)

      # FIXME at the moment, monkey patching the graph object to contain extra
      # data about the plane score. It may be better to calculate this
      # externally once the requested number of search iterations are complete
      g.plane_score = plane_score

      result.append(g)

    return result

  def _fit_crystal_model(self, graph):

    vertices = graph.vertices

    # Reciprocal lattice points of the observations
    sel=flex.size_t([e['spot_id'] for e in vertices])
    reference = self.spots['rlp'].select(sel)

    # Ideal relps from the known cell
    other = flex.vec3_double([e['rlp_datum'] for e in vertices])

    # Add the origin to both sets of points
    origin = flex.vec3_double(1)
    reference.extend(origin)
    other.extend(origin)

    # Find U matrix that takes ideal relps to the reference
    fit = superpose.least_squares_fit(reference, other)

    # Calculate RMSD of the fit
    other_rotated = fit.r.elems * other
    rms = reference.rms_difference(other_rotated)

    # Construct a crystal model
    UB = fit.r * self.Bmat
    xl = Crystal(A=UB, space_group_symbol="P1")

    if self.params.low_res_spot_match.bootstrap_crystal:

      # attempt to index the low resolution spots
      from dials_algorithms_indexing_ext import AssignIndices
      phi = self.spots['xyzobs.mm.value'].parts()[2]
      UB_matrices = flex.mat3_double([xl.get_A(),])
      result = AssignIndices(self.spots['rlp'], phi, UB_matrices, tolerance=0.3)
      hkl = result.miller_indices()
      sel = hkl != (0,0,0)
      hkl_vec = hkl.as_vec3_double().select(sel)

      reference = self.spots['rlp'].select(sel)
      other = self.Bmat.elems * hkl_vec

      # Add origin
      reference.extend(origin)
      other.extend(origin)

      # Find U matrix that takes ideal relps to the reference
      fit = superpose.least_squares_fit(reference, other)

      # Calculate RMSD of the fit
      other_rotated = fit.r.elems * other
      rms = reference.rms_difference(other_rotated)

      # Construct a crystal model
      UB = fit.r * self.Bmat
      xl = Crystal(A=UB, space_group_symbol="P1")

    xl.rms = rms

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
