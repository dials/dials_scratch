from dials.array_family import flex
import cPickle as pickle
import math
import sys
from scitbx import matrix
import scitbx.math as smath
import json

def nint(a):
  return int(round(a))

from scitbx import simplex

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
    self.n = len(values)
    self.x = values
    self.starting_simplex = generate_start(values, offset)
    self.fcount = 0

    optimizer = simplex.simplex_opt(dimension=self.n,
                                    matrix=self.starting_simplex,
                                    evaluator=evaluator,
                                    tolerance=1e-10,
                                    max_iter=max_iter)

    self.x = optimizer.get_solution()
    return

  def get_solution(self):
    return self.x

  def target(self, vector):
    score = scorify(vector)
    return score

class Apple(object):

  def __init__(self, reflection_file, experiment_file):

    data = pickle.load(open(reflection_file, 'rb'))
    print '%d reflections' % data.size()

    self.data = data.select(data['intensity.sum.variance'] > 0)

    i = data['intensity.sum.value']
    v = data['intensity.sum.variance']
    s = flex.sqrt(v)
    self.i_s = i/s

    from dxtbx.model.experiment.experiment_list import ExperimentListFactory
    expt = ExperimentListFactory.from_json_file(experiment_file)
    panel = expt.detectors()[0][0]
    crystal = expt.crystals()[0]
    self.s0 = matrix.col(expt.beams()[0].get_s0())
    wavelength = expt.beams()[0].get_wavelength()

    # make a list of observed q positions

    self.qobs = []
    for j in range(data.size()):
      x, y, z = data['xyzobs.px.value'][j]
      p = matrix.col(panel.get_pixel_lab_coord((x, y)))
      q = p.normalize() / wavelength - self.s0
      self.qobs.append(q)

    self.wavelength = wavelength
    self.panel = panel
    self.beam = expt.beams()[0]

    # slurp data from $somewhere

    imageset = expt.imagesets()[0]
    self.raw_data = imageset.get_raw_data(0)[0]
    self.imageset = imageset

    from dials.algorithms.refinement.parameterisation.crystal_parameters \
      import CrystalUnitCellParameterisation, \
      CrystalOrientationParameterisation

    self.crystal = crystal
    self.cucp = CrystalUnitCellParameterisation(crystal)
    self.cop = CrystalOrientationParameterisation(crystal)

    self.zero()

    # 0-point and deltas
    values = flex.double(self.cucp.get_param_vals() +
                         self.cop.get_param_vals())
    offset = flex.double([0.01 * v for v in self.cucp.get_param_vals()] +
                         [0.1, 0.1, 0.1])

    initial = crystal.get_unit_cell()
    self.cells = []
    self.best_score = 1e99
    initial_score = self.target(values)
    doohicky = simple_simplex(values, offset, self, 2000)
    best = doohicky.get_solution()
    print 'Initial cell:', initial
    print 'Final cell:  ', crystal.get_unit_cell()
    print 'Score change', initial_score, self.target(best, do_print=False)
    self.best = best

  def plot_map(self, map, filename):
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot
    data = map.as_numpy_array()
    fig = pyplot.gcf()
    pyplot.imshow(data, cmap='gray_r')
    pyplot.colorbar()
    pyplot.savefig(filename, dpi=400)
    pyplot.clf()
    return

  def plot_log_map(self, map, filename):
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot

    negative = (map.as_1d() <= 0)
    map.as_1d().set_selected(negative, 1)
    logmap = flex.log10(map.as_double())

    data = logmap.as_numpy_array()
    fig = pyplot.gcf()
    pyplot.imshow(data, cmap='gray_r')
    pyplot.colorbar()
    pyplot.savefig(filename, dpi=400)
    pyplot.clf()
    return

  def render_distance(self):
    distance_map = flex.double(flex.grid(self.raw_data.focus()))
    origin = self.panel.get_origin()
    fast = self.panel.get_fast_axis()
    slow = self.panel.get_slow_axis()
    nfast, nslow = self.panel.get_image_size()

    UB = matrix.sqr(self.crystal.get_A())
    UBi = UB.inverse()

    from dials_scratch import q_map
    distance_map = q_map(self.panel, self.beam, UB, 1)
    return distance_map

  def target(self, vector, do_print=False):
    cell_parms = self.cucp.get_param_vals()
    orientation_parms = self.cop.get_param_vals()
    assert len(vector) == len(cell_parms) + len(orientation_parms)
    tst_cell = vector[:len(cell_parms)]
    tst_orientation = vector[len(cell_parms):len(cell_parms) +
                             len(orientation_parms)]

    self.cucp.set_param_vals(tst_cell)
    self.cop.set_param_vals(tst_orientation)

    from scitbx import matrix

    if do_print:
      print 'Cell: %.3f %.3f %.3f %.3f %.3f %.3f' % \
        tuple(self.crystal.get_unit_cell().parameters())
      print 'Phi(1,2,3): %.3f %.3f %.3f' % tuple(tst_orientation)

    UB = matrix.sqr(self.crystal.get_A())

    score = self.score(UB)

    if score < self.best_score:
      self.best_score = score
      self.cells.append(self.crystal.get_unit_cell().parameters())
    return score

  def score_old(self, RUB):
    score = 0.0
    for j in range(self.data.size()):
      hkl = self.data['miller_index'][j]
      q = RUB * hkl
      score += self.i_s[j] * abs(
        (q + self.s0).length() - matrix.col(self.data['s1'][j]).length())
    return score

  def score(self, UB):
    score = 0.0
    for j in range(self.data.size()):
      hkl = self.data['miller_index'][j]
      q = UB * hkl
      qo = self.qobs[j]
      score += (q - qo).length() ** 2
    return score

  def plotify(self):
    vector = self.best
    cell_parms = self.cucp.get_param_vals()
    orientation_parms = self.cop.get_param_vals()
    assert len(vector) == len(cell_parms) + len(orientation_parms)
    tst_cell = vector[:len(cell_parms)]
    tst_orientation = vector[len(cell_parms):len(cell_parms) +
                             len(orientation_parms)]

    self.cucp.set_param_vals(tst_cell)
    self.cop.set_param_vals(tst_orientation)

    from scitbx import matrix

    UB = matrix.sqr(self.crystal.get_A())
    data = self.data
    self.maxq = 0
    for j in range(data.size()):
      hkl = data['miller_index'][j]
      q = UB * hkl
      qo = self.qobs[j]
      print (q - qo).length(), self.i_s[j], self.dq0[j]
      if (q - qo).length() > self.maxq:
        self.maxq = (q - qo).length()

    return

  def get_maxq(self):
    vector = self.best
    cell_parms = self.cucp.get_param_vals()
    orientation_parms = self.cop.get_param_vals()
    assert len(vector) == len(cell_parms) + len(orientation_parms)
    tst_cell = vector[:len(cell_parms)]
    tst_orientation = vector[len(cell_parms):len(cell_parms) +
                             len(orientation_parms)]

    self.cucp.set_param_vals(tst_cell)
    self.cop.set_param_vals(tst_orientation)

    from scitbx import matrix

    UB = matrix.sqr(self.crystal.get_A())
    data = self.data
    self.maxq = 0
    for j in range(data.size()):
      hkl = data['miller_index'][j]
      q = UB * hkl
      qo = self.qobs[j]
      if (q - qo).length() > self.maxq:
        self.maxq = (q - qo).length()

    return self.maxq

  def zero(self):
    from scitbx import matrix

    UB = matrix.sqr(self.crystal.get_A())
    data = self.data
    self.dq0 = []
    for j in range(data.size()):
      hkl = data['miller_index'][j]
      q = UB * hkl
      qo = self.qobs[j]
      self.dq0.append((q - qo).length())

    return

  def get_signal_mask(self, scale=2):
    if hasattr(self, 'signal_mask'):
      return self.signal_mask
    distance_map = self.render_distance()
    maxq = self.get_maxq()
    self.signal_mask = (distance_map.as_1d() < (2 * maxq))
    self.signal_mask.reshape(self.raw_data.accessor())
    return self.signal_mask

  def make_background(self):
    import copy
    if hasattr(self, 'background'):
      return self.background

    # raw background data
    background = copy.deepcopy(self.raw_data).as_double()

    # mask out the signal areas
    mask = self.get_signal_mask()
    background.as_1d().set_selected(mask.as_1d(), 0.0)
    inv_mask = (~mask).as_1d().as_int()
    inv_mask.reshape(self.raw_data.accessor())

    from dials.algorithms.image.filter import summed_area
    from dials.array_family import flex

    summed_background = summed_area(background, (5,5))
    summed_mask = summed_area(inv_mask, (5,5))
    mean_background = (summed_background /
                       summed_mask.as_double())
    background.as_1d().set_selected(mask.as_1d(), mean_background.as_1d())

    self.background = background
    return background

  def get_background_subtracted_spots(self):
    if hasattr(self, 'background_subtracted_spots'):
      return self.background_subtracted_spots
    mask = self.get_signal_mask()
    background = self.make_background()

    import copy
    background_subtracted_spots = self.raw_data.as_double() - background
    background_subtracted_spots.as_1d().set_selected(~mask.as_1d(), 0)
    self.background_subtracted_spots = background_subtracted_spots
    return background_subtracted_spots

  def integrate(self):
    from scitbx.array_family import flex
    from scitbx import matrix
    nslow, nfast = self.raw_data.focus()

    binary_map = self.get_signal_mask().as_1d().as_int()
    binary_map.reshape(flex.grid(1, nslow, nfast))

    # find connected regions of spots - hacking code for density modification
    # this is used to determine the integration masks for the reflections
    from cctbx import masks
    from cctbx import uctbx
    uc = uctbx.unit_cell((1, nslow, nfast, 90, 90, 90))
    flood_fill = masks.flood_fill(binary_map, uc)
    binary_map = binary_map.as_1d()

    coms = flood_fill.centres_of_mass()

    # now iterate through these blobs, find the intensity and error, and
    # find the Miller index

    UB = matrix.sqr(self.crystal.get_A())
    UBi = UB.inverse()

    winv = 1 / self.beam.get_wavelength()

    data = self.raw_data.as_double()
    background = self.make_background()

    from dials.array_family import flex

    reflections = flex.reflection_table()

    num_pixels_foreground = flex.int()
    background_mean = flex.double()
    background_sum_value = flex.double()
    background_sum_variance = flex.double()
    intensity_sum_value = flex.double()
    intensity_sum_variance = flex.double()
    miller_index = flex.miller_index()
    xyzcal_px = flex.vec3_double()
    bbox = flex.int6()

    fast = flex.int(self.raw_data.size(), -1)
    fast.reshape(self.raw_data.accessor())
    slow = flex.int(self.raw_data.size(), -1)
    slow.reshape(self.raw_data.accessor())

    nslow, nfast = fast.focus()
    for j in range(nslow):
      for i in range(nfast):
        fast[(j, i)] = i
        slow[(j, i)] = j

    for j in range(flood_fill.n_voids()):
      xy = coms[j][2], coms[j][1]
      p = matrix.col(self.panel.get_pixel_lab_coord(xy)).normalize() * winv
      q = p - matrix.col(self.beam.get_s0())
      hkl = UBi * q
      ihkl = [int(round(h)) for h in hkl]
      sel = binary_map == (j + 2)
      pixels = data.select(sel)
      if flex.min(pixels) < 0:
        continue
      n = pixels.size()
      d = flex.sum(pixels)
      b = flex.sum(background.select(sel))
      s = d - b

      # puzzle out the bounding boxes - hack here, we have maps with the
      # fast and slow positions in; select from these then find max, min of
      # this selection
      fs = fast.select(sel)
      f_min = flex.min(fs)
      f_max = flex.max(fs)

      ss = slow.select(sel)
      s_min = flex.min(ss)
      s_max = flex.max(ss)

      bbox.append((f_min, f_max+1, s_min, s_max+1, 0, 1))

      num_pixels_foreground.append(n)
      background_mean.append(b/n)
      background_sum_value.append(b)
      background_sum_variance.append(b)
      intensity_sum_value.append(s)
      intensity_sum_variance.append(d+b)
      miller_index.append(ihkl)
      xyzcal_px.append((xy[0], xy[1], 0.0))

    reflections['num_pixels.foreground'] = num_pixels_foreground
    reflections['background.mean'] = background_mean
    reflections['background.sum.value'] = background_sum_value
    reflections['background.sum.variance'] = background_sum_variance
    reflections['intensity.sum.value'] = intensity_sum_value
    reflections['intensity.sum.variance'] = intensity_sum_variance
    reflections['miller_index'] = miller_index
    reflections['xyzcal.px'] = xyzcal_px
    reflections['id'] = flex.int(miller_index.size(), 0)
    reflections['panel'] = flex.size_t(miller_index.size(), 0)
    reflections['bbox'] = bbox

    from dials.algorithms.shoebox import MaskCode

    reflections["shoebox"] = flex.shoebox(
      reflections["panel"],
      reflections["bbox"],
      allocate=True)

    reflections.extract_shoeboxes(self.imageset, verbose=True)

    # now iterate through these (how?!) and assign mask values
    fg = self.get_signal_mask()
    for reflection in reflections:
      s = reflection['shoebox']
      b = reflection['bbox']
      dz, dy, dx = s.mask.focus()
      for j in range(dy):
        for i in range(dx):
          _i = b[0]
          _j = b[2]
          if fg[(j+_j,i+_i)]:
            m = MaskCode.Valid|MaskCode.Foreground
          else:
            m = MaskCode.Valid|MaskCode.Background
          s.mask[(0,j,i)] = m

    return reflections

if __name__ == '__main__':
  apple = Apple(sys.argv[1], sys.argv[2])
  hklout = sys.argv[3]
  distance_map = apple.render_distance()

  # FIXME at this point subtract background from every pixel - estimate the
  # background from a summed area table - will need to mash around the
  # definitions of foreground and background to do this, and will have to have
  # some idea of the typical size of spots (in order to be able to assign a
  # sensible kernel size)

  mask = apple.get_signal_mask()
  background = apple.make_background()
  spot = apple.get_background_subtracted_spots()

  apple.plot_log_map(spot, 'spot.png')
  apple.plot_log_map(background, 'background.png')
  apple.plot_map(mask, 'mask.png')

  # FIXME at this stage iterate over the image discovering all of the connected
  # components (also known as spots) - integrate them and then determine the
  # Miller index; create an integrated.pickle; shoebox etc.

  reflections = apple.integrate()
  reflections.as_pickle(hklout)
