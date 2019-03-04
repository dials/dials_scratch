
from __future__ import division
from __future__ import print_function


def predict(experiments, bandpass=0.035):
  '''
  Predict all reflections within the limiting Ewald spheres given a fractional
  bandpass (i.e. smin = s0 / (1 + bandpass), smax = s0 / (1 - bandpass))

  '''

  from dials.algorithms.spot_prediction import IndexGenerator
  from scitbx import matrix
  from dials.array_family import flex

  print("Predicting reflections")

  assert len(experiments) == 1
  assert len(experiments[0].detector) == 1

  # Get the resolution at the corners of the detector
  d_min = experiments[0].detector.get_max_resolution(
    experiments[0].beam.get_s0())

  # Get the unit cell and space group
  unit_cell = experiments[0].crystal.get_unit_cell()
  space_group = experiments[0].crystal.get_space_group().type()
  A = matrix.sqr(experiments[0].crystal.get_A())

  # Generate indices
  index_generator = IndexGenerator(unit_cell, space_group, d_min)

  # Get the possible miller indices
  indices = index_generator.to_array()

  # Setup the limiting spheres
  s0 = matrix.col(experiments[0].beam.get_s0())
  smin = s0 / (1 + bandpass)
  smax = s0 / (1 - bandpass)

  # Check if the index is within the limiting spheres
  miller_index = flex.miller_index()
  rlp = flex.vec3_double()
  s1 = flex.vec3_double()
  xyzcalpx = flex.vec3_double()
  xyzcalmm = flex.vec3_double()
  panel = experiments[0].detector[0]
  xsize, ysize = panel.get_image_size()
  print(len(indices))
  for h in indices:
    p = A * matrix.col(h)

    d0 = (p + smin).length()
    d1 = (p + smax).length()
    R = 0.00
    i0 = max(d0-R, smin.length())
    i1 = min(d1+R, smax.length())
    if i0 <= i1:
    #if (p + smin).length() >= smin.length() and (p + smax).length() <= smax.length():
      #t = -p.length_sq() / (2*p.dot(s0))
      #s =  t*s0 + p
      s = s0 + p

      xypx = panel.get_ray_intersection_px(s)
      xymm = panel.pixel_to_millimeter(xypx)

      if xypx[0] >= 0 and xypx[0] < xsize and xypx[1] >= 0 and xypx[1] < ysize:
        miller_index.append(h)
        rlp.append(p)
        s1.append(s)
        xyzcalpx.append((xypx[0], xypx[1], 0))
        xyzcalmm.append((xymm[0], xymm[1], 0))


  # Construct the reflection table
  reflections = flex.reflection_table()
  reflections['id'] = flex.size_t(len(miller_index), 0)
  reflections['panel'] = flex.size_t(len(miller_index), 0)
  reflections['miller_index'] = miller_index
  reflections['rlp'] = rlp
  reflections['xyzcal.px'] = xyzcalpx
  reflections['xyzcal.mm'] = xyzcalmm

  from matplotlib import pylab
  X, Y, _ = reflections['xyzcal.px'].parts()
  pylab.scatter(X,Y)
  pylab.show()


  print("Predicted %d reflections" % len(reflections))

  # Return the reflections
  return reflections


def label_pixels(experiments, mask_distance=0.3):
  '''
  Assign each pixel to a miller index

  '''
  from dials.algorithms.spot_prediction import PixelToMillerIndex
  from math import sqrt, floor
  from dials.array_family import flex
  from collections import defaultdict
  from scitbx import matrix

  experiment = experiments[0]
  detector = experiment.detector
  assert len(detector) == 1
  panel = detector[0]

  xsize, ysize = panel.get_image_size()

  # A class to map pixels to miller indices
  transform = PixelToMillerIndex(
    experiment.beam,
    experiment.detector,
    experiment.crystal)

  s0 = matrix.col(experiments[0].beam.get_s0())
  smin = s0 * (1 - 0.035)
  smax = s0 * (1 + 0.035)
  A = matrix.sqr(experiment.crystal.get_A())

  mask = flex.bool(flex.grid(ysize, xsize), False)
  reflections = defaultdict(list)

  # For each pixel, assign to a miller index and also set whether it is
  # foreground or background
  for j in range(ysize):
    for i in range(xsize):
      h = transform.h(0, i, j)
      h0 = tuple(map(lambda x: int(floor(x+0.5)), h))
      d = sqrt(sum(map(lambda x,y: (x-y)**2, h, h0)))

      if d < mask_distance:
        mask[j,i] = True
      if d < mask_distance * 2:
        reflections[h0].append((j,i))

  return reflections, mask

def extract_shoeboxes(experiments, reflections, labels, mask):
  from dials.array_family import flex
  from dials.model.data import Shoebox
  from dials.algorithms.shoebox import MaskCode

  # Get the image data
  data = experiments[0].imageset.get_raw_data(0)[0]

  print("Extracting shoeboxes")
  bbox = flex.int6()
  shoebox = flex.shoebox()
  selection = flex.size_t()
  for index, h in enumerate(reflections['miller_index']):
    if h not in labels:
      continue

    #assert h in labels
    pixels = labels[h]

    b_sum = 0
    f_sum = 0
    b_cnt = 0
    f_cnt = 0
    Y, X = zip(*pixels)
    x0, x1, y0, y1 = min(X), max(X)+1, min(Y), max(Y)+1

    sbox = Shoebox(0, (x0, x1, y0, y1, 0, 1))
    sbox.allocate(0)
    n_foreground = 0
    for y, x in zip(Y, X):
      #print x, y, reflections[index]['xyzcal.px']
      j = y - y0
      i = x - x0
      sbox.data[0,j,i] = data[y,x]
      if data[y,x] >= 0:
        sbox.mask[0,j,i] = MaskCode.Valid
      if mask[y,x]:
        sbox.mask[0,j,i] |= MaskCode.Foreground
        n_foreground += 1
      else:
        sbox.mask[0,j,i] |= MaskCode.Background

    if n_foreground == 0:
      print("No foreground in ", h)
      continue
    selection.append(index)
    bbox.append((x0, x1, y0, y1, 0, 1))
    shoebox.append(sbox)

  print("Selecting %d reflections in labels" % len(selection))

  reflections = reflections.select(selection)
  print(len(reflections))

  reflections['bbox'] = bbox
  reflections['shoebox'] = shoebox

  return reflections

def integrate(experiments, reflections, mask_distance=0.3):
  '''
  Do a crude integration around all the potential spots

  '''

  print("Labelling pixels")
  labels, mask = label_pixels(experiments, mask_distance=mask_distance)

  reflections = extract_shoeboxes(experiments, reflections, labels, mask)

  reflections.compute_background(experiments)
  # reflections.compute_centroid(experiments)
  reflections.compute_summed_intensity()

  return reflections


def compute_wavelength(experiments, reflections):
  from scitbx import matrix
  from dials.array_family import flex
  wavelength = flex.double()
  s0 = matrix.col(experiments[0].beam.get_s0())
  for r in reflections:
    p = matrix.col(r['rlp'])
    t = -p.length_sq() / (2*p.dot(s0))
    wavelength.append(1/(s0.length() * t))
  reflections['wavelength'] = wavelength
  return reflections


def refine(experiments, bandpass=0.035, isigi_threshold=5, mask_distance=0.3):
  from dials.array_family import flex
  from math import floor

  # Overpredict reflections with inflated bandpass
  reflections = predict(experiments, bandpass=bandpass)

  # Integrate the reflections
  reflections = integrate(experiments, reflections, mask_distance=mask_distance)

  # Select integrated
  selection = reflections.get_flags(reflections.flags.integrated_sum)
  subset = reflections.select(selection)
  print("Integrated %d reflections" % len(subset))

  # Select variance > 0
  variance = subset['intensity.sum.variance']
  selection = variance > 0
  subset = subset.select(selection)
  print("Selecting %d reflections with variance > 0" % len(subset))

  # Select I/Sig(I) > 5
  intensity = subset['intensity.sum.value']
  variance = subset['intensity.sum.variance']
  i_over_s = intensity / flex.sqrt(variance)
  # from matplotlib import pylab
  # pylab.hist(i_over_s, bins=10)
  # pylab.show()
  selection = i_over_s > isigi_threshold
  subset = subset.select(selection)
  print("Selecting %d strong reflections" % len(subset))


  subset = compute_wavelength(experiments, subset)

  l = experiments[0].beam.get_wavelength()
  min_wavelength = l * (1 - bandpass)
  max_wavelength = l * (1 + bandpass)
  nbins = 10
  wl_step = (max_wavelength - min_wavelength) / (nbins-1)

  hist = [0] * nbins
  for r in subset:
    wl = r['wavelength']
    index = int(floor((wl - min_wavelength) / wl_step))
    assert index >= 0 and index < len(hist)
    hist[index] += 1

  from matplotlib import pylab
  pylab.hist(subset['wavelength'])
  pylab.show()

  from matplotlib import pylab
  X, Y, _ = reflections['xyzcal.px'].parts()
  pylab.scatter(X,Y)
  pylab.show()

  print("Saving reflections to refined.pickle")
  reflections.as_pickle("refined.pickle")


if __name__ == '__main__':

  from dxtbx.model.experiment_list import ExperimentListFactory
  import sys

  experiments_filename = sys.argv[1]

  bandpass = 0.05
  isigi_threshold = 5
  mask_distance=0.15

  # Read the experiments
  experiments = ExperimentListFactory.from_json_file(experiments_filename)

  # Do the refinement
  refine(experiments,
         bandpass        = bandpass,
         isigi_threshold = isigi_threshold,
         mask_distance   = mask_distance)
