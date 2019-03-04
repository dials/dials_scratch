from __future__ import division
from __future__ import print_function
import dials

def read_experiments(filename):
  from dxtbx.model.experiment_list import ExperimentListFactory
  print("Reading %s" % filename)
  return ExperimentListFactory.from_json_file(filename)

def read_reference(filename):
  import cPickle as pickle
  print("Reading %s" % filename)
  return pickle.load(open(filename))

def predict_reflections(experiments):
  from dials.array_family import flex

  print("Predicting reflections")

  # Predict
  reflections = flex.reflection_table.from_predictions_multi(experiments)

  # Compute some reflection properties
  reflections.compute_zeta_multi(experiments)
  reflections.compute_d(experiments)
  reflections.compute_bbox(experiments)

  # Filter the reflections by zeta
  mask = flex.abs(reflections['zeta']) < 0.05
  num_ignore = mask.count(True)
  reflections.set_flags(mask, reflections.flags.dont_integrate)

  print("Predicted %d reflections" % len(reflections))
  print("Ignoring %d reflections" % num_ignore)

  return reflections


def read_images(experiments):
  from dials.array_family import flex

  xsize, ysize = experiments[0].detector[0].get_image_size()
  zsize = experiments[0].scan.get_num_images()

  data = flex.double(flex.grid(zsize, ysize, xsize))
  mask = flex.int(flex.grid(zsize, ysize, xsize))

  iset = experiments[0].imageset

  for i in range(len(iset)):
    print("Reading image %d" % i)
    d = iset.get_raw_data(i)[0].as_double()
    m = iset.get_mask(i)[0].as_1d().as_int()
    m.reshape(flex.grid(ysize, xsize))
    d.reshape(flex.grid(1, ysize, xsize))
    m.reshape(flex.grid(1, ysize, xsize))
    data[i:i+1,:,:] = d
    mask[i:i+1,:,:] = m

  return data, mask


# def get_pixel_list(experiments, reflections):
#   from collections import defaultdict
#   from dials.array_family import flex
#   print "Get reflection pixel list"

#   xsize, ysize = experiments[0].detector[0].get_image_size()
#   zsize = experiments[0].scan.get_num_images()

#   pixels = defaultdict(list)
#   for i in range(len(reflections)):

#     subset = reflections.select(flex.size_t([i]))
#     subset['shoebox'] = flex.shoebox(
#       subset['panel'],
#       subset['bbox'],
#       allocate=True)
#     subset.compute_mask(experiments)

#     x0, x1, y0, y1, z0, z1 = subset['bbox'][0]
#     mask = subset['shoebox'][0].mask
#     for z in range(mask.all()[0]):
#       for y in range(mask.all()[1]):
#         for x in range(mask.all()[2]):
#           if (z0 + z >= 0 and z0 + z < zsize and
#               y0 + y >= 0 and y0 + y < ysize and
#               x0 + x >= 0 and x0 + x < xsize):
#             if mask[z,y,x] == 4:
#               pixels[(x0 + x, y0 + y, z0 + z)].append(i)

#   print "Reflections recorded on %d pixels" % len(pixels)
#   return pixels

# def find_overlaps(pixels):
#   from collections import defaultdict
#   print "Finding overlaps"
#   overlaps = defaultdict(list)

#   for k, v in pixels.iteritems():
#     if len(v) > 1:
#       for j in range(len(v)-1):
#         for i in range(j+1, len(v)):
#           overlaps[v[j]].append(v[i])
#           overlaps[v[i]].append(v[j])
#   overlaps = dict((k, set(v)) for k, v in overlaps.iteritems())
#   for k, v in overlaps.iteritems():
#     if len(v) > 0:
#       print "Reflection %d operlaps with %s" % (k, list(v))

#   return overlaps


# def update_mask(mask, pixels):
#   print "Updating the mask with reflection info"
#   for k, v in pixels.iteritems():
#     x, y, z = k
#     mask[z, y, x] = mask[z, y, x] | 4

#   # from matplotlib import pylab
#   # pylab.imshow(mask.as_numpy_array()[0,:,:], interpolation='none')
#   # pylab.show()
#   return mask


# def compute_background(data, mask, reflections):
#   from dials.algorithms.background.glm import RobustPoissonMean
#   from dials.array_family import flex
#   print "Computing the background"
#   zsize, ysize, xsize = data.all()

#   background = flex.double()
#   success = flex.bool()

#   for i in range(len(reflections)):
#     x0, x1, y0, y1, z0, z1 = reflections[i]['bbox']
#     x0 = max(x0, 0)
#     y0 = max(y0, 0)
#     z0 = max(z0, 0)
#     x1 = min(x1, xsize)
#     y1 = min(y1, ysize)
#     z1 = min(z1, zsize)
#     pixels = flex.double()
#     for z in range(z0, z1):
#       for y in range(y0, y1):
#         for x in range(x0, x1):
#           if mask[z,y,x] == 1:
#             pixels.append(data[z,y,x])

#     if pixels.all_eq(0):
#       background.append(0)
#       success.append(0)
#     else:
#       try:
#         median = sorted(pixels)[len(pixels)//2]
#         median = max(1, median)
#         b = RobustPoissonMean(pixels, median)
#         background.append(b.mean())
#         success.append(True)
#       except Exception:
#         success.append(False)

#   return background, success


def load_sampler(experiments, reference):
  from dials.algorithms.profile_model.modeller import CircleSampler
  assert len(reference[0]) % 9 == 0
  num_scan_points = len(reference[0]) // 9

  sampler = CircleSampler(
    experiments[0].detector[0].get_image_size(),
    experiments[0].scan.get_array_range(),
    num_scan_points)


  return sampler

# def load_transform(experiments, reference):
#   from dials.algorithms.profile_model.gaussian_rs.transform import TransformSpec

#   n_sigma = 4.0
#   grid_size = 25
#   assert reference[0][0].all() == (2*grid_size+1, 2*grid_size+1, 2*grid_size+1)
#   spec = TransformSpec(
#     experiments[0].beam,
#     experiments[0].detector,
#     experiments[0].goniometer,
#     experiments[0].scan,
#     experiments[0].profile.sigma_b(deg=False),
#     experiments[0].profile.sigma_m(deg=False),
#     n_sigma,
#     grid_size)

#   return spec

# def get_reverse_pixel_list(pixels):
#   from collections import defaultdict
#   reverse_pixels = defaultdict(list)

#   for k, v in pixels.iteritems():
#     for item in v:
#       reverse_pixels[item].append(k)

#   return reverse_pixels


# def get_profile_on_detector(experiments, transform_spec, reflection, sampler, reference):
#   from dials.algorithms.profile_model.gaussian_rs.transform import TransformReverse
#   from dials.algorithms.profile_model.gaussian_rs.transform import TransformForward
#   from dials.algorithms.profile_model.gaussian_rs.transform import TransformReverseNoModel
#   from dials.algorithms.profile_model.gaussian_rs import CoordinateSystem

#   m2 = experiments[0].goniometer.get_rotation_axis()
#   s0 = experiments[0].beam.get_s0()

#   s1 = reflection['s1']
#   phi = reflection['xyzcal.mm'][2]
#   xyz = reflection['xyzcal.px']
#   bbox = reflection['bbox']
#   panel = reflection['panel']

#   index = sampler.nearest(0, xyz)

#   profile = reference[0][index]
#   cs = CoordinateSystem(m2, s0, s1, phi)

#   transform = TransformReverseNoModel(
#     transform_spec,
#     cs,
#     bbox,
#     panel,
#     profile)

#   return transform.profile()


# def get_data_on_detector(data, mask, reflections, i, pixel_list,
#                          reverse_pixel_list):

#   from dials.array_family import flex

#   x0, x1, y0, y1, z0, z1 = reflections[i]['bbox']

#   if (x0 < 0 or x1 > data.all()[2] or
#       y0 < 0 or y1 > data.all()[1] or
#       z0 < 0 or z1 > data.all()[0]):
#     return None, None, None


#   r_data = data[z0:z1,y0:y1,x0:x1]
#   r_bgrd = flex.double(r_data.accessor(), reflections[i]['background'])
#   r_mask = mask[z0:z1,y0:y1,x0:x1]

#   r_mask = r_mask.as_1d()
#   r_mask.set_selected(r_mask == 5, 4)
#   r_mask.reshape(r_data.accessor())

#   for p in reverse_pixel_list[i]:
#     if len(pixel_list[p]) == 1:
#       x, y, z = p
#       zz = z-z0
#       yy = y-y0
#       xx = x-x0
#       # print i, p, pixel_list[p]
#       # print xx, yy, zz, reflections[i]['bbox'], r_mask.all()
#       r_mask[zz,yy,xx] = 5


#   return r_data, r_bgrd, r_mask



def integrate(experiments, reflections, reference):
  from dials.array_family import flex
  from dials_scratch.jmp.fitting import PixelList

  # Load the reference sampler
  sampler = load_sampler(experiments, reference)

  # Read all the image data
  data, mask = read_images(experiments)

  print("Computing pixel lookup")
  from time import time
  st = time()
  pixel_list = PixelList(
    reflections,
    experiments[0].beam,
    experiments[0].detector,
    experiments[0].goniometer,
    experiments[0].scan,
    experiments[0].profile.sigma_b(deg=False),
    experiments[0].profile.sigma_m(deg=False))
  print("Time taken: ", time() - st)

  print("Updating mask")
  st = time()
  pixel_list.update_mask(mask)
  print("Time taken: ", time() - st)

  print("Computing background")
  st = time()
  success = pixel_list.compute_background(data, mask, reflections)
  print("Time taken: ", time() - st)

  background = reflections['background']
  print("Mean background %f" % flex.mean(background))
  print("%d failed during background modelling" % (success.count(False)))

  print("Computing intensity")
  st = time()
  success = pixel_list.compute_intensity(data, mask, reflections, reference[0], sampler)
  print("Time taken: ", time() - st)
  print("%d reflection failed to integrate" % success.count(False))

  return reflections
  exit(0)


  # # Get a list of foreground pixels with associated reflections
  # pixel_list = get_pixel_list(experiments, reflections)

  # # Get pixel list associating pixels with reflections
  # reverse_pixel_list = get_reverse_pixel_list(pixel_list)

  # # Find overlaps
  # overlaps = find_overlaps(pixel_list)

  # Update mask
  # mask = update_mask(mask, pixel_list)

  # Compute the background for each reflection
  # background, success = compute_background(data, mask, reflections)
  # reflections['background'] = background
  # num_ignore = success.count(False)
  # reflections.set_flags(~success, reflections.flags.dont_integrate)
  # Isum = flex.double(len(reflections))
  # Vsum = flex.double(len(reflections))
  # Iprf = flex.double(len(reflections))
  # Vprf = flex.double(len(reflections))
  # partiality = flex.double(len(reflections))
  # success = flex.bool(len(reflections), False)
  # selection_to_integrate = ~reflections.get_flags(reflections.flags.dont_integrate)


  # print "Computing intensity"
  # for i in range(len(reflections)):

  #   if selection_to_integrate[i] == False:
  #     continue

  #   # Get the detector profile
  #   profile = get_profile_on_detector(experiments, transform_spec,
  #                                     reflections[i], sampler, reference)

  #   r_data, r_bgrd, r_mask = get_data_on_detector(
  #     data,
  #     mask,
  #     reflections,
  #     i,
  #     pixel_list,
  #     reverse_pixel_list)

  #   if r_data is None:
  #     continue

  #   partiality[i] = flex.sum(profile.as_1d().select(r_mask.as_1d() == 5))

  #   from dials.algorithms.integration.sum import integrate_by_summation
  #   summation = integrate_by_summation(r_data, r_bgrd, r_mask)
  #   Isum[i] = summation.intensity()
  #   Vsum[i] = summation.variance()

  #   try:
  #     from dials.algorithms.integration.fit import ProfileFittingDouble as ProfileFitting
  #     fitting = ProfileFitting(profile, r_mask == 5, r_data, r_bgrd)
  #     assert fitting.niter() < 100
  #     Iprf[i] = fitting.intensity()
  #     Vprf[i] = fitting.variance()
  #     success[i] = True
  #   except Exception, e:
  #     print e
  #     Iprf[i] = 0
  #     Vprf[i] = -1
  #     success[i] = False
  #     pass

  #   print i, Isum[i], Iprf[i], partiality[i]

  # print "%d reflection failed to integrate" % success.count(False)

  # reflections['partiality'] = partiality
  # reflections['intensity.prf.value'] = Iprf
  # reflections['intensity.prf.variance'] = Vprf
  # reflections['intensity.sum.value'] = Isum
  # reflections['intensity.sum.variance'] = Vsum
  # reflections.set_flags(success, reflections.flags.integrated_prf)
  # reflections.set_flags(success, reflections.flags.integrated_sum)

  # return reflections

  # from matplotlib import pylab
  # pylab.imshow(mask.as_numpy_array()[0,:,:], interpolation='none')
  # pylab.show()

if __name__ == '__main__':
  from os.path import join
  import sys

  experiments_filename = sys.argv[1]
  reference_filename = sys.argv[2]

  experiments = read_experiments(experiments_filename)
  reference = read_reference(reference_filename)

  experiments[0].imageset = experiments[0].imageset[0:1]
  experiments[0].scan = experiments[0].imageset.get_scan()

  reflections = predict_reflections(experiments)

  reflections = integrate(experiments, reflections, reference)

  reflections.as_pickle("integrated.pickle")
