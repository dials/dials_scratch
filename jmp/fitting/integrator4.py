from __future__ import division
import dials

def read_experiments(filename):
  from dxtbx.model.experiment_list import ExperimentListFactory
  return ExperimentListFactory.from_json_file(filename)

def read_reference(filename):
  import cPickle as pickle
  return pickle.load(open(filename))

def read_reflections(filename):
  from dials.array_family import flex
  reflections = []
  for f in filename:
    print "Reading %s" % f
    r = flex.reflection_table.from_pickle(f)
    reflections.append(r)
    del r['shoebox']
  return reflections


class MaskCalculatorFactory(object):

  @classmethod
  def build(Class, experiments):
    from dials_scratch.jmp.fitting import MaskCalculator
    return MaskCalculator(
      experiments[0].beam,
      experiments[0].detector,
      experiments[0].goniometer,
      experiments[0].scan,
      experiments[0].profile.delta_b(deg=False),
      experiments[0].profile.delta_m(deg=False))


class BackgroundCalculatorFactory(object):

  @classmethod
  def build(Class, experiments):
    from dials_scratch.jmp.fitting import BackgroundCalculator
    from dials.algorithms.background.glm import Creator
    return BackgroundCalculator(
      Creator.model.constant3d,
      1.345,
      100,
      10)

class IntensityCalculatorFactory(object):

  @classmethod
  def build(self, experiments, reference, grid_size=5, detector_space=False,
            deconvolution=False):
    from dials_scratch.jmp.fitting import IntensityCalculator
    from dials_scratch.jmp.fitting import Reference
    from dials.algorithms.profile_model.modeller import CircleSampler
    from dials.array_family import flex
    from dials.algorithms.profile_model.gaussian_rs.transform import TransformSpec
    from dials.algorithms.profile_model.gaussian_rs import CoordinateSystem

    assert len(reference) % 9 == 0
    num_scan_points = len(reference) // 9

    sampler = CircleSampler(
      experiments[0].detector[0].get_image_size(),
      experiments[0].scan.get_array_range(),
      num_scan_points)


    spec = TransformSpec(
      experiments[0].beam,
      experiments[0].detector,
      experiments[0].goniometer,
      experiments[0].scan,
      experiments[0].profile.sigma_b(deg=False),
      experiments[0].profile.sigma_m(deg=False),
      experiments[0].profile.n_sigma() * 1.5,
      grid_size)

    temp = reference

    reference = Reference()
    for d, m in temp:
      reference.append(d, m)
    print detector_space, deconvolution
    return IntensityCalculator(
      reference,
      sampler,
      spec,
      detector_space,
      deconvolution)


# class IntensityCalculator(object):

#   def __init__(self, experiments, reference, grid_size=5, detector_space=False):
#     self.experiments = experiments
#     self.reference = reference
#     self.detector_space = detector_space
#     self.grid_size = grid_size

#   def __call__(self, reflection):
#     from dials.algorithms.profile_model.modeller import CircleSampler
#     from dials.array_family import flex
#     from dials.algorithms.profile_model.gaussian_rs.transform import TransformReverse
#     from dials.algorithms.profile_model.gaussian_rs.transform import TransformForward
#     from dials.algorithms.profile_model.gaussian_rs.transform import TransformReverseNoModel
#     from dials.algorithms.profile_model.gaussian_rs.transform import TransformSpec
#     from dials.algorithms.profile_model.gaussian_rs import CoordinateSystem
#     from dials.algorithms.integration.fit import ProfileFitter
#     from dials.array_family import flex
#     from dials.model.data import make_image
#     from dials.array_family import flex
#     from dials_scratch.jmp.fitting import get_reflection
#     reflections = reflection.to_table()

#     grid_size = self.grid_size
#     detector_space = self.detector_space

#     assert len(self.reference) % 9 == 0
#     num_scan_points = len(self.reference) // 9

#     sampler = CircleSampler(
#       self.experiments[0].detector[0].get_image_size(),
#       self.experiments[0].scan.get_array_range(),
#       num_scan_points)


#     spec = TransformSpec(
#       self.experiments[0].beam,
#       self.experiments[0].detector,
#       self.experiments[0].goniometer,
#       self.experiments[0].scan,
#       self.experiments[0].profile.sigma_b(deg=False),
#       self.experiments[0].profile.sigma_m(deg=False),
#       self.experiments[0].profile.n_sigma() * 1.5,
#       grid_size)

#     m2 = self.experiments[0].goniometer.get_rotation_axis()
#     s0 = self.experiments[0].beam.get_s0()

#     Iprf = flex.double(len(reflections))
#     Vprf = flex.double(len(reflections))
#     Cprf = flex.double(len(reflections))
#     Fprf = flex.bool(len(reflections))
#     Part = reflections['partiality']

#     reflections['intensity.prf_old.value'] = reflections['intensity.prf.value']
#     reflections['intensity.prf_old.variance'] = reflections['intensity.prf.variance']

#     selection = reflections.get_flags(reflections.flags.integrated_prf)

#     reflections.unset_flags(~Fprf, reflections.flags.integrated_prf)

#     for i, r in enumerate(reflections):

#       if selection[i] == False:
#         continue

#       s1 = r['s1']
#       phi = r['xyzcal.mm'][2]
#       xyz = r['xyzcal.px']
#       bbox = r['bbox']
#       panel = r['panel']
#       image = r['shoebox'].data.as_double()
#       background = r['shoebox'].background.as_double()
#       mask = (r['shoebox'].mask.as_1d() == 5)#| (r['shoebox'].mask.as_1d() == 3)
#       mask.reshape(image.accessor())
#       cs = CoordinateSystem(m2, s0, s1, phi)

#       index = sampler.nearest(0, xyz)

#       profile, profile_mask = self.reference[index]

#       # print flex.sum(profile)
#       # print r['partiality']

#       if detector_space:

#         transform = TransformReverseNoModel(
#           spec,
#           cs,
#           bbox,
#           panel,
#           profile)
#         p = transform.profile()
#         d = image
#         m = mask
#         b = background
#         Part[i] = flex.sum(p)
#         try:
#           fit = ProfileFitter(d, b, m, p, 1e-3, 100)
#           assert fit.niter() < 100
#           Iprf[i] = fit.intensity()
#           Vprf[i] = fit.variance()
#           Cprf[i] = fit.correlation()
#           Fprf[i] = True
#           print r['miller_index'], i, fit.intensity(),  r['intensity.sum.value'],  r['intensity.prf_old.value'], Part[i], fit.niter()
#         except Exception, e:
#           print e
#           pass

#       else:

#         try:

#           transform = TransformForward(
#             spec,
#             cs,
#             bbox,
#             panel,
#             image,
#             background,
#             mask)

#           p = profile
#           d = transform.profile()
#           b = transform.background()
#           m = transform.mask() & profile_mask

#           fit = ProfileFitter(d, b, m, p, 1e-3, 100)
#           assert fit.niter() < 100
#           Iprf[i] = fit.intensity()[0]
#           Vprf[i] = fit.variance()[0]
#           Cprf[i] = fit.correlation()
#           Fprf[i] = True
#           print r['miller_index'], i, fit.intensity(), r['intensity.prf_old.value']
#         except Exception, e:
#           pass

#     reflections['intensity.prf.value'] = Iprf
#     reflections['intensity.prf.variance'] = Vprf
#     reflections['intensity.prf.correlation'] = Cprf
#     reflections.set_flags(Fprf, reflections.flags.integrated_prf)
#     return get_reflection(reflections, 0)



# class Integrator(object):

#   def __init__(self,
#                reflections,
#                imageset,
#                compute_mask,
#                compute_background,
#                compute_intensity):
#     from dials.array_family import flex
#     from dials.model.data import make_image

#     frame0, frame1 = imageset.get_scan().get_array_range()

#     reflections['shoebox'] = flex.shoebox(reflections['panel'],
#                                           reflections['bbox'], allocate=True)
#     extractor = flex.ShoeboxExtractor(reflections, 1, frame0, frame1)

#     for i in range(len(imageset)):
#       print "Reading image %d" % i
#       data = imageset.get_raw_data(i)
#       mask = imageset.get_mask(i)
#       extractor.next(make_image(data, mask))

#     compute_mask(reflections)

#     compute_background(reflections)

#     print "Computing centroid"
#     reflections.compute_centroid(experiments)

#     print "Computing summed intensity"
#     reflections.compute_summed_intensity()

#     compute_intensity(reflections)

#     del reflections['shoebox']

#     self._reflections = reflections

#   def reflections(self):
#     return self._reflections


def integrate_job(block, experiments, reflections, reference, grid_size=5,
                  detector_space=False, deconvolution=False):
  from dials_scratch.jmp.fitting import Integrator
  from dials.array_family import flex

  reflections["intensity.prf_old.value"] = reflections["intensity.prf.value"]
  reflections["intensity.prf_old.variance"] = reflections["intensity.prf.variance"]
  reflections["intensity.prf.value"] = flex.double(len(reflections))
  reflections["intensity.prf.variance"] = flex.double(len(reflections))

  integrator = Integrator(
    reflections        = reflections,
    imageset           = experiments[0].imageset[block[0]:block[1]],
    compute_mask       = MaskCalculatorFactory.build(experiments),
    compute_background = BackgroundCalculatorFactory.build(experiments),
    compute_intensity  = IntensityCalculatorFactory.build(experiments, reference,
                                             grid_size=grid_size,
                                             detector_space=detector_space,
                                             deconvolution=deconvolution),
    nthreads           = 1,
    debug = False
  )

  reflections = integrator.reflections()

  # from dials.algorithms.shoebox import MaskCode
  # from matplotlib import pylab
  # for sbox in reflections["shoebox"]:
  #   if sbox.count_mask_values(MaskCode.Overlapped) > 0:
  #     for i in range(sbox.mask.all()[0]):
  #       print "slice", i
  #       pylab.imshow(sbox.mask.as_numpy_array()[i,:,:], interpolation='none')
  #       pylab.show()

  # del reflections["shoebox"]

  # exit(0)

  return reflections


def integrate(experiments, reflections, reference, grid_size=5,
              detector_space=False, deconvolution=False):
  from dials.array_family import flex
  from dials.array_family import flex
  from dials.model.data import make_image

  block_size = 20
  blocks = []
  start = 0
  while start < len(experiments[0].imageset) - block_size // 2:
    blocks.append((start, start+block_size))
    start += block_size // 2
  print blocks
  assert len(blocks) == len(reflections)

  result = flex.reflection_table()
  for i in range(len(blocks)):

    reflections[i].compute_bbox(experiments)

    for j in range(len(reflections[i])):
      x0, x1, y0, y1, z0, z1 = reflections[i]["bbox"][j]
      if z0 < blocks[i][0]:
        z0 = blocks[i][0]
      if z1 > blocks[i][1]:
        z1 = blocks[i][1]
      reflections[i]["bbox"][j] = (x0, x1, y0, y1, z0, z1)

    r = integrate_job(blocks[i], experiments, reflections[i], reference, grid_size,
                  detector_space, deconvolution)
    result.extend(r)

  return result


if __name__ == '__main__':
  from os.path import join
  import sys

  experiments_filename = sys.argv[1]
  reference_filename = sys.argv[2]
  grid_size = int(sys.argv[3])
  if sys.argv[4] == "True":
    detector_space = True
  elif sys.argv[4] == "False":
    detector_space = False
  else:
    raise RuntimeError(sys.argv[4])
  reflections_filename = sys.argv[5:]

  experiments = read_experiments(experiments_filename)
  reflections = read_reflections(reflections_filename)
  reference = read_reference(reference_filename)

  experiments[0].profile._sigma_b *= 2

  print "Read %d reflections" % len(reflections)

  deconvolution = True

  from time import time
  st = time()
  reflections = integrate(experiments, reflections, reference[0],
                          grid_size=grid_size, detector_space=detector_space,
                          deconvolution=deconvolution)
  print "Num profile fitted", reflections.get_flags(reflections.flags.integrated_prf).count(True)
  print "Time taken: ", time() - st

  reflections.as_pickle("integrated.pickle")
