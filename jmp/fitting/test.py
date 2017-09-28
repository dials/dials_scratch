from __future__ import division
import dials

def read_experiments(filename):
  from dxtbx.model.experiment_list import ExperimentListFactory
  return ExperimentListFactory.from_json_file(filename)

def read_reflections(filename):
  from dials.array_family import flex
  r = flex.reflection_table()
  for f in filename:
    print "Reading %s" % f
    r.extend(flex.reflection_table.from_pickle(f))
  return r

def select_strong(reflections):
  from dials.array_family import flex
  selection1 = reflections.get_flags(reflections.flags.indexed)
  selection2 = reflections.get_flags(reflections.flags.centroid_outlier)
  selection3 = flex.abs(reflections['zeta']) > 0.05
  selection4 = reflections['partiality'] > 0.99
  selection = (selection1) & (~selection2) & (selection3) & (selection4)
  return reflections.select(selection)

def compute_reference(experiments, reflections):
  from dials.algorithms.profile_model.modeller import CircleSampler
  from dials.array_family import flex
  from dials.algorithms.profile_model.gaussian_rs.transform import TransformForward
  from dials.algorithms.profile_model.gaussian_rs.transform import TransformSpec
  from dials.algorithms.profile_model.gaussian_rs import CoordinateSystem

  reflections = select_strong(reflections)
  print "Selected %d strong spots" % len(reflections)

  sampler = CircleSampler(
    experiments[0].detector[0].get_image_size(),
    experiments[0].scan.get_array_range(),
    1)


  n_sigma = 4.0
  grid_size = 25
  spec = TransformSpec(
    experiments[0].beam,
    experiments[0].detector,
    experiments[0].goniometer,
    experiments[0].scan,
    experiments[0].profile.sigma_b(deg=False),
    experiments[0].profile.sigma_m(deg=False),
    n_sigma,
    grid_size)

  m2 = experiments[0].goniometer.get_rotation_axis()
  s0 = experiments[0].beam.get_s0()

  reference = [flex.double(flex.grid(1+2*grid_size, 1+2*grid_size,
                                     1+2*grid_size)) for i in
               range(len(sampler))]
  count = [0] * len(sampler)

  for r in reflections:
    s1 = r['s1']
    phi = r['xyzcal.mm'][2]
    xyz = r['xyzcal.px']
    bbox = r['bbox']
    panel = r['panel']
    image = r['shoebox'].data.as_double()
    mask = r['shoebox'].mask.as_1d() == 5
    mask.reshape(image.accessor())
    cs = CoordinateSystem(m2, s0, s1, phi)

    try:
      transform = TransformForward(
        spec,
        cs,
        bbox,
        panel,
        image,
        mask)
      d = transform.profile()

      d /= flex.sum(d)

      index = sampler.nearest(0, xyz)
      indices = sampler.nearest_n(0, xyz)
      for i in indices:
        w = sampler.weight(i, 0, xyz)
        reference[i] += w*d
        count[i] += 1
    except Exception:
      pass

  for i in range(len(reference)):
    r = reference[i]
    if flex.sum(r) > 0:
      print flex.max(r)
      g = r.accessor()
      r = r.as_1d()
      s = r > 0.02 * flex.max(r)
      r.set_selected(~s, flex.double(len(r), 0))
      r = r / flex.sum(r)
      r.reshape(g)
      reference[i] = r

  for i in range(len(reference)):
    from matplotlib import pylab
    print count[i]
    r = reference[i]
    d = r.as_numpy_array()[11,:,:]
    pylab.imshow(d, interpolation='None')
    pylab.show()

  return reference


def integrate(experiments, reflections, reference):
  from dials.algorithms.profile_model.modeller import CircleSampler
  from dials.array_family import flex
  from dials.algorithms.profile_model.gaussian_rs.transform import TransformReverse
  from dials.algorithms.profile_model.gaussian_rs.transform import TransformForward
  from dials.algorithms.profile_model.gaussian_rs.transform import TransformReverseNoModel
  from dials.algorithms.profile_model.gaussian_rs.transform import TransformSpec
  from dials.algorithms.profile_model.gaussian_rs import CoordinateSystem

  selection = reflections.get_flags(reflections.flags.integrated_sum)
  reflections = reflections.select(selection)
  print "Selected %d reflections to integrate" % len(reflections)

  sampler = CircleSampler(
    experiments[0].detector[0].get_image_size(),
    experiments[0].scan.get_array_range(),
    1)


  n_sigma = 4.0
  grid_size = 25
  spec = TransformSpec(
    experiments[0].beam,
    experiments[0].detector,
    experiments[0].goniometer,
    experiments[0].scan,
    experiments[0].profile.sigma_b(deg=False),
    experiments[0].profile.sigma_m(deg=False),
    n_sigma,
    grid_size)

  m2 = experiments[0].goniometer.get_rotation_axis()
  s0 = experiments[0].beam.get_s0()

  Iprf = flex.double(len(reflections))
  Vprf = flex.double(len(reflections))
  Cprf = flex.double(len(reflections))
  Fprf = flex.bool(len(reflections))

  for i, r in enumerate(reflections):
    s1 = r['s1']
    phi = r['xyzcal.mm'][2]
    xyz = r['xyzcal.px']
    bbox = r['bbox']
    panel = r['panel']
    image = r['shoebox'].data.as_double()
    background = r['shoebox'].background.as_double()
    mask = (r['shoebox'].mask.as_1d() == 5)#| (r['shoebox'].mask.as_1d() == 3)
    mask.reshape(image.accessor())
    cs = CoordinateSystem(m2, s0, s1, phi)

    index = sampler.nearest(0, xyz)

    profile = reference[index]

    # print flex.sum(profile)
    # print r['partiality']

    if False:
      from dials.algorithms.integration.maximum_likelihood import  ProfileFittingDouble as ProfileFitting

      transform = TransformReverseNoModel(
        spec,
        cs,
        bbox,
        panel,
        profile)
      p = transform.profile()
      d = image
      m = mask
      b = background
      # print flex.sum(p)

      ysize, xsize = p.all()[1:3]

      p1 = flex.double(flex.grid(1, ysize , xsize))
      d1 = flex.double(flex.grid(1, ysize , xsize))
      b1 = flex.double(flex.grid(1, ysize , xsize))
      m1 = flex.double(flex.grid(1, ysize , xsize))
      for k in range(p.all()[0]):
        p1 += p[k:k+1,:,:]
        d1 += d[k:k+1,:,:]
        b1 += b[k:k+1,:,:]
        m1 = m[k:k+1,:,:]


      try:


        fit = ProfileFitting(p1, m1, d1, b1, 1e-3, 1000)
        assert fit.niter() < 1000
        Iprf[i] = fit.intensity()
        Vprf[i] = fit.variance()
        Cprf[i] = fit.correlation()
        Fprf[i] = True
        print i, fit.intensity(), flex.sum(p1)
        # from matplotlib import pylab
        # pylab.imshow(p1.as_numpy_array()[0,:,:], interpolation='none')
        # pylab.show()
      except Exception:
        pass

    else:
      from dials.algorithms.integration.fit import ProfileFittingDouble as ProfileFitting

      try:

        transform = TransformForward(
          spec,
          cs,
          bbox,
          panel,
          image,
          background,
          mask)

        index = sampler.nearest(0, xyz)

        p = reference[index]
        d = transform.profile()
        b = transform.background()
        m = p > 0


        fit = ProfileFitting(p, m, d, b, 1e-3, 1000)
        assert fit.niter() < 1000
        Iprf[i] = fit.intensity()
        Vprf[i] = fit.variance()
        Cprf[i] = fit.correlation()
        Fprf[i] = True
        print i, fit.intensity(), flex.sum(p)
        # from matplotlib import pylab
        # pylab.imshow(p1.as_numpy_array()[0,:,:], interpolation='none')
        # pylab.show()
      except Exception:
        pass


  reflections['intensity.prf.value'] = Iprf
  reflections['intensity.prf.variance'] = Vprf
  reflections['intensity.prf.correlation'] = Cprf
  reflections.set_flags(Fprf, reflections.flags.integrated_prf)

  return reflections


if __name__ == '__main__':
  from os.path import join
  import sys

  experiments_filename = sys.argv[1]
  reflections_filename = sorted(sys.argv[2:])

  experiments = read_experiments(experiments_filename)
  reflections = read_reflections(reflections_filename)

  print "Read %d reflections" % len(reflections)

  #reference = compute_reference(experiments, reflections)

  import cPickle as pickle
  #pickle.dump(reference, open("reference.pickle", "w"))
  reference = pickle.load(open("reference.pickle"))

  reflections = integrate(experiments, reflections, reference)

  del reflections['shoebox']

  reflections.as_pickle("integrated.pickle")
