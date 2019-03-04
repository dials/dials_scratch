from __future__ import print_function


def compute_profile(experiments, reflection, reference, N):
  from dials.array_family import flex
  from dials.algorithms.profile_model.gaussian_rs import CoordinateSystem
  from dials.algorithms.profile_model.modeller import GridSampler
  from dials_scratch.jmp.sim import compute_profile_internal
  from random import uniform

  sbox = reflection['shoebox']
  bbox = sbox.bbox
  zs = sbox.zsize()
  ys = sbox.ysize()
  xs = sbox.xsize()

  profile = flex.double(flex.grid(zs, ys, xs))

  m2 = experiments[0].goniometer.get_rotation_axis_datum()
  s0 = experiments[0].beam.get_s0()
  s1 = reflection['s1']
  phi = reflection['xyzcal.mm'][2]
  detector = experiments[0].detector
  scan = experiments[0].scan

  cs = CoordinateSystem(m2, s0, s1, phi)

  scan_range = scan.get_array_range()
  image_size = detector[0].get_image_size()
  grid_size = (3, 3, 40)
  assert grid_size[0]*grid_size[1]*grid_size[2] == len(reference[0])

  sampler = GridSampler(image_size, scan_range, grid_size)

  xyz = reflection['xyzcal.px']
  index = sampler.nearest(0, xyz)

  for g in reference[0]:
    assert abs(flex.sum(g) - 1.0) < 1e-7

  grid = reference[0][index]

  sigma_d = experiments[0].profile.sigma_b(deg=False)
  sigma_m = experiments[0].profile.sigma_m(deg=False)
  delta_d = 3.0 * sigma_d
  delta_m = 3.0 * sigma_m

  profile = compute_profile_internal(
    grid,
    bbox,
    zs,
    ys,
    xs,
    N,
    delta_d,
    delta_m,
    detector,
    scan,
    cs)

  # from dials_scratch.jmp.viewer import show_image_stack_multi_view
  # show_image_stack_multi_view(profile.as_numpy_array(), vmax=max(profile))
  sum_p = flex.sum(profile)
  print("Partiality: %f" % sum_p)
  try:
    assert sum_p > 0, "sum_p == 0"
  except Exception as e:
    print(e)
    return None

  return profile




def write_profiles(experiments, reflections, reference, N, filename):

  profiles = []

  for i, refl in enumerate(reflections):

    int_flag = refl['flags'] & reflections.flags.integrated_prf

    if int_flag:
      print("Computing profile for %d / %d" % (i, len(reflections)))
      s = compute_profile(experiments, refl, reference, N)
    else:
      s = None

    profiles.append(s)

  import cPickle as pickle
  pickle.dump(profiles, open(filename, "w"))



from libtbx.phil import parse

phil_scope = parse('''

  experiments = None
    .type = str

  reflections = NOne
    .type = str

  reference = None
    .type = str

  N = 125
    .type = int

  output = "profile.pickle"
    .type = str

''')


if __name__ == '__main__':

  from dials.array_family import flex
  from dxtbx.model.experiment_list import ExperimentListFactory
  from dials.algorithms.integration.report import IntegrationReport
  import cPickle as pickle
  import sys

  # Parse command line
  interpretor = phil_scope.command_line_argument_interpreter()
  user_phils = interpretor.process_args(sys.argv[1:])
  working_phil = phil_scope.fetch(sources=user_phils)
  params = working_phil.extract()

  assert params.experiments is not None
  assert params.reflections is not None
  assert params.reference is not None

  # Extract experiments and reflections
  experiments = ExperimentListFactory.from_json_file(params.experiments)
  reflections = flex.reflection_table.from_pickle(params.reflections)
  reference = pickle.load(open(params.reference))

  write_profiles(experiments, reflections, reference, params.N, params.output)
