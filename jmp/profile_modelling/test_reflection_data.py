

def gaussian(A, B, N, mu, sigma):
  from math import sqrt, erf
  from numpy.random import poisson
  G = []
  for a, b in zip(A,B):
    G.append(
      poisson(
        N*0.5*(
          erf((b-mu)/(sqrt(2)*sigma)) -
          erf((a-mu)/(sqrt(2)*sigma)))))
  return G



def select_reflections(reference, experiments):
  ''' Load the reference spots. '''
  from dials.array_family import flex
  from time import time
  from libtbx.utils import Sorry
  if reference is None:
    return None, None
  st = time()
  assert("miller_index" in reference)
  assert("id" in reference)
  print('Processing reference reflections')
  print(' read %d strong spots' % len(reference))
  mask = reference.get_flags(reference.flags.indexed)
  rubbish = reference.select(mask == False)
  if mask.count(False) > 0:
    reference.del_selected(mask == False)
    print(' removing %d unindexed reflections' %  mask.count(False))
  if len(reference) == 0:
    raise Sorry('''
      Invalid input for reference reflections.
      Expected > %d indexed spots, got %d
    ''' % (0, len(reference)))
  mask = reference.get_flags(reference.flags.centroid_outlier)
  if mask.count(True) > 0:
    rubbish.extend(reference.select(mask))
    reference.del_selected(mask)
    print(' removing %d reflections marked as centroid outliers' %  mask.count(True))
  mask = reference['miller_index'] == (0, 0, 0)
  if mask.count(True) > 0:
    rubbish.extend(reference.select(mask))
    reference.del_selected(mask)
    print(' removing %d reflections with hkl (0,0,0)' %  mask.count(True))
  mask = reference['id'] < 0
  if mask.count(True) > 0:
    raise Sorry('''
      Invalid input for reference reflections.
      %d reference spots have an invalid experiment id
    ''' % mask.count(True))
  print(' using %d indexed reflections' % len(reference))
  print(' found %d junk reflections' % len(rubbish))
  print(' time taken: %g' % (time() - st))

  predicted = flex.reflection_table.from_predictions_multi(
    experiments)

  # Match reference with predicted
  matched, reference, unmatched = predicted.match_with_reference(reference)
  assert(len(matched) == len(predicted))
  assert(matched.count(True) <= len(reference))
  return reference

if __name__ == '__main__':
  from dials.array_family import flex
  from dxtbx.model.experiment.experiment_list import ExperimentListFactory
  import sys

  from truncated_normal import compute_all_derivatives, estimate
  from normal import compute_all_derivatives, estimate
  from normal_known_scale import compute_all_derivatives, estimate

  print "Load reflections"
  reflections = flex.reflection_table.from_pickle(sys.argv[1])

  print "Load experiments"
  experiments = ExperimentListFactory.from_json_file(sys.argv[2],
                                                     check_format=False)

  beam = experiments[0].beam
  detector = experiments[0].detector
  goniometer = experiments[0].goniometer
  scan = experiments[0].scan

  reflections = select_reflections(reflections, experiments)

  # selection = reflections.get_flags(reflections.flags.used_in_refinement)
  # reflections = reflections.select(selection)
  reflections.compute_zeta(experiments[0])

  phi = reflections['xyzcal.mm'].parts()[2]
  sbox = reflections['shoebox']
  zeta = reflections['zeta']

  dphi2 = scan.get_oscillation(deg=False)[1] / 2.0
  # selection = flex.abs(zeta) > 0.05
  # print selection.count(False)
  # reflections = reflections.select(selection)


  #print "Num Refl: ", len(reflections)

  reflections.compute_background(experiments)
  print "Get data"
  a_list = []
  b_list = []
  n_list = []
  i_list = []
  for i, (p, s, z) in enumerate(zip(phi, sbox, zeta)):
    print i, len(phi)
    z0 = s.bbox[4]
    z1 = s.bbox[5]
    phi0 = scan.get_angle_from_array_index(z0, deg=False)
    phi1 = scan.get_angle_from_array_index(z1, deg=False)
    # if phi0 > p or phi1 < p or (z1 -z0) == 1:
    #   print "Skipping"
    #   continue

    # if s.data.all()[0] == 1:
    #   continue

    a_temp = []
    b_temp = []
    n_temp = []
    for k in range(s.data.all()[0]):
      phi0 = scan.get_angle_from_array_index(z0+k, deg=False)
      phi1 = scan.get_angle_from_array_index(z0+k+1, deg=False)

      t = (phi1 + phi0)/2.0 - p
      b = (t + dphi2) * abs(z)
      a = (t - dphi2) * abs(z)
      sum_frames = 0
      for j in range(s.data.all()[1]):
        for i in range(s.data.all()[2]):
          if s.mask[k,j,i] != 0:
            sum_frames += s.data[k,j,i]# - s.background[k,j,i]
            #n_list.append(s.data[k,j,i])
      # if sum_frames <= 0:
      #   sum_frames=1
      n_temp.append(sum_frames)
      a_temp.append(a)
      b_temp.append(b)
      assert a < b

    # a_temp = list(reversed(a_temp))
    # b_temp = list(reversed(b_temp))

    assert all(abs(bb-aa) < 1e-7 for aa, bb in zip(a_temp[1:], b_temp[:-1]))
      # phi0 = scan.get_angle_from_array_index(z0+k, deg=False)
      # phi1 = scan.get_angle_from_array_index(z0+k+1, deg=False)

      # if z < 0:
      #   a = (p - phi0) * z#abs(z)
      #   b = (p - phi1) * z#abs(z)
      # else:
      #   b = (p - phi0) * z#abs(z)
      #   a = (p - phi1) * z#abs(z)

      # sum_frames = 0
      # for j in range(s.data.all()[1]):
      #   for i in range(s.data.all()[2]):
      #     if s.mask[k,j,i] != 0:
      #       sum_frames += s.data[k,j,i]
      #       #n_list.append(s.data[k,j,i])
      # n_temp.append(sum_frames)
      # a_temp.append(a)
      # b_temp.append(b)
      # assert a < b

    # if z > 0:
      # a_temp = list(reversed(a_temp))
      # b_temp = list(reversed(b_temp))

    # assert all(abs(bb-aa) < 1e-7 for aa, bb in zip(a_temp[1:], b_temp[:-1]))

    a_list.extend(a_temp)
    b_list.extend(b_temp)
    n_list.extend(n_temp)
    i_list.append(s.data.all()[0])

  a_list = flex.double(a_list)
  b_list = flex.double(b_list)
  n_list = flex.double(n_list)
  i_list = flex.size_t(i_list)

  print len(a_list)

  assert sum(i_list) == len(a_list)
  from math import pi
  mu = 0
  sigma = 1.0 * pi / 180


  A = a_list
  B = b_list
  N = n_list
  I = i_list

  if True:

    #USE = [True] * len(A)

    X = []
    Y = []
    DY = []
    D2Y = []

    min_sigma = 0.01 * pi / 180.0
    max_sigma = 0.3 * pi / 180.0
    num = 100
    for j in range(num):
      USE = [True] * len(A)
      print j
      sigma = min_sigma + j * (max_sigma - min_sigma) / (num - 1)

      X.append(sigma*180/pi)


      ntot = sum(N)

      L, dL, d2L = compute_all_derivatives(A, B, N, I, mu, sigma, USE)

      # Y.append(vtot)
      # DY.append(dvtot)
      # Y.append(sum_lnvi)
      # DY.append(dsum_lnvi)
      Y.append(L)
      DY.append(dL)
      D2Y.append(d2L)

    # with open("sum.txt", "w") as outfile:
    #   for x, y in zip(X,Y):
    #     print >>outfile, x, y
    # exit(0)

    print USE.count(False)

    from numpy import gradient

    DY2 = gradient(Y, gradient(X))
    D2Y2 = gradient(DY, gradient(X))

    from matplotlib import pylab

    D = [abs(dy) for dy in DY]
    minimum = X[Y.index(min(Y))]
    print minimum

    # pylab.plot(N)
    # pylab.ylim((0, max(N)))
    # pylab.show()

    pylab.plot(X, Y, color='black')
    pylab.axvline(minimum)
    # pylab.plot(X, DY, color='blue')
    # pylab.plot(X, DY2, color='orange')
    # pylab.plot(X, D2Y, color='red')
    # pylab.plot(X, D2Y2, color='purple')
    pylab.show()
  exit(0)
  print estimate(A, B, N, I, mu, 0.01*pi/180, 2.0*pi/180) * 180.0/pi
