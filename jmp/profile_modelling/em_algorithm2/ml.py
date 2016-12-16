
from scipy.integrate import romberg

# def romberg(function, a, b, max_steps=10, accuracy=1e-3):
#    Rp = [0] * max_steps
#    Rc = [0] * max_steps
#    h = b-a
#    Rp[0] = (function(a) + function(b))*h*0.5

#    for i in range(1, max_steps):
#       h /= 2.0
#       ep = 1 << (i-1)
#       c = sum(function(a + (2*j-1)*h) for j in range(1, ep+1))
#       Rc[0] = h*c + .5*Rp[0]

#       for j in range(1, i+1):
#          n_k = pow(4, j)
#          Rc[j] = (n_k*Rc[j-1] - Rp[j-1])/(n_k-1)

#       if i > 1 and abs(Rp[i-1] - Rc[i]) < accuracy:
#         return Rc[i-1]

#       Rp, Rc = Rc, Rp

#    return Rp[max_steps-1]



def compute_L(a_list, b_list, n_list, i_list, mu, sigma, USE=[]):
  from math import log

  if len(USE) == 0:
    for i in range(len(a_list)):
      USE.append(True)

  def function(x):
    from math import pi, sqrt, exp
    return (1.0/(sqrt(2*pi)*sigma)) * exp(-(x - mu)**2 / (2 * sigma**2))

  def integral(a, b):
    from math import pi, sqrt, erf
    return 0.5 * (erf((b - mu)/(sqrt(2)*sigma)) - erf((a - mu)/(sqrt(2)*sigma)))

  redo = True
  while redo:
    redo = False
    lnL = 0
    i0 = 0
    for i in range(len(i_list)):
      i1 = i0 + i_list[i]
      A = a_list[i0:i1]
      B = b_list[i0:i1]
      N = n_list[i0:i1]
      sum_p = 0
      tempL = 0
      sum_n = sum(N)
      Pr = integral(A[0], B[-1])#1.0#romberg(function, A[0], B[-1])
      Pv = 1.0 - Pr
      # M = sum_n * res / tot
      # norm = sum_n * (1.0 + res / tot)
      # norm = sum_n / tot
      vtot = sum_n * (1.0 + Pv / Pr)
      # if tot < 1e-3:
      #   print sum_n, list(A), list(B), tot, norm
      for j, (a, b, n) in enumerate(zip(A, B, N)):
        if USE[i0+j]:
          p = vtot * integral(a, b)#romberg(function, a, b)
          if p < 1e-10:
            print "NOO"
            USE[i0+j] = False
            redo = True
            break
          tempL += n * log(p)
          # sum_p += p
      tempL -= vtot
      # tempL -= sum_p
      lnL += tempL
      i0 = i1
  return -lnL

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

  reflections = flex.reflection_table.from_pickle(sys.argv[1])
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

  # selection = flex.abs(zeta) > 0.05
  # print selection.count(False)
  # reflections = reflections.select(selection)


  #print "Num Refl: ", len(reflections)

  a_list = []
  b_list = []
  n_list = []
  i_list = []
  for p, s, z in zip(phi, sbox, zeta):
    z0 = s.bbox[4]

    # if s.data.all()[0] == 1:
    #   continue

    a_temp = []
    b_temp = []
    n_temp = []
    for k in range(s.data.all()[0]):
      phi0 = scan.get_angle_from_array_index(z0+k, deg=False)
      phi1 = scan.get_angle_from_array_index(z0+k+1, deg=False)

      if z < 0:
        a = (p - phi0) * z#abs(z)
        b = (p - phi1) * z#abs(z)
      else:
        b = (p - phi0) * z#abs(z)
        a = (p - phi1) * z#abs(z)

      sum_frames = 0
      for j in range(s.data.all()[1]):
        for i in range(s.data.all()[2]):
          if s.mask[k,j,i] != 0:
            sum_frames += s.data[k,j,i]
            #n_list.append(s.data[k,j,i])
      n_temp.append(sum_frames)
      a_temp.append(a)
      b_temp.append(b)
      assert a < b

    if z > 0:
      a_temp = list(reversed(a_temp))
      b_temp = list(reversed(b_temp))

    assert all(abs(bb-aa) < 1e-7 for aa, bb in zip(a_temp[1:], b_temp[:-1]))

    a_list.extend(a_temp)
    b_list.extend(b_temp)
    n_list.extend(n_temp)
    i_list.append(s.data.all()[0])

  a_list = flex.double(a_list)
  b_list = flex.double(b_list)
  n_list = flex.double(n_list)
  i_list = flex.size_t(i_list)

  assert sum(i_list) == len(a_list)

  mid = (a_list + b_list) / 2.0

  index = sorted(range(len(mid)), key=lambda x: mid[x])

  K = 1.5
  N = len(index)
  I1 = index[int(N*0.25)]
  I3 = index[int(N*0.75)]
  Q1 = mid[I1]
  Q3 = mid[I3]
  IQR = Q3 - Q1
  L = Q1 - K*IQR
  R = Q3 + K*IQR

  # USE = [m > L and m < R for m in mid]
  USE = [True for m in mid]


  ###########

  from math import pi
  mu = 0
  x = []
  y = []
  min_x = 0.05* pi / 180.0
  max_x = 0.4 * pi / 180.0
  num = 100
  for i in range(num):
    sigma = min_x + i * (max_x - min_x) / float(num-1)
    L = compute_L(a_list, b_list, n_list, i_list, mu, sigma, USE)
    print sigma * 180/pi, L, i_list.count(1), i_list.count(2), i_list.count(3), len(i_list)
    x.append(sigma)
    y.append(L)

  from matplotlib import pylab
  pylab.plot([xx*180.0/pi for xx in x], y)
  pylab.show()
