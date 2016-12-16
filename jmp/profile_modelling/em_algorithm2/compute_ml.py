
from scipy.integrate import romberg



def minimize(a_list, b_list, n_list, i_list, mu, sigma):
  from math import pi, exp, sqrt, erf

  # def derivative_vtot(ntot, A, B, mu, sigma):
  #   X1 = -2 * ntot * (sqrt(2) / (sqrt(pi) * sigma**2))
  #   X2 = (A - mu)*exp(-(A - mu)**2 / (2*sigma**2))
  #   X3 = (B - mu)*exp(-(B - mu)**2 / (2*sigma**2))
  #   X4 = erf((B - mu)/(sqrt(2)*sigma))
  #   X5 = erf((A - mu)/(sqrt(2)*sigma))
  #   return X1 * (X2 - X3) / (X4 - X5)


  # def derivative_f(ntot, A, B, mu, sigma):
  #   X1 = (1.0 / (sqrt(2*pi) * sigma**2))
  #   X2 = (A - mu)*exp(-(A - mu)**2 / (2*sigma**2))
  #   X3 = (B - mu)*exp(-(B - mu)**2 / (2*sigma**2))
  #   return X1*(X2 - X3)

  # def derivative_vi(ntot, A, B, mu, sigma):
  #   dvtot = derivative_vtot(ntot, A, B, mu, sigma)
  #   vtot = function_vtot(ntot, A, B, mu, sigma)
  #   df = derivative_f(ntot, A, B, mu, sigma)
  #   f = function_f(ntot, A, B, mu, sigma)
  
  # def function_vtot(ntot, A, B, mu, sigma):
  #   Pr = integral(A[0], B[-1])
  #   Pv = 1.0 - Pr
  #   return ntot * (1.0 + Pv / Pr)

  # def function_vi(ntot, A, B, mu, sigma):
  #   return function_vtot(ntot, A, B, mu, sigma) * 0.5 * (erf((B - mu)/(sqrt(2)*sigma)) - erf((A - mu)/(sqrt(2)*sigma)))
  
  # def derivative_lnvi(ntot, A, B, mu, sigma):
  #   return derivative_vi(ntot, A, B, mu, sigma) / function_vi(ntot, A, B, mu, sigma)

  while True:
    restart = False
    i0 = 0
    sum_dvtot = 0
    sum_dlnvi = 0
    for i in range(len(i_list)):
      i1 = i0 + i_list[i]
      A = a_list[i0:i1]
      B = b_list[i0:i1]
      N = n_list[i0:i1]


      erf_A = erf((A[0] - mu) / (sqrt(2) * sigma))
      erf_B = erf((B[-1] - mu) / (sqrt(2) * sigma))
      exp_A = exp(-(A[0] - mu)**2 / (2 * sigma**2))
      exp_B = exp(-(B[-1] - mu)**2 / (2 * sigma**2))
  
      ntot = sum(N)
      vtot = ntot / (0.5 * (erf_B - erf_A))
      dvtot = 2.0*ntot*(sqrt(2/pi)/sigma**2)*((B[-1]-mu)*exp_B - (A[0]-mu)*exp_A) / (erf_B - erf_A)**2
      
      sum_dvtot += dvtot

      for j, (a, b, n) in enumerate(zip(A, B, N)):
        if USE[i0+j]:
          erf_a = erf((a - mu) / (sqrt(2) * sigma))
          erf_b = erf((b - mu) / (sqrt(2) * sigma))
          exp_a = exp(-(a - mu)**2 / (2 * sigma**2))
          exp_b = exp(-(b - mu)**2 / (2 * sigma**2))
          p = 0.5 * (erf_b - erf_a)
          dp = -(1.0 / (sqrt(2*pi) *sigma**2))*((b-mu)*exp_b - (a-mu)*exp_a)
          vi = vtot * p
          dvi = vtot * dp + p * dvtot
          if vi < 1e-10:
            USE[i0+j] = False
            restart = True
            break
          dlnvi = dvi / vi

          sum_dlnvi += n * dlnvi

      i0 = i1

    if restart:
      continue

    df = sum_dlnvi - sum_dvtot
    print sum_dlnvi, sum_dvtot, df

    sigma += -0.1 * df
      
    print sigma * 180 / pi

  return -lnL

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

  from math import pi

  minimize(a_list, b_list, n_list, i_list, 0, 0.1 * pi / 180.0)

  ###########

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

