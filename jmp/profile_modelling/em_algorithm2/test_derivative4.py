

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


def compute_derivatives(A, B, N, mu, sigma, USE):
  from math import erf, sqrt, exp, pi, log
  erf_A = erf((A[0] - mu) / (sqrt(2) * sigma))
  erf_B = erf((B[-1] - mu) / (sqrt(2) * sigma))
  exp_A = exp(-(A[0] - mu)**2 / (2 * sigma**2))
  exp_B = exp(-(B[-1] - mu)**2 / (2 * sigma**2))
  Z = erf_B - erf_A
  if Z < 1e-10:
    print "ARG"
    return -200, 0, 0

  ntot = sum(N)
  # vtot = ntot / (0.5 * (erf_B - erf_A))
  # dvtot = 2.0*ntot*(sqrt(2/pi)/sigma**2)*((B[-1]-mu)*exp_B - (A[0]-mu)*exp_A) / (erf_B - erf_A)**2

  restart = True
  while restart:
    sum_lnvi = 0
    dsum_lnvi = 0
    dsum_tot = 0
    dsum2 = 0
    dsum_tot2 = 0
    sum_p = 0
    restart = False
    for j, (a, b, n) in enumerate(zip(A, B, N)):
      if USE[j]:
        erf_a = erf((a - mu) / (sqrt(2) * sigma))
        erf_b = erf((b - mu) / (sqrt(2) * sigma))
        exp_a = exp(-(a - mu)**2 / (2 * sigma**2))
        exp_b = exp(-(b - mu)**2 / (2 * sigma**2))
        Z = 0.5 * (erf_B - erf_A)
        kt = sum(N)

        zi = erf_b - erf_a

        if zi < 1e-10:
          print "NOOO"
          zi = 1e-10

        A1 = sqrt(2/pi)*(((b-mu)*exp_b - (a-mu)*exp_a)**2 / sigma**4) / zi**2
        A2 = (2.0*((b-mu)*exp_b - (a-mu)*exp_a)/(sigma**3)) / zi
        A3 = (((b-mu)**3 * exp_b - (a-mu)**3 * exp_a) / sigma**5) / zi

        sum_lnvi += n * (log(kt) + log(zi) - log(erf_B - erf_A))
        dsum_lnvi += n* sqrt(2/pi)*(1/sigma**2)*((b-mu)*exp_b - (a-mu)*exp_a)/zi
        dsum2 += n * sqrt(2/pi) * (A1 - A2 + A3) 

  zi = erf_B - erf_A
  A1 = sqrt(2/pi)*(((B[-1]-mu)*exp_B - (A[0]-mu)*exp_A)**2 / sigma**4) / zi**2
  A2 = (2.0*((B[-1]-mu)*exp_B - (A[0]-mu)*exp_A)/(sigma**3)) / zi
  A3 = (((B[-1]-mu)**3 * exp_B - (A[0]-mu)**3 * exp_A) / sigma**5) / zi

  dsum_tot2 = ntot * sqrt(2/pi) * (A1 - A2 + A3) 
  dsum_tot = ntot* sqrt(2/pi)*(1/sigma**2)*((B[-1]-mu)*exp_B - (A[0]-mu)*exp_A)/zi
  L = sum_lnvi - kt
  dL = dsum_lnvi - dsum_tot
  d2L = dsum2 - dsum_tot2

  return -L, dL, d2L
  
def compute_all_derivatives(A, B, N, I, mu, sigma, USE):
  i0 = 0
  DL = 0
  D2L = 0
  L = 0
  for j in range(len(I)):

    i1 = i0 + I[j]
    AA = A[i0:i1]
    BB = B[i0:i1]
    NN = N[i0:i1]
    UU = USE[i0:i1]

    l,dl,d2l = compute_derivatives(AA, BB, NN, mu, sigma, UU)
    L += l
    DL += dl
    D2L += d2l
    USE[i0:i1] = UU

    i0 = i1
  return L, DL, D2L

def estimate(A, B, N, I, mu, sigma):
  from scipy.optimize import minimize
  
  USE = [True] * len(A)

  def func(sigma):
    L, DL, D2L = compute_all_derivatives(A, B, N, I, mu, sigma, USE)
    return L, DL, D2L


  for i in range(100):
    
    L, DL, D2L = func(sigma)

    delta = DL / abs(D2L)

    from math import pi
    print sigma * 180 /pi, (sigma - delta) * 180 / pi
    sigma -= delta

    if sigma < 0.0001:
       sigma = 0.0001



  # print minimize(func, sigma)

  # print DL, sigma
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

  from math import pi 
  mu = 0
  sigma = 1.0 * pi / 180

 
  A = a_list
  B = b_list
  N = n_list
  I = i_list
  USE = [True] * len(A)

  # X = []
  # Y = []
  # DY = []
  # D2Y = []

  # min_sigma = 2.0
  # max_sigma = 20.0
  # num = 1000
  # for j in range(num):

  #   sigma = min_sigma + j * (max_sigma - min_sigma) / (num - 1)

  #   X.append(sigma)


  #   ntot = sum(N)

  #   L, dL, d2L = compute_all_derivatives(A, B, N, I, mu, sigma, USE)

  #   # Y.append(vtot)
  #   # DY.append(dvtot)
  #   # Y.append(sum_lnvi)
  #   # DY.append(dsum_lnvi)
  #   Y.append(L)
  #   DY.append(dL)
  #   D2Y.append(d2L)

  # print USE.count(False)

  # from numpy import gradient

  # DY2 = gradient(Y, gradient(X))
  # D2Y2 = gradient(DY, gradient(X))

  # from matplotlib import pylab

  # D = [abs(dy) for dy in DY]
  # minimum = X[D.index(min(D))]
  # print minimum

  # # pylab.plot(N)
  # # pylab.ylim((0, max(N)))
  # # pylab.show()

  # pylab.plot(X, Y, color='black')
  # pylab.axvline(minimum)
  # pylab.plot(X, DY, color='blue')
  # # pylab.plot(X, DY2, color='orange')
  # pylab.plot(X, D2Y, color='red')
  # # pylab.plot(X, D2Y2, color='purple')
  # pylab.show()
  
  estimate(A, B, N, I, mu, sigma)
