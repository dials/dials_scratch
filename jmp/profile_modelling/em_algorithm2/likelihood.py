
def regrid(a_list, b_list, n_list):
  from math import floor, ceil
  min_a = min(a_list)
  max_b = max(b_list)
  num_bins = 100
  s1 = (num_bins-1) / (max_b - min_a)
  s0 = -s1 * min_a

  x1 = [min_a + i * (max_b-min_a) / num_bins for i in range(num_bins)]
  x2 = [min_a + (i+1) * (max_b-min_a) / num_bins for i in range(num_bins)]
  y = [0] * num_bins
  for a, b, n in zip(a_list, b_list, n_list):
    i0 = s0 + s1 * a
    i1 = s0 + s1 * b
    assert i0 < i1
    for i in range(int(floor(i0)), int(ceil(i1))):
      j0 = max(i0, i)
      j1 = min(i1, i+1)
      if j0 < j1:
        fraction = (j1 - j0) / (i1 - i0)
        y[i] += n*fraction
  assert abs(sum(y) - sum(n_list)) < 1e-7
  return x1, x2, y

def compute_centroid(a, b, p):
  from dials.array_family import flex
  from math import sqrt
  x = (a + b) / 2
  xc = flex.sum(p * x) / flex.sum(p)
  xv = flex.sum(p * (x - xc)**2) / (flex.sum(p)-1)
  return xc, sqrt(xv)


def compute_likelihood(a, b, n, index, mean0, sigma0, USE=None):
  from math import sqrt, pi, erf, exp, log

  TINY = 1e-10
  d = [bb - aa for aa, bb in zip(a, b)]
  min_d = min(d)
  X1 = max(b)
  X0 = min(a)
  fi0 = min_d / (X1 - X0)
  pi0 = TINY / fi0
  pi1 = 1.0 - pi0

  def expectation(a, b, mu, si):
    e1 = erf((b-mu)/(sqrt(2)*si))
    e2 = erf((a-mu)/(sqrt(2)*si))
    e3 = exp(-(a-mu)**2 / (2 * si**2)) / (sqrt(2.0*pi)*si)
    e4 = exp(-(b-mu)**2 / (2 * si**2)) / (sqrt(2.0*pi)*si)
    exp0 = 0.5 * (e1 - e2)
    exp1 = exp0 * mu + (si*si)*(e3 - e4)
    exp2 = exp0 * si*si + si*si*((a-mu)*e3 - (b-mu)*e4)
    return exp0, exp1, exp2

  if USE is None:
    USE = [True] * len(a)

  AA = a
  BB = b
  NN = n

  mu = mean0
  si = sigma0

  mu_sum = 0
  va_sum = 0
  c0_sum = 0
  c1_sum = 0
  m_sum = 0
  lnL_sum = 0
  i0 = 0
  for l in range(len(index)):
    i1 = i0 + index[l]
    a = AA[i0:i1]
    b = BB[i0:i1]
    n = NN[i0:i1]
    assert all(aa < bb for aa, bb in zip(a, b))
    assert all(abs(bb - aa) < 1e-7 for aa, bb in zip(a[1:], b[:-1]))

    E0j, E1j, E2j = [], [], []
    for j in range(len(a)):
      e0, e1, e2 = expectation(a[j], b[j], mu, si)
      E0j.append(e0)
      E1j.append(e1)
      E2j.append(e2)

    m = n

    Pj0 = [pi0 * (b[j]-a[j]) / (X1 - X0) for j in range(len(a))]
    Pj1 = [pi1 * E0j[j] for j in range(len(a))]
    Pj = [pj0 + pj1 for pj0, pj1 in zip(Pj0, Pj1)]
    #Pj = Pj1

    for j in range(len(a)):
      if USE[i0+j] and Pj[j] < 1e-12:
        #print "RESTART"
        USE[i0+j] = False
        exit(0)

    U = [USE[i0+j] for j in range(len(a))]
    if any(U):

      P0 = sum(Pj0[j] for j in range(len(a)) if USE[i0+j])
      P1 = sum(Pj1[j] for j in range(len(a)) if USE[i0+j])
      P = sum(Pj[j] for j in range(len(a)) if USE[i0+j])
      N = sum(n[j] for j in range(len(a)) if USE[i0+j])

      assert 1.0 - P > -1e-7
      if P > 1.0:
        P = 1.0
      Pv = 1.0 - P
      Pv1 = 1.0 - P1
      Pv0 = 1.0 - P0
      assert Pv >= 0
      mv = N * Pv / P
      E1v = mu - sum(E1j)
      E2v = si*si - sum(E2j)

      for j in range(len(a)):
        if USE[i0+j]:
          c0_sum += m[j] * pi0 * Pj0[j] / Pj[j]
          c1_sum += m[j] * pi1 * Pj1[j] / Pj[j]
          mu_sum += m[j] * pi1 * E1j[j] / Pj[j]
          va_sum += m[j] * pi1 * E2j[j] / Pj[j]
          m_sum += m[j]
          lnL_sum += n[j] * log(Pj[j])
      lnL_sum -= N * log(P)

      c0_sum += N*pi0*Pv0 / P
      c1_sum += N*pi1*Pv1 / P #mv*pi1*Pv1 / Pv
      mu_sum += N*pi1*E1v / P #mv*pi1*E1v / Pv
      va_sum += N*pi1*E2v / P #mv*pi1*E2v / Pv
      m_sum += mv

    i0 = i1
  assert i0 == len(USE)


  #print c0_sum, c1_sum

  c0 = c0_sum
  c1 = c1_sum
  pi0 = c0 / m_sum
  pi1 = c1 / m_sum
  mu = mu_sum / c1_sum
  va = va_sum / c1_sum
  si = sqrt(va)

  # c0 = sum(mm * pi0 * pj0 / pj for mm, pj0, pj in zip(m, Pj0, Pj))
  # c1 = sum(mm * pi1 * pj1 / pj for mm, pj1, pj in zip(m, Pj1, Pj))
  # pi0 = c0 / sum(m)
  # pi1 = c1 / sum(m)
  # mu = sum(m[j] * pi1 * E1j[j] / Pj[j] for j in range(len(a))) / c1
  # va = sum(m[j] * pi1 * E2j[j] / Pj[j] for j in range(len(a))) / c1
  # si = sqrt(va)
  lnL = lnL_sum

  #print lnL0, lnL, lnL - lnL0, mu, si * 180 / pi
  #if (i > 0 and abs((lnL - lnL0) / lnL0) < tolerance):
  #  break
  #print P, pi0, pi1, lnL, lnL0, lnL - lnL0, mu, si * 180 / pi

  return lnL, USE

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

  selection = reflections.get_flags(reflections.flags.used_in_refinement)
  reflections = reflections.select(selection)
  reflections.compute_zeta(experiments[0])

  phi = reflections['xyzcal.mm'].parts()[2]
  sbox = reflections['shoebox']
  zeta = reflections['zeta']

  #print "Num Refl: ", len(reflections)

  a_list = []
  b_list = []
  n_list = []
  i_list = []
  for p, s, z in zip(phi, sbox, zeta):
    z0 = s.bbox[4]

#     if s.data.all()[0] == 1:
#       continue

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


  from dials.algorithms.statistics import BinnedGMMSingle1DFixedMean, BinnedGMMSingle1D

  from time import time

  # a_list, b_list, n_list = regrid(a_list, b_list, n_list)
  # i_list = [len(a_list)]
  # a_list = flex.double(a_list)
  # b_list = flex.double(b_list)
  # n_list = flex.double(n_list)
  # i_list = flex.size_t(i_list)

  # mean0, sigma0 = compute_centroid(a_list, b_list, n_list)
  # print sigma0 * 180.0 / pi
  # exit(0)

  x = []
  y = []
  min_x = 0.08* pi / 180.0
  max_x = 0.15 * pi / 180.0
  num = 10
  USE = None
  for i in range(num):
    sigma = min_x + i * (max_x - min_x) / float(num-1)
    L, USE = compute_likelihood(a_list, b_list, n_list, i_list, 0.0, sigma,
                                USE=None)
    print sigma * 180/pi, L, USE.count(True), len(USE), i_list.count(1), len(i_list)
    x.append(sigma)
    y.append(L)

  from matplotlib import pylab
  pylab.plot([xx*180.0/pi for xx in x], y)
  pylab.show()