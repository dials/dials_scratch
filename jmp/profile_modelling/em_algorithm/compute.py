
def compute_centroid(a, b, p):
  from dials.array_family import flex
  from math import sqrt
  x = (a + b) / 2
  xc = flex.sum(p * x) / flex.sum(p)
  xv = flex.sum(p * (x - xc)**2) / (flex.sum(p)-1)
  return xc, sqrt(xv)


def estimate_parameters(a, b, n, index, mean0, sigma0i, tolerance=1e-7):
  from math import sqrt, pi, erf, exp, log

  TINY = 1e-10
  d = [bb - aa for aa, bb in zip(a, b)]
  min_d = min(d)
  X1 = max(b)
  X0 = min(a)
  fi0 = min_d / (X1 - X0)
  #pi0 = TINY / fi0
  pi1 = 1.0# - pi0

  def expectation(a, b, mu, si):
    e1 = erf((b-mu)/(sqrt(2)*si)) 
    e2 = erf((a-mu)/(sqrt(2)*si))
    e3 = exp(-(a-mu)**2 / (2 * si**2)) / (sqrt(2.0*pi)*si)
    e4 = exp(-(b-mu)**2 / (2 * si**2)) / (sqrt(2.0*pi)*si)
    exp0 = 0.5 * (e1 - e2)
    exp1 = exp0 * mu + (si*si)*(e3 - e4)
    exp2 = exp0 * si*si + si*si*((a-mu)*e3 - (b-mu)*e4)
    return exp0, exp1, exp2


  USE = [True] * len(a)

  AA = a
  BB = b
  NN = n

  mu = mean0
  si = sigma0
  lnL0 = 0
  for i in range(1000):

    mu_sum = 0
    va_sum = 0
    #c0_sum = 0
    c1_sum = 0
    m_sum = 0
    lnL_sum = 0
    i0 = 0
    for l in range(len(index)):
      i1 = i0 + index[l]
      a = AA[i0:i1]
      b = BB[i0:i1]
      n = NN[i0:i1]

      E0j, E1j, E2j = [], [], []
      for j in range(len(a)):
        e0, e1, e2 = expectation(a[j], b[j], mu, si)
        E0j.append(e0)
        E1j.append(e1)
        E2j.append(e2)

      m = n

      #Pj0 = [pi0 * (b[j]-a[j]) / (X1 - X0) for j in range(len(a))]
      Pj1 = [pi1 * E0j[j] for j in range(len(a))]
      #Pj = [pj0 + pj1 for pj0, pj1 in zip(Pj0, Pj1)]
      Pj = Pj1

      for j in range(len(a)):
        if USE[i0+j] and Pj[j] < 1e-10:
          #print "RESTART"
          USE[i0+j] = False


      U = [USE[i0+j] for j in range(len(a))]
      if not any(U):
        continue
      
      P = sum(Pj[j] for j in range(len(a)) if USE[i0+j])
      N = sum(n[j] for j in range(len(a)) if USE[i0+j])

      if P > 1.0:
        P = 1.0
      Pv = 1.0 - P
      Pv1 = Pv
      assert Pv >= 0
      mv = N * Pv / P
      E1v = mu - sum(E1j)
      E2v = si*si - sum(E2j)
      #c0_sum += sum(mm * pi0 * pj0 / pj for mm, pj0, pj in zip(m, Pj0, Pj))# +mv*pi0*(1.0-sum(Pj0))/Pv

      for j in range(len(a)):
        if USE[i0+j]:
          c1_sum += m[j] * pi1 * Pj1[j] / Pj[j]
          mu_sum += m[j] * pi1 * E1j[j] / Pj[j]
          va_sum += m[j] * pi1 * E2j[j] / Pj[j]
          m_sum += m[j]
          lnL_sum += n[j] * log(Pj[j])
      lnL_sum -= N * log(P)

      if (Pv > 1e-10):
        c1_sum += mv*pi1*Pv1 / Pv
        mu_sum += mv*pi1*E1v / Pv
        va_sum += mv*pi1*E2v / Pv
        m_sum += mv
      
      i0 = i1

    #print c0_sum, c1_sum

    #c0 = c0_sum
    c1 = c1_sum
    #pi0 = c0 / m_sum
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
    if (i > 0 and abs((lnL - lnL0) / lnL0) < tolerance):
      break
    #print P, pi0, pi1, lnL, lnL0, lnL - lnL0, mu, si * 180 / pi
    
    lnL0 = lnL
  
  return mu, si

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
    
    for k in range(s.data.all()[0]):
      phi0 = scan.get_angle_from_array_index(z0+k, deg=False)
      phi1 = scan.get_angle_from_array_index(z0+k+1, deg=False)

      a = (p - phi0) * z#abs(z)
      b = (p - phi1) * z#abs(z)
      if a > b:
        a, b = b, a

      sum_frames = 0
      for j in range(s.data.all()[1]):
        for i in range(s.data.all()[2]):
          if s.mask[k,j,i] != 0:
            sum_frames += s.data[k,j,i]
            #n_list.append(s.data[k,j,i])
      n_list.append(sum_frames)
      a_list.append(a)
      b_list.append(b)
      assert a < b

    i_list.append(s.data.all()[0])

  a_list = flex.double(a_list)
  b_list = flex.double(b_list)
  n_list = flex.double(n_list)
  i_list = flex.size_t(i_list)

  assert sum(i_list) == len(a_list)

  from math import pi


  from dials.algorithms.statistics import BinnedGMMSingle1DFixedMean, BinnedGMMSingle1D

  from time import time



  mean0, sigma0 = compute_centroid(a_list, b_list, n_list)
  #print "Mean0, Sigma0", mean0, sigma0 * 180 / pi

  mean = mean0
  sigma = sigma0

  mean, sigma = estimate_parameters(a_list, b_list, n_list, i_list, mean0, sigma0)

  #print "Mean1, Sigma1", mean, sigma * 180 / pi
  print sigma0 * 180.0 / pi, sigma * 180.0 / pi


