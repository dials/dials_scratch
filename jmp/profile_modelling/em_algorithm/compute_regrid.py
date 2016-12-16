
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
  import pickle
  from math import pi

  A, B, N = pickle.load(open(sys.argv[1]))
  A = flex.double(A)
  B = flex.double(B)
  N = flex.double(N)

  mean0, sigma0 = compute_centroid(A, B, N)
  #print "Mean0, Sigma0", mean0, sigma0 * 180 / pi

  mean = mean0
  sigma = sigma0

  mean, sigma = estimate_parameters(A, B, N, [len(A)], mean0, sigma0)

  #print "Mean1, Sigma1", mean, sigma * 180 / pi
  print sigma0, sigma


