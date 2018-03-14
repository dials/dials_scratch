from __future__ import division
import numpy.random
from scitbx import matrix
from math import sqrt, pi, sin, cos, log, exp, ceil
from dials_scratch.jmp.potato.util.simplex import SimpleSimplex
from dials_scratch.jmp.potato.util.generate_simple import generate_from_reflections
from dials_scratch.jmp.potato.util.generate_simple import generate_from_reflections_binned
from dials_scratch.jmp.potato.model import compute_change_of_basis_operation
from dxtbx.model.experiment_list import ExperimentListFactory
from dials.array_family import flex


def log_likelihood(params, s0, s2_list, xbar_list, ctot_list, Sobs_list, test=0):


  M = matrix.sqr((
    params[0], 0, 0,
    params[1], params[2], 0,
    params[3], params[4], params[5]))

  sigma = M*M.transpose()

  lnL = 0
  for i in range(len(s2_list)):

    s2 = s2_list[i]
    xbar = xbar_list[i]
    ctot = ctot_list[i]
    Sobs = Sobs_list[i]

    R = compute_change_of_basis_operation(s0, s2)
    S = R*sigma*R.transpose()
    S11 = matrix.sqr((
      S[0], S[1],
      S[3], S[4]))
    S12 = matrix.col((S[2], S[5]))
    S21 = matrix.col((S[6], S[7])).transpose()
    S22 = S[8]
    mu = R*s2
    assert abs(1 - mu.normalize().dot(matrix.col((0, 0, 1)))) < 1e-7
    mu1 = matrix.col((mu[0], mu[1]))
    mu2 = mu[2]

    mubar = mu1 + S12*(1/S22)*(s0.length()-mu2)
    Sbar = S11 - S12*(1/S22)*S21
    try:
      Sobs = matrix.sqr(Sobs)
      B1 = log(S22)
      B2 = (1/S22)*(s0.length()-mu2)**2
      A1 = log(Sbar.determinant())*ctot
      A2 = ((Sbar.inverse())*(Sobs)).trace()
      A3 = ((Sbar.inverse())*(ctot*(xbar-mubar)*(xbar-mubar).transpose())).trace()
      if test:
        lnL += -0.5*(A1+A2+A3+ctot*(B1+B2))
      else:
        lnL += -0.5*(A1+A2+A3+(B1+B2))
    except Exception:
      lnL += -1e15

  print tuple(sigma), lnL
  return lnL


class Target(object):

  def __init__(self,
               s0,
               s2_list,
               xbar_list,
               ctot_list,
               Sobs_list,
               test=1):
    self.s0 = s0
    self.s2_list = s2_list
    self.xbar_list = xbar_list
    self.ctot_list = ctot_list
    self.Sobs_list = Sobs_list
    self.test = test

  def target(self, params):

    lnL = log_likelihood(
      params,
      self.s0,
      self.s2_list,
      self.xbar_list,
      self.ctot_list,
      self.Sobs_list,
      test=self.test)

    score = -lnL
    return score


def generate_observations2(experiments, reflections, sigma):

  A = matrix.sqr(experiments[0].crystal.get_A())
  s0 = matrix.col(experiments[0].beam.get_s0())

  s2_obs = flex.vec3_double()
  for i in range(len(reflections)):

    h = matrix.col(reflections[i]['miller_index'])

    r = A * h
    s2 = s0 + r

    s2_obs.append(s2)

  reflections['s2'] = s2_obs
  return reflections


def tst_ideal():

  numpy.random.seed(100)

  # Ensure we have a data block
  experiments = ExperimentListFactory.from_json_file("experiments.json")
  experiments[0].scan.set_oscillation((0, 1.0), deg=True)
  experiments[0].beam.set_s0((0,0,-1))

  s0 = matrix.col(experiments[0].beam.get_s0())

  # The predicted reflections
  reflections = flex.reflection_table.from_predictions_multi(experiments, padding=4)
  print len(reflections)


  sigma = matrix.sqr((1e-6, 0, 0,
                      0, 2e-6, 0,
                      0, 0, 3e-6))

  reflections = generate_observations2(experiments, reflections, sigma)

  s2_list, ctot_list, xbar_list, Sobs_list = generate_from_reflections(s0, sigma, reflections)

  print "Using %d reflections: " % len(s2_list)

  values = flex.double((
    sqrt(1.1e-6),
    0, sqrt(2.1e-6),
    0, 0, sqrt(3.1e-6)))
  offset = flex.double(
    [sqrt(1e-7)  for v in values])


  optimizer = SimpleSimplex(
    values,
    offset,
    Target(
      s0,
      s2_list,
      xbar_list,
      ctot_list,
      Sobs_list,
      test=0), 2000)
  params = optimizer.get_solution()

  M = matrix.sqr((
    params[0], 0, 0,
    params[1], params[2], 0,
    params[3], params[4], params[5]))

  sigma = M*M.transpose()
  print sigma

  expected = matrix.sqr((
    1.00037823148e-06, -5.33165381576e-10, -1.13490868834e-09,
    -5.33165381576e-10, 1.999416147e-06, 2.13056916858e-09,
    -1.13490868834e-09, 2.13056916858e-09, 3.00917159468e-06))

  assert all(1e6*abs(a-b) < 1e-7 for a, b in zip(sigma, expected))

  print 'OK'


def tst_binned():

  numpy.random.seed(100)

  # Ensure we have a data block
  experiments = ExperimentListFactory.from_json_file("experiments.json")
  experiments[0].scan.set_oscillation((0, 1.0), deg=True)
  experiments[0].beam.set_s0((0,0,-1))

  s0 = matrix.col(experiments[0].beam.get_s0())

  # The predicted reflections
  reflections = flex.reflection_table.from_predictions_multi(experiments, padding=4)
  print len(reflections)


  sigma = matrix.sqr((1e-6, 0, 0,
                      0, 2e-6, 0,
                      0, 0, 3e-6))

  reflections = generate_observations2(experiments, reflections, sigma)

  s2_list, ctot_list, xbar_list, Sobs_list = generate_from_reflections_binned(s0, sigma, reflections)

  print "Using %d reflections: " % len(s2_list)

  values = flex.double((
    sqrt(1e-6),
    0, sqrt(2e-6),
    0, 0, sqrt(3e-6)))
  offset = flex.double(
    [sqrt(1e-7)  for v in values])


  optimizer = SimpleSimplex(
    values,
    offset,
    Target(
      s0,
      s2_list,
      xbar_list,
      ctot_list,
      Sobs_list,
      test=0), 2000)
  params = optimizer.get_solution()

  M = matrix.sqr((
    params[0], 0, 0,
    params[1], params[2], 0,
    params[3], params[4], params[5]))

  sigma = M*M.transpose()
  print sigma

  expected = matrix.sqr((
    1.06700290634e-06, -1.338495946e-10, -1.78808654488e-09,
    -1.338495946e-10, 2.09354554783e-06, 1.72284152176e-09,
    -1.78808654488e-09, 1.72284152176e-09, 3.17207510387e-06))

  assert all(1e6*abs(a-b) < 1e-7 for a, b in zip(sigma, expected))

  print 'OK'

if __name__ == '__main__':
  tst_ideal()
  tst_binned()
