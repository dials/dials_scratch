from dials.array_family import flex
import cPickle as pickle
import math
import sys
from scitbx import matrix
import scitbx.math as smath
import json

def scorify(params):
  rx, ry, rz = params
  R = matrix.sqr(smath.euler_angles_as_matrix((rx, ry, rz), deg=True))
  score = 0.0
  for j in range(data.size()):
    if i_s[j] < 3:
      continue
    hkl = data['miller_index'][j]
    q = R * UB * hkl
    score += i_s[j] * abs((q + s0).length() - matrix.col(data['s1'][j]).length())
  return score

from scitbx import simplex

def generate_start(values, offset):
  assert len(values) == len(offset)
  start = [values]
  for j, o in enumerate(offset):
    next = values.deep_copy()
    next[j] += o
    start.append(next)
  return start

class simple_simplex(object):
  def __init__(self, values, offset, evaluator, max_iter):
    self.n = len(values)
    self.x = values
    self.starting_simplex = generate_start(values, offset)
    self.fcount = 0

    optimizer = simplex.simplex_opt(dimension=self.n,
                                    matrix=self.starting_simplex,
                                    evaluator=evaluator,
                                    tolerance=1e-10,
                                    max_iter=max_iter)

    self.x = optimizer.get_solution()
    return

  def get_solution(self):
    return self.x

  def target(self, vector):
    score = scorify(vector)
    return score

class Blueberry(object):

  def __init__(self, reflection_file, experiment_file):
    # make a pie

    data = pickle.load(open(reflection_file, 'rb'))
    print '%d reflections' % data.size()

    self.data = data.select(data['intensity.sum.variance'] > 0)

    i = data['intensity.sum.value']
    v = data['intensity.sum.variance']
    s = flex.sqrt(v)
    self.i_s = i/s

    from dxtbx.model.experiment_list import ExperimentListFactory
    expt = ExperimentListFactory.from_json_file(experiment_file)
    crystal = expt.crystals()[0]
    self.s0 = matrix.col(expt.beams()[0].get_s0())

    from dials.algorithms.refinement.parameterisation.crystal_parameters \
      import CrystalUnitCellParameterisation, \
      CrystalOrientationParameterisation

    self.crystal = crystal
    self.cucp = CrystalUnitCellParameterisation(crystal)
    self.cop = CrystalOrientationParameterisation(crystal)

    # 0-point and deltas
    values = flex.double(self.cucp.get_param_vals() +
                         self.cop.get_param_vals())
    offset = flex.double([0.01 * v for v in self.cucp.get_param_vals()] +
                         [0.1, 0.1, 0.1])

    initial = crystal.get_unit_cell()
    initial_score = self.target(values)
    doohicky = simple_simplex(values, offset, self, 2000)
    best = doohicky.get_solution()
    print 'Initial cell:', initial
    print 'Final cell:  ', crystal.get_unit_cell()
    print 'Score change', initial_score, self.target(best, do_print=False)
    self.best = best

  def target(self, vector, do_print=False):
    cell_parms = self.cucp.get_param_vals()
    orientation_parms = self.cop.get_param_vals()
    assert len(vector) == len(cell_parms) + len(orientation_parms)
    tst_cell = vector[:len(cell_parms)]
    tst_orientation = vector[len(cell_parms):len(cell_parms) +
                             len(orientation_parms)]

    self.cucp.set_param_vals(tst_cell)
    self.cop.set_param_vals(tst_orientation)

    from scitbx import matrix

    if do_print:
      print 'Cell: %.3f %.3f %.3f %.3f %.3f %.3f' % \
        tuple(self.crystal.get_unit_cell().parameters())
      print 'Phi(1,2,3): %.3f %.3f %.3f' % tuple(tst_orientation)

    UB = matrix.sqr(self.crystal.get_A())

    return self.score(UB)

  def score_old(self, RUB):
    score = 0.0
    for j in range(self.data.size()):
      hkl = self.data['miller_index'][j]
      q = RUB * hkl
      score += self.i_s[j] * abs(
        (q + self.s0).length() - matrix.col(self.data['s1'][j]).length())
    return score

  def score(self, RUB):
    score = 0.0
    for j in range(self.data.size()):
      hkl = self.data['miller_index'][j]
      q = RUB * hkl
      score += abs((q + self.s0).length() - self.s0.length())
    return score

  def plotify(self):
    vector = self.best
    cell_parms = self.cucp.get_param_vals()
    orientation_parms = self.cop.get_param_vals()
    assert len(vector) == len(cell_parms) + len(orientation_parms)
    tst_cell = vector[:len(cell_parms)]
    tst_orientation = vector[len(cell_parms):len(cell_parms) +
                             len(orientation_parms)]

    self.cucp.set_param_vals(tst_cell)
    self.cop.set_param_vals(tst_orientation)

    from scitbx import matrix

    UB = matrix.sqr(self.crystal.get_A())
    data = self.data
    for j in range(data.size()):
      hkl = data['miller_index'][j]
      q = UB * hkl
      print (q + self.s0).length() - matrix.col(data['s1'][j]).length(), self.i_s[j]

    return

blueberry = Blueberry(sys.argv[1], sys.argv[2])
blueberry.plotify()
