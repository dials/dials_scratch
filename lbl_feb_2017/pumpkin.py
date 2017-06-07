from dials.array_family import flex
import cPickle as pickle
import math
import sys
from scitbx import matrix
import scitbx.math as smath
import json

filename = sys.argv[1]

data = pickle.load(open(filename, 'rb'))
print '%d reflections' % data.size()

# keep = data.select(data.get_flags(data.flags.integrated))
# print keep.size()

data = data.select(data['intensity.sum.variance'] > 0)

delpsi = data['delpsical.rad'] * 180.0 / math.pi
i = data['intensity.sum.value']
v = data['intensity.sum.variance']
s = flex.sqrt(v)
i_s = i/s

from dxtbx.model.experiment_list import ExperimentListFactory
expt = ExperimentListFactory.from_json_file(sys.argv[2])
crystal = expt.crystals()[0]
UB = matrix.sqr(crystal.get_A())
s0 = matrix.col(expt.beams()[0].get_s0())

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

def plotify(params):
  rx, ry, rz = params
  R = matrix.sqr(smath.euler_angles_as_matrix((rx, ry, rz), deg=True))
  score = 0.0
  for j in range(data.size()):
    hkl = data['miller_index'][j]
    q = R * UB * hkl
    print (q + s0).length() - matrix.col(data['s1'][j]).length(), i_s[j]

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
  def __init__(self, values, offset, max_iter):
    self.n = len(values)
    self.x = values
    self.starting_simplex = generate_start(values, offset)
    self.fcount = 0

    optimizer = simplex.simplex_opt(dimension=self.n,
                                    matrix=self.starting_simplex,
                                    evaluator=self,
                                    tolerance=1e-10,
                                    max_iter=max_iter)

    self.x = optimizer.get_solution()
    return

  def get_solution(self):
    return self.x

  def target(self, vector):
    score = scorify(vector)
    return score

doohicky = simple_simplex(
  flex.double((0, 0, 0)), flex.double((0.1, 0.1, 0.1)), 2000)
best = doohicky.get_solution()
plotify(best)
