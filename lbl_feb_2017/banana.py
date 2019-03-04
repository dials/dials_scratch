from __future__ import print_function
from dials.array_family import flex
import cPickle as pickle
import sys

reflection_file = sys.argv[1]
data = pickle.load(open(reflection_file, 'rb'))

i = data['intensity.sum.value']
v = data['intensity.sum.variance']
obs = data['xyzobs.px']
cal = data['xyzcal.px']
s = flex.sqrt(v)
i_s = i / s
dq = data['dq']

def histogram():
  h = flex.histogram(i_s, -5, 20, n_slots=50)
  c = h.slot_centers()
  d = h.slots()

  for _c, _d in zip(c, d):
    print(_c, _d)

def scatter():
  for j in range(i_s.size()):
    print(dq[j], i[j], i_s[j], cal[j][0], cal[j][1], obs[j][0], obs[j][1])

scatter()
