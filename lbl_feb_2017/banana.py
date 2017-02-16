from dials.array_family import flex
import cPickle as pickle
import sys

reflection_file = sys.argv[1]
data = pickle.load(open(reflection_file, 'rb'))

i = data['intensity.sum.value']
v = data['intensity.sum.variance']
s = flex.sqrt(v)
i_s = i / s
dq = data['dq']

def histogram():
  h = flex.histogram(i_s, -5, 20, n_slots=50)
  c = h.slot_centers()
  d = h.slots()

  for _c, _d in zip(c, d):
    print _c, _d

def scatter():
  for j in range(i_s.size()):
    print i_s[j], dq[j]

scatter()
