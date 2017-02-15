from dials.array_family import flex
import cPickle as pickle
import sys

reflection_file = sys.argv[1]
data = pickle.load(open(reflection_file, 'rb'))
print '%d reflections' % data.size()

i = data['intensity.sum.value']
v = data['intensity.sum.variance']
s = flex.sqrt(v)
i_s = i / s

h = flex.histogram(i_s, -5, 20, n_slots=50)
c = h.slot_centers()
d = h.slots()

for _c, _d in zip(c, d):
  print _c, _d
