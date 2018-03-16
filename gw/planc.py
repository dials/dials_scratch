from __future__ import division, print_function
from dials.array_family import flex
import cPickle as pickle

data = pickle.load(open('strong.pickle'))
x, y, z = data['xyzobs.px.value'].parts()
xy = flex.vec2_double(x, y)

zs = range(int(flex.min(z)), int(flex.max(z)) + 1)
n = len(zs)

sij = flex.double(n * n, 0.0)
sij.reshape(flex.grid(n, n))

for z0 in zs:
  s0 = z == (z0 + 0.5)
  xy0 = xy.select(s0)
  n0 = xy0.size()
  if n0 < 5:
    continue
  from annlib_ext import AnnAdaptor as ann_adaptor
  ann = ann_adaptor(xy0.as_double().as_1d(), 2)
  for z1 in zs:
    if z1 >= z0:
      break
    s1 = z == (z1 + 0.5)
    xy1 = xy.select(s1)
    n1 = xy1.size()
    if n1 < 5:
      continue
    ann.query(xy1.as_double().as_1d())
    d1 = flex.sqrt(ann.distances)
    m01 = (d1 < 5.0).count(True)
    s = m01 / (0.5 * (n0 + n1))
    sij[z0, z1] = s
    sij[z1, z0] = s

pickle.dump(sij, open('sij.pickle', 'w'))

hij = flex.histogram(sij.as_1d(), data_min=0, data_max=1, n_slots=100)
for _c, _s in zip(hij.slot_centers(), hij.slots()):
  print(_c, _s)
