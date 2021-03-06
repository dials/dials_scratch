from __future__ import division
from __future__ import print_function

from scitbx.array_family import flex
from dials.model.data import Reflection, ReflectionList
from dials.algorithms import shoebox

rl = ReflectionList()
r1 = Reflection()
r1.shoebox = (10, 20, 10, 20, 10, 20)
r2 = Reflection()
r2.shoebox = (15, 25, 15, 25, 15, 25)
r3 = Reflection()
r3.shoebox = (20, 30, 20, 30, 20, 30)
rl.append(r1)
rl.append(r2)
rl.append(r3)
overlapping = shoebox.find_overlapping(rl)

for e in overlapping.edges():
    print("Edge: ", overlapping.edge_vertices(e))

for v in overlapping.vertices():
    print("Vertex: ", v, " => ", [a for a in overlapping.adjacent_vertices(v)])
