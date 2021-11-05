import sys

import numpy as np

from matplotlib import pyplot as plt

from dials.array_family import flex

data = [flex.reflection_table.from_file(arg) for arg in sys.argv[1:]]
data = [d.select(d.get_flags(d.flags.scaled)) for d in data]


def cc(r1, r2):
    m12 = r1.match(r2)

    r1 = r1.select(m12[0])
    r2 = r2.select(m12[1])

    i1 = r1["intensity.scale.value"]
    i2 = r2["intensity.scale.value"]

    lc = flex.linear_correlation(i1, i2)
    return lc.coefficient()


l = len(data)
mat = np.empty(shape=(l, l))
for a in range(l):
    for b in range(l):
        mat[a, b] = cc(data[a], data[b])
        # print(a, b, cc(data[a], data[b]))

plt.imshow(mat)
plt.show()
