from scitbx.array_family import flex
from scitbx.random import variate, uniform_distribution, poisson_distribution
import math

rate = 100
nn = 1000000
ss = 3

scale = variate(uniform_distribution(min=-ss, max=ss))
intensity = variate(poisson_distribution(mean=rate))

d = flex.double(nn)

for j in range(nn):
  x = next(scale)
  d[j] = math.exp(- x * x) * next(intensity)

h = flex.histogram(d, data_min=0, data_max=2*rate, n_slots=100)

total = 0
for c, s in zip(h.slot_centers(), h.slots()):
  total += s
  print c, s, total / nn
