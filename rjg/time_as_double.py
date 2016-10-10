import time
from dials.array_family import flex
n = 100
a = flex.int(flex.grid(2527, 2463))
t0 = time.time()
for i in range(n):
  b = a.as_double()
  #c = b.iround()
  #c = flex.double(flex.grid(2527, 2463))
  #c = flex.random_double(a.size())
t1 = time.time()
t = t1-t0
print "%.2fs (%.3fs/call)" %(t, t/n)

import numpy as np
a = np.random.randint(100, size=a.all())
t0 = time.time()
for i in range(n):
  b = a.astype(np.float64)
t1 = time.time()
t = t1-t0
print "%.2fs (%.3fs/call)" %(t, t/n)
