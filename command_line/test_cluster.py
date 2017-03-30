
def test_method(method, n):
  from dials.util.mp import parallel_map
  from dials_scratch.jmp.cluster_func_test import func
  from time import time
  result = parallel_map(
    func      = func,
    iterable  = list(range(n)),
    processes = n,
    nslots    = 8,
    method    = method)

def test_method_n_times(method, m, n):
  from time import time
  st = time()
  for i in range(n):
    test_method(method, m)
  t = time() - st
  return t / n

if __name__ == '__main__':
  for m in [10, 50, 100, 500, 1000]:
    print "EASY_MP: %d jobs -> %f seconds" % (m, test_method_n_times("sge", m, 1))
  for m in [10, 50, 100, 500, 1000]:
    print "DRMAA:   %d jobs -> %f seconds" % (m, test_method_n_times("drmaa", m, 1))
