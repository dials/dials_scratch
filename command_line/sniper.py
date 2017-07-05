def allowed(perm):
  m = perm[0] + perm[1]
  if m != (perm[-2] + perm[-1]):
    return None
  for j in range(1, len(perm) - 1, 2):
    if sum(perm[j:j+3]) != m:
      return None
  return perm

def csolver(n):
  import itertools
  permutations = 0
  solutions = []
  for perm in itertools.permutations(range(1, 2 * n)):
    permutations += 1
    if allowed(perm):
      if not tuple(reversed(perm)) in solutions:
        solutions.append(perm)
  print 'Of %d permutations %d solutions:' % (permutations, len(solutions))
  for s in solutions: print s

if __name__ == '__main__':
  import sys, time
  t0 = time.time()
  csolver(int(sys.argv[1]))
  t1 = time.time()
  print 'Took %.2fs' % (t1 - t0)
