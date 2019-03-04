from __future__ import print_function
def allowed(perm):
  m = perm[0] + perm[1]
  if m != (perm[-2] + perm[-1]): return None
  for j in range(1, len(perm) - 1, 2):
    if sum(perm[j:j+3]) != m: return None
  return perm

def solver(n):
  import itertools
  perms = 0
  solns = []
  for perm in itertools.permutations(range(1, 2 * n)):
    perms += 1
    if allowed(perm):
      if not tuple(reversed(perm)) in solns: solns.append(perm)
  return perms, solns

if __name__ == '__main__':
  import sys
  perms, solns = solver(int(sys.argv[1]))
  print('Of %d permutations %d unique solutions:' % (perms, len(solns)))
  for s in solns: print(s)
