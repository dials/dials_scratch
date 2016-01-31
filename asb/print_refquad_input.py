from __future__ import division
import os, sys, random

paths = []
n_subset = None
for arg in sys.argv[1:]:
  if os.path.isdir(arg):
    paths.append(arg)
  else:
    try:
      n_subset = int(arg)
    except Exception, e:
      pass

all_exp = []
all_ref = []

for path in paths:
  for filename in os.listdir(path):
    if "indexed" in filename:
      exp_path = os.path.join(path, filename.rstrip("_indexed.pickle") + "_refined_experiments.json")
      if not os.path.exists(exp_path): continue
      all_exp.append(exp_path)
      all_ref.append(os.path.join(path, filename))

if n_subset is not None:
  subset_all_exp = []
  subset_all_ref = []
  n_picked = 0

  while n_picked < n_subset:
    idx = random.randint(0, len(all_exp)-1)
    subset_all_exp.append(all_exp.pop(idx))
    subset_all_ref.append(all_ref.pop(idx))
    n_picked += 1
  all_exp = subset_all_exp
  all_ref = subset_all_ref

for exp_path, ref_path in zip(all_exp, all_ref):
  print "input {"
  print "  experiments =", exp_path
  print "  reflections =", ref_path
  print "}"

