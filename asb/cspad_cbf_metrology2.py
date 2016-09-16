from __future__ import division
import os
import libtbx.load_env
from libtbx import easy_run

data_path = "/net/dials/raid1/aaron/cxim1416/strong_006"
indexing_phil = "/net/dials/raid1/aaron/cxim1416/indexing/index.phil"
refinement_phil = "/net/dials/raid1/aaron/cxim1416/metrology/refine.phil"

indexing_program = os.path.join(libtbx.env.find_in_repositories("dials_scratch"), "asb", "mp_index.py")

steps = {}
steps[0] = [2, 3]
steps[1] = steps[0] + [0, 1]
steps[2] = steps[1] + [14, 15]
steps[3] = steps[2] + [6, 7]
steps[4] = steps[3] + [4, 5]
steps[5] = steps[4] + [12, 13]
steps[6] = steps[5] + [8, 9]
steps[7] = steps[6] + [10, 11]

for s, panels in steps.iteritems():
  rest = []
  for p in panels:
    rest.append(p+16)
    rest.append(p+32)
    rest.append(p+48)
  panels.extend(rest)

levels = {0: 1}
for i in xrange(7):
  levels[i+1] = 2

cwd = os.getcwd()

for d in ["indexing", "metrology"]:
  if not os.path.exists(d):
    os.makedirs(d)

for i in xrange(8):
  os.chdir("indexing")

  trial = "%03d"%i
  if not os.path.exists(trial):
    os.makedirs(trial)
  command = "libtbx.python %s %s mp.nproc=64 %s output_dir=%s"%(indexing_program, data_path, indexing_phil, trial)
  if i > 0:
    prev_level = levels[i-1]
    prev_trial = "%03d"%(i-1)
    geom_path = os.path.join(cwd, "metrology", prev_trial, "t%s_refined_experiments_level%d.json"%(prev_trial, prev_level))
    command += " reference_geometry=%s"%geom_path

  print command
  easy_run.fully_buffered(command).show_stdout()

  os.chdir("../metrology")
  if not os.path.exists(trial):
    os.makedirs(trial)
  os.chdir(trial)

  command = "cspad.cbf_metrology %s tag=t%s refine_to_hierarchy_level=%d rmsd_filter.enable=False panel_filter=%s %s"%(
    os.path.join(cwd, "indexing", trial), trial, levels[i], ",".join(["%d"%p for p in steps[i]]), refinement_phil)

  print command
  easy_run.fully_buffered(command).show_stdout()

  os.chdir(cwd)
