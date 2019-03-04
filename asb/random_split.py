from __future__ import division
import sys, os, random
from libtbx import easy_pickle
from dials.array_family import flex
from libtbx.phil import parse
from libtbx.utils import Sorry

"""
Jiffy script to split a reflection table into two random halves.
Usage: libtbx.python random_split.py reflections.pickle
"""

phil_scope = parse("""
  output_selection = None
    .type = str
    .help = If not None, save the selection in the path indicated for future use
  use_selection = None
    .type = str
    .help = If not None, use the random split selection from a previous run
""")

def run(args):
  user_phil = []
  for arg in args:
    if os.path.isfile(arg):
      filename = arg
    else:
      try:
        user_phil.append(parse(arg))
      except Exception as e:
        raise Sorry("Unrecognized argument: %s"%arg)

  params = phil_scope.fetch(sources = user_phil).extract()

  name = os.path.basename(filename)
  base, ext = os.path.splitext(name)
  filea = base + "_a" + ext
  fileb = base + "_b" + ext

  data = easy_pickle.load(filename)

  if params.use_selection is None:
    sel = flex.random_permutation(len(data))
  else:
    sel = easy_pickle.load(params.use_selection)
    assert len(sel) == len(data), "Length of selection doesn't match length of input"

  data_a = data.select(sel[:len(data)//2])
  data_b = data.select(sel[len(data)//2:])
  data_a = data_a.select(flex.sort_permutation(data_a['id']))
  data_b = data_b.select(flex.sort_permutation(data_b['id']))

  assert len(data_a) + len(data_b) == len(data)

  easy_pickle.dump(filea, data_a)
  easy_pickle.dump(fileb, data_b)

  if params.output_selection is not None:
    easy_pickle.dump(params.output_selection, sel)

if __name__ == "__main__":
  run(sys.argv[1:])
