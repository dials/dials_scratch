from __future__ import division
try:
  import boost.python
except Exception:
  ext = None
else:
  ext = boost.python.import_ext("dials_scratch_ext", optional = False)

if not ext is None:
  from dials_scratch_ext import *
