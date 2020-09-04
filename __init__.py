from __future__ import division

try:
    import boost_adaptbx.boost.python
except Exception:
    ext = None
else:
    ext = boost_adaptbx.boost.python.import_ext("dials_scratch_ext", optional=True)

if not ext is None:
    from dials_scratch_ext import *
