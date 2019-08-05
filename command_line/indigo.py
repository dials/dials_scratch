# LIBTBX_SET_DISPATCHER_NAME dev.dials.indigo
from __future__ import absolute_import, division

from dials_scratch.dgw.indigo import run

if __name__ == "__main__":
    import sys

    run(sys.argv[1:])
