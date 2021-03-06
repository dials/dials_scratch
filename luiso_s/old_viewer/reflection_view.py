from __future__ import absolute_import

#
#  DIALS viewer
#
#  Copyright (C) 2014 Diamond Light Source
#
#  Author: Luis Fuentes-Montero (Luiso)
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package."

# should be called with:
# dials.python new_viewer_tst.py ../../../../../../example/indexed.pickle
# or dials.python new_viewer_tst.py PATH/TO/MY/pickle/file
#
import wx
from .viewer_frame import ReflectionFrame


class viewer_App(wx.App):
    def OnInit(self):
        self.frame = ReflectionFrame(
            None, -1, "DIALS Reflections Viewer", wx.DefaultPosition, (550, 200)
        )

        self.SetTopWindow(self.frame)
        self.frame.Show(True)
        return True

    def table_in(self, loc_tabl):
        self.frame.tabl_to_frame(loc_tabl)


if __name__ == "__main__":

    import cPickle as pickle
    import sys
    from dials.array_family import flex

    table = flex.reflection_table.from_pickle(sys.argv[1])

    pos_of_best_example_so_far = """
  .... dials_regression/refinement_test_data/radiation_damaged_thaumatin/
  indexed.pickle
  """
    # print "num of ref =", len(table)

    app = viewer_App()
    app.table_in(table)

    app.MainLoop()
