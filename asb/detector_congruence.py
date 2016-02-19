#!/usr/bin/env python
#
# detector_congruence.py
#
#  Copyright (C) 2016 Lawrence Berkeley National Laboratory (LBNL)
#
#  Author: Aaron Brewster
#
#  This code is distributed under the X license, a copy of which is
#  included in the root directory of this package.
#
from __future__ import division
from scitbx.matrix import col
from matplotlib import pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.colors import Normalize
from matplotlib import cm

help_message = '''

This program is used to calculate statisical measurements of consistency
between two detectors

Example:

  libtbx.python asic_stats.py experiment1.json experiment2.json
'''

def iterate_detector_at_level(item, depth, level):
  if level == depth:
    yield item
  else:
    for child in item:
      for subitem in iterate_detector_at_level(child, depth+1, level):
        yield subitem

def iterate_panels(panelgroup):
  if hasattr(panelgroup, 'children'):
    for child in panelgroup:
      for subitem in iterate_panels(child):
        yield subitem
  else:
    yield panelgroup

class Script(object):
  ''' Class to parse the command line options. '''

  def __init__(self):
    ''' Set the expected options. '''
    from dials.util.options import OptionParser
    from libtbx.phil import parse
    import libtbx.load_env
    # Create the phil parameters
    phil_scope = parse('''
    hierarchy_level=0
      .type=int
    ''')

    # Create the option parser
    usage = "usage: %s [options] /path/to/refined/json/file" % libtbx.env.dispatcher_name
    self.parser = OptionParser(
      usage=usage,
      sort_options=True,
      phil=phil_scope,
      read_experiments=True,
      read_datablocks=True,
      epilog=help_message)

  def run(self):
    ''' Parse the options. '''
    from dials.util.options import flatten_experiments, flatten_datablocks
    # Parse the command line arguments
    params, options = self.parser.parse_args(show_diff_phil=True)
    experiments = flatten_experiments(params.input.experiments)
    datablocks = flatten_datablocks(params.input.datablock)

    detectors = []
    detectors.extend(experiments.detectors())
    dbs = []
    for datablock in datablocks:
      dbs.extend(datablock.unique_detectors())
    detectors.extend(dbs)

    # Verify inputs
    if len(detectors) != 2:
      print "Please provide two experiments and or datablocks for comparison"
      return

    """
    print "Testing iterate_detector_at_level"
    for level in xrange(4):
      print "iterating at level", level
      for panelg in iterate_detector_at_level(detectors[0].hierarchy(), 0, level):
        print panelg.get_name()

    print "Testing iterate_panels"
    for level in xrange(4):
      print "iterating at level", level
      for panelg in iterate_detector_at_level(detectors[0].hierarchy(), 0, level):
        for panel in iterate_panels(panelg):
          print panel.get_name()
    """

    angles = {}
    for pg1, pg2 in zip(iterate_detector_at_level(detectors[0].hierarchy(), 0, params.hierarchy_level),
                        iterate_detector_at_level(detectors[1].hierarchy(), 0, params.hierarchy_level)):
      angle = col(pg1.get_normal()).angle(col(pg2.get_normal()), deg=True)
      for p1, p2 in zip(iterate_panels(pg1), iterate_panels(pg2)):
        assert p1.get_name() == p2.get_name()
        angles[p1.get_name()] = angle

    self.detector_plot(detectors[0], angles)

  def detector_plot(self, detector, data):
    norm = Normalize(vmin=min(data.values()), vmax=max(data.values()))
    sm = cm.ScalarMappable(norm=norm, cmap=cm.RdYlGn_r)

    fig = plt.figure()
    ax = fig.add_subplot(111, aspect='equal')
    max_dim = 0
    for panel_id, panel in enumerate(detector):
      size = panel.get_image_size()
      p0 = col(panel.get_pixel_lab_coord((0,0)))
      p1 = col(panel.get_pixel_lab_coord((size[0]-1,0)))
      p2 = col(panel.get_pixel_lab_coord((size[0]-1,size[1]-1)))
      p3 = col(panel.get_pixel_lab_coord((0,size[1]-1)))

      v1 = p1-p0
      v2 = p3-p0
      vcen = ((v2/2) + (v1/2)) + p0

      ax.add_patch(Polygon((p0[0:2],p1[0:2],p2[0:2],p3[0:2]), closed=True, color=sm.to_rgba(data[panel.get_name()]), fill=True))
      ax.annotate("%d %.2f"%(panel_id, data[panel.get_name()]), vcen[0:2], ha='center')

      for p in [p0, p1, p2, p3]:
        for c in p[0:2]:
          if abs(c) > max_dim:
            max_dim = abs(c)
    ax.set_xlim((-max_dim,max_dim))
    ax.set_ylim((-max_dim,max_dim))
    plt.show()


if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
