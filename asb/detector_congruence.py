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
from scitbx.array_family import flex
from scitbx.matrix import col
from matplotlib import pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.colors import Normalize
from matplotlib import cm
import numpy as np
from libtbx.phil import parse

help_message = '''

This program is used to calculate statisical measurements of consistency
between two detectors.

Example:

  libtbx.python asic_stats.py experiment1.json experiment2.json
'''

# Create the phil parameters
phil_scope = parse('''
hierarchy_level=0
  .type=int
  .help=Provide congruence statistics for detector modules at the given hierarchy level.
colormap=RdYlGn_r
  .type=str
  .help=matplotlib color map. See e.g.: \
        http://matplotlib.org/examples/color/colormaps_reference.html
show_plots=True
  .type=bool
  .help=Whether to show congruence plots
''')

def iterate_detector_at_level(item, depth = 0, level = 0):
  """
  Iterate through all panel groups or panels of a detector object at a given
  hierarchy level
  @param item panel group or panel. Use detector.hierarchy().
  @param depth current dept for recursion. Should be 0 for initial call.
  @param level iterate groups at this level
  @return next panel or panel group object
  """
  if level == depth:
    yield item
  else:
    for child in item:
      for subitem in iterate_detector_at_level(child, depth+1, level):
        yield subitem

def iterate_panels(panelgroup):
  """
  Find and iterate all panels in the given panel group, regardless of the hierarchly level
  of this panelgroup
  @param panelgroup the panel group of interest
  @return the next panel
  """
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
    import libtbx.load_env

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
    self.params = params
    experiments = flatten_experiments(params.input.experiments)
    datablocks = flatten_datablocks(params.input.datablock)

    # Find all detector objects
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

    # These lines exercise the iterate_detector_at_level and iterate_panels functions
    # for a detector with 4 hierarchy levels
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

    # Iiterate through the detectors, computing the congruence statistics
    normal_angles = {}
    z_angles = {}
    xy_deltas = {}
    z_deltas = {}
    all_normal_angles = flex.double()
    all_z_angles = flex.double()
    all_xy_deltas = flex.double()
    all_z_deltas = flex.double()
    table_header = ["PanelG","Normal","Z rot","Delta","Delta"]
    table_header2 = ["Id","Angle","Angle","XY","Z"]
    table_data = []
    table_data.append(table_header)
    table_data.append(table_header2)

    for pg_id, (pg1, pg2) in enumerate(zip(iterate_detector_at_level(detectors[0].hierarchy(), 0, params.hierarchy_level),
                                           iterate_detector_at_level(detectors[1].hierarchy(), 0, params.hierarchy_level))):
      norm_angle = col(pg1.get_normal()).angle(col(pg2.get_normal()), deg=True)
      z_angle = col(pg1.get_fast_axis()[0:2]).angle(col(pg2.get_fast_axis()[0:2]), deg=True)
      if hasattr(pg1, 'children'):
        v1 = col(pg1.get_origin())
        v2 = col(pg2.get_origin())
      else:
        s = pg1.get_image_size()
        assert s == pg2.get_image_size()
        v1 = col(pg1.get_pixel_lab_coord((s[0]/2, s[1]/2)))
        v2 = col(pg2.get_pixel_lab_coord((s[0]/2, s[1]/2)))

      delta = v1 - v2
      xyd = col(delta[0:2]).length()
      zd = abs(delta[2])

      all_normal_angles.append(norm_angle)
      all_z_angles.append(z_angle)
      all_xy_deltas.append(xyd)
      all_z_deltas.append(zd)
      table_data.append(["%d"%pg_id, "%.4f"%norm_angle, "%.4f"%z_angle, "%.4f"%xyd, "%.4f"%zd])

      for p1, p2 in zip(iterate_panels(pg1), iterate_panels(pg2)):
        assert p1.get_name() == p2.get_name()
        normal_angles[p1.get_name()] = norm_angle
        z_angles[p1.get_name()] = z_angle
        xy_deltas[p1.get_name()] = xyd
        z_deltas[p1.get_name()] = zd

    table_data.append(["Mean", "%.4f"%flex.mean(all_normal_angles),
                               "%.4f"%flex.mean(all_z_angles),
                               "%.4f"%flex.mean(all_xy_deltas),
                               "%.4f"%flex.mean(all_z_deltas)])
    table_data.append(["Stddev", "%.4f"%flex.mean_and_variance(all_normal_angles).unweighted_sample_standard_deviation(),
                                 "%.4f"%flex.mean_and_variance(all_z_angles).unweighted_sample_standard_deviation(),
                                 "%.4f"%flex.mean_and_variance(all_xy_deltas).unweighted_sample_standard_deviation(),
                                 "%.4f"%flex.mean_and_variance(all_z_deltas).unweighted_sample_standard_deviation()])
    from libtbx import table_utils
    print table_utils.format(table_data,has_header=2,justify='center',delim=" ")

    if params.show_plots:
      # Plot the results
      self.detector_plot(detectors[0], normal_angles, u"Angle between normal vectors (\N{DEGREE SIGN})", u"%.2f\N{DEGREE SIGN}")
      self.detector_plot(detectors[0], z_angles, u"Z rotation angle between panels (\N{DEGREE SIGN})", u"%.2f\N{DEGREE SIGN}")
      self.detector_plot(detectors[0], xy_deltas, u"XY displacements between panels (mm)", u"%.2fmm")
      self.detector_plot(detectors[0], z_deltas, u"Z displacements between panels (mm)", u"%.2fmm")

  def detector_plot(self, detector, data, title, units_str):
    """
    Use matplotlib to plot a detector, color coding panels according to data
    @param detector detector reference detector object
    @param data python dictionary of panel names as keys and numbers as values
    @param title title string for plot
    @param units_str string with a formatting statment for units on each panel
    """
    # initialize the color map
    values = flex.double(data.values())
    norm = Normalize(vmin=flex.min(values), vmax=flex.max(values))
    cmap = plt.cm.get_cmap(self.params.colormap)
    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array(np.arange(flex.min(values), flex.max(values), (flex.max(values)-flex.min(values))/20)) # needed for colorbar

    fig = plt.figure()
    ax = fig.add_subplot(111, aspect='equal')
    max_dim = 0
    for panel_id, panel in enumerate(detector):
      # get panel coordinates
      size = panel.get_image_size()
      p0 = col(panel.get_pixel_lab_coord((0,0)))
      p1 = col(panel.get_pixel_lab_coord((size[0]-1,0)))
      p2 = col(panel.get_pixel_lab_coord((size[0]-1,size[1]-1)))
      p3 = col(panel.get_pixel_lab_coord((0,size[1]-1)))

      v1 = p1-p0
      v2 = p3-p0
      vcen = ((v2/2) + (v1/2)) + p0

     # add the panel to the plot
      ax.add_patch(Polygon((p0[0:2],p1[0:2],p2[0:2],p3[0:2]), closed=True, color=sm.to_rgba(data[panel.get_name()]), fill=True))
      ax.annotate("%d %s"%(panel_id, units_str%data[panel.get_name()]), vcen[0:2], ha='center')

      # find the plot maximum dimensions
      for p in [p0, p1, p2, p3]:
        for c in p[0:2]:
          if abs(c) > max_dim:
            max_dim = abs(c)

    # plot the results
    ax.set_xlim((-max_dim,max_dim))
    ax.set_ylim((-max_dim,max_dim))
    ax.set_xlabel("mm")
    ax.set_ylabel("mm")
    fig.colorbar(sm)
    plt.title(title)
    plt.show()

if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)