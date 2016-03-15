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
tag = None
  .type = str
  .help = Used in the plot titles
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

def id_from_name(detector, name):
  """ Jiffy function to get the id of a panel using its name
  @param detector detector object
  @param name panel name
  @return index of panel in detector
  """
  return [p.get_name() for p in detector].index(name)

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
      read_reflections=True,
      epilog=help_message)

  def run(self):
    ''' Parse the options. '''
    from dials.util.options import flatten_experiments, flatten_datablocks, flatten_reflections
    # Parse the command line arguments
    params, options = self.parser.parse_args(show_diff_phil=True)
    self.params = params
    experiments = flatten_experiments(params.input.experiments)
    datablocks = flatten_datablocks(params.input.datablock)
    reflections = flatten_reflections(params.input.reflections)

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

    # Iterate through the detectors, computing the congruence statistics
    delta_normals = {}
    z_angles = {}
    xy_deltas = {}
    z_deltas = {}
    refl_counts = {}
    all_delta_normals = flex.double()
    all_z_angles = flex.double()
    all_xy_deltas = flex.double()
    all_z_deltas = flex.double()
    all_refls_count = flex.int()

    all_normal_angles = flex.double()
    all_rot_z = flex.double()
    pg_bc_dists = flex.double()
    all_bc_dist = flex.double()
    all_z_offsets = flex.double()
    all_weights = flex.double()

    precision_table_data = []
    detector_table_data = []
    root1 = detectors[0].hierarchy()
    root2 = detectors[1].hierarchy()

    root_normal_1 = col(root1.get_normal())
    root_normal_2 = col(root2.get_normal())

    s0 = col(flex.vec3_double([col(b.get_s0()) for b in experiments.beams()]).mean())

    for pg_id, (pg1, pg2) in enumerate(zip(iterate_detector_at_level(root1, 0, params.hierarchy_level),
                                           iterate_detector_at_level(root2, 0, params.hierarchy_level))):
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
      xyd = col(delta[0:2]).length()*1000
      zd = abs(delta[2])*1000

      all_delta_normals.append(norm_angle)
      all_z_angles.append(z_angle)
      all_xy_deltas.append(xyd)
      all_z_deltas.append(zd)

      total_refls = 0
      total_refls_1 = 0
      total_refls_2 = 0
      for p1, p2 in zip(iterate_panels(pg1), iterate_panels(pg2)):
        assert p1.get_name() == p2.get_name()
        z_angles[p1.get_name()] = z_angle
        xy_deltas[p1.get_name()] = xyd
        z_deltas[p1.get_name()] = zd
        r1 = len(reflections[0].select(reflections[0]['panel'] == id_from_name(detectors[0], p1.get_name())))
        r2 = len(reflections[1].select(reflections[1]['panel'] == id_from_name(detectors[1], p2.get_name())))
        total_refls += r1 + r2
        total_refls_1 += r1
        total_refls_2 += r2
        refl_counts[p1.get_name()] = r1 + r2

      all_refls_count.append(total_refls)

      # Compute distances between panel groups and beam center
      # Also compute offset along Z axis
      dists = flex.double()
      z_offsets = flex.double()
      for pg, r in zip([pg1, pg2], [root1, root2]):
        bc = col(pg.get_beam_centre_lab(s0))
        if hasattr(pg, 'children'):
          ori = col(pg.get_origin())
        else:
          x, y = pg.get_image_size_mm()
          offset = col((x, y))/2
          ori = col(pg.get_lab_coord(offset))

        dists.append((ori-bc).length())

        rori = col(r.get_origin())
        delta_ori = ori-rori
        r_norm = col(r.get_normal())
        z_offsets.append(r_norm.dot(delta_ori))

      stats = flex.mean_and_variance(dists, flex.double([total_refls_1, total_refls_2]))
      dist_m = stats.mean()
      dist_s = stats.gsl_stats_wsd()
      pg_bc_dists.append(dist_m)
      all_bc_dist.extend(dists)
      all_z_offsets.extend(z_offsets)
      precision_table_data.append(["%d"%pg_id, "%5.1f"%dist_m, "%.4f"%dist_s, "%.4f"%norm_angle, "%.4f"%z_angle, "%4.1f"%xyd, "%4.1f"%zd, "%6d"%total_refls])

      # Compute angle between detector normal and panel group normal
      norm_angle_1 = root_normal_1.angle(col(pg1.get_normal()), deg=True)
      norm_angle_2 = root_normal_2.angle(col(pg2.get_normal()), deg=True)
      stats = flex.mean_and_variance(flex.double([norm_angle_1, norm_angle_2]), flex.double([total_refls_1, total_refls_2]))
      norm_angle_m = stats.mean()
      norm_angle_s = stats.gsl_stats_wsd()
      all_normal_angles.append(norm_angle_1)
      all_normal_angles.append(norm_angle_2)
      all_weights.append(total_refls_1)
      all_weights.append(total_refls_2)

      # Compute rotation of panel group around detector normal
      pg_rotz = flex.double()
      for pg, r in zip([pg1, pg2], [root1, root2]):
        pgf = col(pg.get_fast_axis())
        rf = col(r.get_fast_axis())
        rs = col(r.get_slow_axis())

        # v is the component of pgf in the rf rs plane
        v = (rf.dot(pgf) * rf) + (rs.dot(pgf) * rs)
        pg_rotz.append(rf.angle(v, deg=True))
        all_rot_z.append(pg_rotz[-1])

      stats = flex.mean_and_variance(pg_rotz, flex.double([total_refls_1, total_refls_2]))
      rotz_m = stats.mean()
      rotz_s = stats.gsl_stats_wsd()

      stats = flex.mean_and_variance(z_offsets, flex.double([total_refls_1, total_refls_2]))
      zo_m = stats.mean()
      zo_s = stats.gsl_stats_wsd()

      detector_table_data.append(["%d"%pg_id, "%5.1f"%dist_m, "%.4f"%dist_s, "%.4f"%norm_angle_m, "%.4f"%norm_angle_s, "%6.2f"%rotz_m, "%.4f"%rotz_s, "%.4f"%zo_m, "%.4f"%zo_s, "%6d"%total_refls])

    table_d = {d:row for d, row in zip(pg_bc_dists, precision_table_data)}
    table_header = ["PanelG","Dist","Dist","Normal","Z rot","Delta","Delta","N"]
    table_header2 = ["Id","","Sigma","Angle","Angle","XY","Z","Refls"]
    table_header3 = ["", "(mm)","(mm)","(deg)","(deg)","(microns)","(microns)",""]
    precision_table_data = [table_header, table_header2, table_header3]
    precision_table_data.extend([table_d[key] for key in sorted(table_d)])

    table_d = {d:row for d, row in zip(pg_bc_dists, detector_table_data)}
    table_header = ["PanelG","Dist","Dist","Normal","Normal","RotZ", "RotZ", "Z Offset", "Z Offset", "N"]
    table_header2 = ["Id","","Sigma","Angle","Angle S","", "Sigma", "", "Sigma", "Refls"]
    table_header3 = ["", "(mm)","(mm)","(deg)","(deg)","(deg)","(deg)","(mm)","(mm)",""]
    detector_table_data = [table_header, table_header2, table_header3]
    detector_table_data.extend([table_d[key] for key in sorted(table_d)])

    r1 = ["Weighted mean"]
    r2 = ["Weighted stddev"]
    r1.append("")
    r2.append("")
    r1.append("")
    r2.append("")
    stats = flex.mean_and_variance(all_delta_normals, all_refls_count.as_double())
    r1.append("%.4f"%stats.mean())
    r2.append("%.4f"%stats.gsl_stats_wsd())
    stats = flex.mean_and_variance(all_z_angles, all_refls_count.as_double())
    r1.append("%.4f"%stats.mean())
    r2.append("%.4f"%stats.gsl_stats_wsd())
    stats = flex.mean_and_variance(all_xy_deltas, all_refls_count.as_double())
    r1.append("%4.1f"%stats.mean())
    r2.append("%4.1f"%stats.gsl_stats_wsd())
    stats = flex.mean_and_variance(all_z_deltas, all_refls_count.as_double())
    r1.append("%4.1f"%stats.mean())
    r2.append("%4.1f"%stats.gsl_stats_wsd())
    r1.append("")
    r2.append("")
    precision_table_data.append(r1)
    precision_table_data.append(r2)
    precision_table_data.append(["Mean", "", "", "", "", "", "", "%6.1f"%flex.mean(all_refls_count.as_double())])

    from libtbx import table_utils
    print "Congruence statistics:"
    print table_utils.format(precision_table_data,has_header=3,justify='center',delim=" ")

    print "PanelG Id: panel group id or panel id, depending on hierarchy_level. For each panel group, statistics are computed between the matching panel groups between the two input experiments."
    print "Dist: distance from center of panel group to the beam center"
    print "Dist Sigma: weighted standard deviation of the measurements used to compute Dist"
    print "Normal angle: angle between the normal vectors of matching panel groups."
    print "Z rot: angle between the XY components of the fast axes of the panel groups."
    print "Delta XY: XY shift between matching panel groups."
    print "Delta Z: Z shift between matching panel groups."
    print "N refls: number of reflections summed between both matching panel groups. This number is used as a weight when computing means and standard deviations."
    print
    print


    stats1 = flex.mean_and_variance(all_normal_angles, all_weights.as_double())
    stats2 = flex.mean_and_variance(all_z_offsets, all_weights.as_double())
    detector_table_data.append(["All", "", "", "%.4f"%stats1.mean(), "%.4f"%stats1.gsl_stats_wsd(), "", "", "%.4f"%stats2.mean(), "%.4f"%stats2.gsl_stats_wsd(), ""])
    detector_table_data.append(["Mean", "", "", "", "", "", "", "", "", "%6.1f"%flex.mean(all_weights.as_double())])

    print "Detector level statistics:"
    print table_utils.format(detector_table_data,has_header=3,justify='center',delim=" ")

    print "PanelG Id: panel group id or panel id, depending on hierarchy_level. For each panel group, statistics are computed using the matching panel groups between the two input experiments."
    print "Dist: distance from center of panel group to the beam center"
    print "Dist Sigma: weighted standard deviation of the measurements used to compute Dist"
    print "Normal Angle: angle between the normal vector of the detector at its root hierarchy level and the normal of the panel group"
    print "Normal Angle S: weighted standard deviation of the measurements used to compute Normal Angle"
    print "RotZ: rotation of each panel group around the detector normal"
    print "RotZ Sigma: weighted standard deviation of the measurements used to compute RotZ"
    print "Z Offset: offset of panel group along the detector normal"
    print "Z Offset Sigma: weighted standard deviation of the measurements used to compute Z Offset"
    print "N refls: number of reflections summed between both matching panel groups. This number is used as a weight when computing means and standard deviations."

    for d_id, d in enumerate(detectors):
      ori = d.hierarchy().get_origin()
      norm = d.hierarchy().get_normal()
      print "Detector", d_id, "origin:", ori
      print "Detector", d_id, "normal:", norm

    if params.tag is None:
      tag = ""
    else:
      tag = "%s "%params.tag

    if params.show_plots:
      # Plot the results
      self.detector_plot_dict(detectors[0], refl_counts, u"%sN reflections"%tag, u"%6d", show=False)
      self.detector_plot_dict(detectors[0], delta_normals, u"%sAngle between normal vectors (\N{DEGREE SIGN})"%tag, u"%.2f\N{DEGREE SIGN}", show=False)
      self.detector_plot_dict(detectors[0], z_angles, u"%sZ rotation angle between panels (\N{DEGREE SIGN})"%tag, u"%.2f\N{DEGREE SIGN}", show=False)
      self.detector_plot_dict(detectors[0], xy_deltas, u"%sXY displacements between panels (microns)"%tag, u"%4.1f", show=False)
      self.detector_plot_dict(detectors[0], z_deltas, u"%sZ displacements between panels (microns)"%tag, u"%4.1f", show=False)
      plt.show()

  def detector_plot_dict(self, detector, data, title, units_str, show=True, reverse_colormap=False):
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
    if reverse_colormap:
      cmap = plt.cm.get_cmap(self.params.colormap + "_r")
    else:
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
    if show:
      plt.show()

if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
