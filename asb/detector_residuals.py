# -*- coding: utf-8 -*-
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
from dials.array_family import flex
from scitbx.matrix import col
from matplotlib import pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.colors import Normalize
from matplotlib import cm
import numpy as np
from libtbx.phil import parse
import math

help_message = '''

This program is used to calculate statisical measurements of consistency
between observed and predicted reflections

Example:

  libtbx.python detector_residuals.py experiment.json reflections.pickle
'''

# Create the phil parameters
phil_scope = parse('''
dot_size = 10
  .type = int
  .help = Size of dots in detector plots
panel_numbers = True
  .type = bool
  .help = Whether to show panel numbers on each panel
residuals {
  plot_max=None
    .type = float
    .help = Maximum residual value to be shown in the detector plot
  histogram_max=None
    .type = float
    .help = Maximum x value to be used computing the histogram
  histogram_xmax=None
    .type = float
    .help = Maximum x value to be shown in the histogram
  histogram_ymax=None
    .type = float
    .help = Maximum y value to be shown in the histogram
}
include scope dials_scratch.asb.detector_congruence.phil_scope
''', process_includes=True)

from detector_congruence import iterate_detector_at_level, iterate_panels, id_from_name
from detector_congruence import Script as DCScript
class Script(DCScript):
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

  def get_normalized_colors(self, data, vmin=None, vmax=None):
    if vmax is None:
      vmax = self.params.residuals.plot_max
    if vmax is None:
      vmax = flex.max(data)
    if vmin is None:
      vmin = flex.min(data)

    # initialize the color map
    norm = Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.cm.get_cmap(self.params.colormap)
    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    color_vals = np.linspace(vmin, vmax, 11)
    sm.set_array(color_vals) # needed for colorbar

    return norm, cmap, color_vals, sm

  def plot_deltas(self, reflections, panel = None, ax = None, bounds = None):
    assert panel is not None and ax is not None and bounds is not None

    data = (reflections['xyzcal.mm']-reflections['xyzobs.mm.value']).norms()
    norm, cmap, color_vals, sm = self.get_normalized_colors(data)
    deltas = (reflections['xyzcal.mm']-reflections['xyzobs.mm.value'])*self.delta_scalar

    x, y = panel.get_image_size_mm()
    offset = col((x, y, 0))/2
    deltas += offset
    mm_panel_coords = flex.vec2_double(deltas.parts()[0], deltas.parts()[1])

    lab_coords = panel.get_lab_coord(mm_panel_coords)

    ax.scatter(lab_coords.parts()[0], lab_coords.parts()[1], c = data, norm=norm, cmap = cmap, linewidths=0, s=self.params.dot_size)

    return sm, color_vals

  def plot_obs_colored_by_deltas(self, reflections, panel = None, ax = None, bounds = None):
    assert panel is not None and ax is not None and bounds is not None

    data = (reflections['xyzcal.mm']-reflections['xyzobs.mm.value']).norms()
    norm, cmap, color_vals, sm = self.get_normalized_colors(data)
    mm_panel_coords = flex.vec2_double(reflections['xyzobs.mm.value'].parts()[0], reflections['xyzobs.mm.value'].parts()[1])
    lab_coords = panel.get_lab_coord(mm_panel_coords)

    ax.scatter(lab_coords.parts()[0], lab_coords.parts()[1], c = data, norm=norm, cmap = cmap, linewidths=0, s=self.params.dot_size)

    return sm, color_vals

  def plot_obs_colored_by_deltapsi(self, reflections, panel = None, ax = None, bounds = None):
    assert panel is not None and ax is not None and bounds is not None
    data = reflections['delpsical.rad'] * (180/math.pi)
    norm, cmap, color_vals, sm = self.get_normalized_colors(data, vmin=-0.1, vmax=0.1)
    deltas = (reflections['xyzcal.mm']-reflections['xyzobs.mm.value'])*self.delta_scalar

    x, y = panel.get_image_size_mm()
    offset = col((x, y, 0))/2
    deltas += offset
    mm_panel_coords = flex.vec2_double(deltas.parts()[0], deltas.parts()[1])

    lab_coords = panel.get_lab_coord(mm_panel_coords)

    ax.scatter(lab_coords.parts()[0], lab_coords.parts()[1], c = data, norm=norm, cmap = cmap, linewidths=0, s=self.params.dot_size)

    return sm, color_vals

  def plot_radial_displacements_vs_deltapsi(self, reflections, panel = None, ax = None, bounds = None):
    assert panel is not None and ax is not None and bounds is not None
    data = reflections['difference_vector_norms']
    norm, cmap, color_vals, sm = self.get_normalized_colors(data)

    a = reflections['delpsical.rad']*180/math.pi
    b = reflections['radial_displacements']

    fake_coords = flex.vec2_double(a, b) * self.delta_scalar

    x, y = panel.get_image_size_mm()
    offset = col((x, y))/2

    lab_coords = fake_coords + panel.get_lab_coord(offset)[0:2]

    ax.scatter(lab_coords.parts()[0], lab_coords.parts()[1], c = data, norm=norm, cmap = cmap, linewidths=0, s=self.params.dot_size)

    return sm, color_vals

  def plot_unitcells(self, experiments):
    if len(experiments) == 1:
      return
    all_a = flex.double()
    all_b = flex.double()
    all_c = flex.double()
    for crystal in experiments.crystals():
      a, b, c = crystal.get_unit_cell().parameters()[0:3]
      all_a.append(a); all_b.append(b); all_c.append(c)

    fig, axes = plt.subplots(nrows=3, ncols=1)
    for ax, axis, data in zip(axes, ['A', 'B', 'C'], [all_a, all_b, all_c]):
      h = flex.histogram(data,n_slots=50)
      ax.plot(h.slot_centers().as_numpy_array(),h.slots().as_numpy_array(),'-')
      ax.set_title("%s axis histogram. Mean: %7.2f Stddev: %7.2f"%(axis,
                                                                 flex.mean(data),
                                                                 flex.mean_and_variance(data).unweighted_sample_standard_deviation()))
      ax.set_ylabel("N lattices")
      ax.set_xlabel(r"$\AA$")
    plt.tight_layout()

  def plot_histograms(self, reflections, panel = None, ax = None, bounds = None):
    data = reflections['difference_vector_norms']
    colors = ['b-', 'g-', 'g--', 'r-', 'b-', 'b--']
    n_slots = 20
    if self.params.residuals.histogram_max is None:
      h = flex.histogram(data, n_slots=n_slots)
    else:
      h = flex.histogram(data.select(data <= self.params.residuals.histogram_max), n_slots=n_slots)

    n = len(reflections)
    rmsd_obs = math.sqrt((reflections['xyzcal.mm']-reflections['xyzobs.mm.value']).sum_sq()/n)
    sigma = mode = h.slot_centers()[list(h.slots()).index(flex.max(h.slots()))]
    mean_obs = flex.mean(data)
    median = flex.median(data)
    mean_rayleigh = math.sqrt(math.pi/2)*sigma
    rmsd_rayleigh = math.sqrt(2)*sigma

    data = flex.vec2_double([(i,j) for i, j in zip(h.slot_centers(), h.slots())])
    n = len(data)
    for i in [mean_obs, mean_rayleigh, mode, rmsd_obs, rmsd_rayleigh]:
      data.extend(flex.vec2_double([(i, 0), (i, flex.max(h.slots()))]))
    data = self.get_bounded_data(data, bounds)
    tmp = [data[:n]]
    for i in xrange(len(colors)):
      tmp.append(data[n+(i*2):n+((i+1)*2)])
    data = tmp

    for d, c in zip(data, colors):
      ax.plot(d.parts()[0], d.parts()[1], c)

    if ax.get_legend() is None:
      ax.legend([r"$\Delta$XY", "MeanObs", "MeanRayl", "Mode", "RMSDObs", "RMSDRayl"])

  def plot_cdf_manually(self, reflections, panel = None, ax = None, bounds = None):
    colors = ['blue', 'green']
    r = (reflections['xyzcal.mm']-reflections['xyzobs.mm.value']).norms()
    h = flex.histogram(r)
    sigma = h.slot_centers()[list(h.slots()).index(flex.max(h.slots()))] # mode

    x_extent = max(r)
    y_extent = len(r)
    xobs = [i/x_extent for i in sorted(r)]
    yobs = [i/y_extent for i in xrange(y_extent)]
    obs = [(x, y) for x, y in zip(xobs, yobs)]

    ncalc = 100
    xcalc = [i/ncalc for i in xrange(ncalc)]
    ycalc = [1-math.exp((-i**2)/(2*(sigma**2))) for i in xcalc]
    calc = [(x, y) for x, y in zip(xcalc, ycalc)]

    data = [flex.vec2_double(obs),
            flex.vec2_double(calc)]
    if bounds is None:
      ax.set_xlim((-1,1))
      ax.set_ylim((-1,1))
      ax.set_title("%s Outlier SP Manually"%self.params.tag)
    if bounds is not None:
      data = [self.get_bounded_data(d, bounds) for d in data]

    if ax is None:
      fig = plt.figure()
      ax = fig.add_subplot(111)

    for subset,c in zip(data, colors):
        ax.plot(subset.parts()[0], subset.parts()[1], '-', c=c)

  def get_bounded_data(self, data, bounds):
    assert len(bounds) == 4
    x = [b[0] for b in bounds]
    y = [b[1] for b in bounds]
    left = sorted(x)[1]
    right = sorted(x)[2]
    top = sorted(y)[2]
    bottom = sorted(y)[1]
    origin = col((left, bottom))
    scale_x = right-left
    scale_y = top-bottom
    scale = min(scale_x, scale_y)

    data_max_x = flex.max(data.parts()[0])
    data_min_x = flex.min(data.parts()[0])
    data_max_y = flex.max(data.parts()[1])
    data_min_y = flex.min(data.parts()[1])
    data_scale_x = data_max_x - data_min_x
    data_scale_y = data_max_y - data_min_y

    if data_scale_x == 0 or data_scale_y == 0:
      print "WARNING bad scale"
      return data

    return flex.vec2_double(data.parts()[0] * (scale/abs(data_scale_x)),
                            data.parts()[1] * (scale/abs(data_scale_y))) + origin

  def run(self):
    ''' Parse the options. '''
    from dials.util.options import flatten_experiments, flatten_reflections
    # Parse the command line arguments
    params, options = self.parser.parse_args(show_diff_phil=True)
    self.params = params
    experiments = flatten_experiments(params.input.experiments)

    # Find all detector objects
    detectors = experiments.detectors()

    # Verify inputs
    if len(params.input.reflections) == len(detectors) and len(detectors) > 1:
      # case for passing in multiple images on the command line
      assert len(params.input.reflections) == len(detectors)
      reflections = flex.reflection_table()
      for expt_id in xrange(len(detectors)):
        subset = params.input.reflections[expt_id].data
        subset['id'] = flex.int(len(subset), expt_id)
        reflections.extend(subset)
    else:
      # case for passing in combined experiments and reflections
      reflections = flatten_reflections(params.input.reflections)[0]

    detector = detectors[0]

    #from dials.algorithms.refinement.prediction import ExperimentsPredictor
    #ref_predictor = ExperimentsPredictor(experiments, force_stills=experiments.all_stills())

    reflections['difference_vector_norms'] = (reflections['xyzcal.mm']-reflections['xyzobs.mm.value']).norms()

    n = len(reflections)
    rmsd = math.sqrt((reflections['xyzcal.mm']-reflections['xyzobs.mm.value']).sum_sq()/n)
    print "Dataset RMSD", rmsd * 1000

    if params.tag is None:
      tag = ''
    else:
      tag = '%s '%params.tag

    # set up delta-psi ratio heatmap
    p = flex.int() # positive
    n = flex.int() # negative
    for i in set(reflections['id']):
      exprefls = reflections.select(reflections['id']==i)
      p.append(len(exprefls.select(exprefls['delpsical.rad']>0)))
      n.append(len(exprefls.select(exprefls['delpsical.rad']<0)))
    plt.hist2d(p, n, bins=30)
    cb = plt.colorbar()
    cb.set_label("N images")
    plt.title(r"%s2D histogram of pos vs. neg $\Delta\Psi$ per image"%tag)
    plt.xlabel(r"N reflections with $\Delta\Psi$ > 0")
    plt.ylabel(r"N reflections with $\Delta\Psi$ < 0")

    self.delta_scalar = 50

    # Iterate through the detectors, computing detector statistics at the per-panel level (IE one statistic per panel)
    # Per panel dictionaries
    rmsds = {}
    refl_counts = {}
    transverse_rmsds = {}
    radial_rmsds = {}
    ttdpcorr = {}
    # per panelgroup flex arrays
    pg_rmsds = flex.double()
    pg_r_rmsds = flex.double()
    pg_t_rmsds = flex.double()
    pg_refls_count = flex.int()
    table_header = ["PG id", "RMSD","Radial", "Transverse", "N refls"]
    table_header2 = ["","","RMSD","RMSD",""]
    table_data = []
    table_data.append(table_header)
    table_data.append(table_header2)

    # Compute a set of radial and transverse displacements for each reflection
    print "Setting up stats..."
    tmp = flex.reflection_table()
    # Need to construct a variety of vectors
    for panel_id, panel in enumerate(detector):
      panel_refls = reflections.select(reflections['panel'] == panel_id)
      bcl = flex.vec3_double()
      tto = flex.double()
      ttc = flex.double()
      # Compute the beam center in lab space (a vector pointing from the origin to where the beam would intersect
      # the panel, if it did intersect the panel)
      for expt_id in set(panel_refls['id']):
        beam = experiments[expt_id].beam
        s0 = beam.get_s0()
        expt_refls = panel_refls.select(panel_refls['id'] == expt_id)
        beam_centre = panel.get_beam_centre_lab(s0)
        bcl.extend(flex.vec3_double(len(expt_refls), beam_centre))
        obs_x, obs_y, _ = expt_refls['xyzobs.px.value'].parts()
        cal_x, cal_y, _ = expt_refls['xyzcal.px'].parts()
        tto.extend(flex.double([panel.get_two_theta_at_pixel(s0, (obs_x[i], obs_y[i])) for i in xrange(len(expt_refls))]))
        ttc.extend(flex.double([panel.get_two_theta_at_pixel(s0, (cal_x[i], cal_y[i])) for i in xrange(len(expt_refls))]))
      panel_refls['beam_centre_lab'] = bcl
      panel_refls['two_theta_obs'] = tto
      panel_refls['two_theta_cal'] = ttc #+ (0.5*panel_refls['delpsical.rad']*panel_refls['two_theta_obs'])
      # Compute obs in lab space
      x, y, _ = panel_refls['xyzobs.mm.value'].parts()
      c = flex.vec2_double(x, y)
      panel_refls['obs_lab_coords'] = panel.get_lab_coord(c)
      # Compute deltaXY in panel space. This vector is relative to the panel origin
      x, y, _ = (panel_refls['xyzcal.mm'] - panel_refls['xyzobs.mm.value']).parts()
      # Convert deltaXY to lab space, subtracting off of the panel origin
      panel_refls['delta_lab_coords'] = panel.get_lab_coord(flex.vec2_double(x,y)) - panel.get_origin()
      tmp.extend(panel_refls)
    reflections = tmp
    # The radial vector points from the center of the reflection to the beam center
    radial_vectors = (reflections['obs_lab_coords'] - reflections['beam_centre_lab']).each_normalize()
    # The transverse vector is orthogonal to the radial vector and the beam vector
    transverse_vectors = radial_vectors.cross(reflections['beam_centre_lab']).each_normalize()
    # Compute the raidal and transverse components of each deltaXY
    reflections['radial_displacements']     = reflections['delta_lab_coords'].dot(radial_vectors)
    reflections['transverse_displacements'] = reflections['delta_lab_coords'].dot(transverse_vectors)

    # Iterate through the detector at the specified hierarchy level
    for pg_id, pg in enumerate(iterate_detector_at_level(detector.hierarchy(), 0, params.hierarchy_level)):
      pg_msd_sum = 0
      pg_r_msd_sum = 0
      pg_t_msd_sum = 0
      pg_refls = 0
      for p in iterate_panels(pg):
        panel_id = id_from_name(detector, p.get_name())
        panel_refls = reflections.select(reflections['panel'] == panel_id)
        n = len(panel_refls)
        pg_refls += n
        refl_counts[p.get_name()] = n

        if n == 0:
          rmsds[p.get_name()] = -1
          radial_rmsds[p.get_name()] = -1
          transverse_rmsds[p.get_name()] = -1
          continue

        delta_x = panel_refls['xyzcal.mm'].parts()[0] - panel_refls['xyzobs.mm.value'].parts()[0]
        delta_y = panel_refls['xyzcal.mm'].parts()[1] - panel_refls['xyzobs.mm.value'].parts()[1]

        tmp = flex.sum((delta_x**2)+(delta_y**2))
        rmsds[p.get_name()] = math.sqrt(tmp/n) * 1000
        pg_msd_sum += tmp

        r = panel_refls['radial_displacements']
        t = panel_refls['transverse_displacements']
        radial_rmsds[p.get_name()]     = math.sqrt(flex.sum_sq(r)/len(r)) * 1000
        transverse_rmsds[p.get_name()] = math.sqrt(flex.sum_sq(t)/len(t)) * 1000
	pg_r_msd_sum += flex.sum_sq(r)
	pg_t_msd_sum += flex.sum_sq(t)

        a = panel_refls['delpsical.rad']*180/math.pi
        b = panel_refls['two_theta_obs'] - panel_refls['two_theta_cal']
        lc = flex.linear_correlation(a, b)
        ttdpcorr[p.get_name()] = lc.coefficient()


      pg_rmsd = math.sqrt(pg_msd_sum/pg_refls) * 1000
      pg_r_rmsd = math.sqrt(pg_r_msd_sum/pg_refls) * 1000
      pg_t_rmsd = math.sqrt(pg_t_msd_sum/pg_refls) * 1000
      pg_rmsds.append(pg_rmsd)
      pg_r_rmsds.append(pg_r_rmsd)
      pg_t_rmsds.append(pg_t_rmsd)
      pg_refls_count.append(pg_refls)
      table_data.append(["%d"%pg_id, "%.1f"%pg_rmsd, "%.1f"%pg_r_rmsd, "%.1f"%pg_t_rmsd, "%6d"%pg_refls])

    r1 = ["Weighted mean"]
    r2 = ["Weighted stddev"]
    stats = flex.mean_and_variance(pg_rmsds, pg_refls_count.as_double())
    r1.append("%.1f"%stats.mean())
    r2.append("%.1f"%stats.gsl_stats_wsd())
    stats = flex.mean_and_variance(pg_r_rmsds, pg_refls_count.as_double())
    r1.append("%.1f"%stats.mean())
    r2.append("%.1f"%stats.gsl_stats_wsd())
    stats = flex.mean_and_variance(pg_t_rmsds, pg_refls_count.as_double())
    r1.append("%.1f"%stats.mean())
    r2.append("%.1f"%stats.gsl_stats_wsd())
    r1.append("")
    r2.append("")
    table_data.append(r1)
    table_data.append(r2)
    table_data.append(["Mean", "", "", "", "%8.1f"%flex.mean(pg_refls_count.as_double())])

    from libtbx import table_utils
    print "Detector statistics.  Angles in degrees, RMSDs in microns"
    print table_utils.format(table_data,has_header=2,justify='center',delim=" ")

    if params.show_plots:
      if self.params.tag is None:
        t = ""
      else:
        t = "%s "%self.params.tag

      # Plots! these are plots with callbacks to draw on individual panels
      self.detector_plot_refls(detector, reflections, reflections['difference_vector_norms'], '%sDifference vector norms (mm)'%tag, show=False, plot_callback=self.plot_obs_colored_by_deltas)
      self.detector_plot_refls(detector, reflections, reflections['difference_vector_norms'], r'%s$\Delta\Psi$'%tag, show=False, plot_callback=self.plot_obs_colored_by_deltapsi, colorbar_units=r"$\circ$")
      self.detector_plot_refls(detector, reflections, reflections['difference_vector_norms'], r'%s$\Delta$XY*%s'%(tag, self.delta_scalar), show=False, plot_callback=self.plot_deltas)
      self.detector_plot_refls(detector, reflections, reflections['difference_vector_norms'], '%sSP Manual CDF'%tag, show=False, plot_callback=self.plot_cdf_manually)
      self.detector_plot_refls(detector, reflections, reflections['difference_vector_norms'], r'%s$\Delta$XY Histograms'%tag, show=False, plot_callback=self.plot_histograms)
      self.detector_plot_refls(detector, reflections, reflections['difference_vector_norms'], r'%sRadial displacements vs. $\Delta\Psi$, colored by $\Delta$XY'%tag, show=False, plot_callback=self.plot_radial_displacements_vs_deltapsi)

      # Plot intensity vs. radial_displacement
      fig = plt.figure()
      panel_id = 15
      panel_refls = reflections.select(reflections['panel'] == panel_id)
      a = panel_refls['radial_displacements']
      b = panel_refls['intensity.sum.value']
      sel = (a > -0.2) & (a < 0.2) & (b < 50000)
      plt.hist2d(a.select(sel), b.select(sel), bins=100)
      plt.title("%s2D histogram of intensity vs. radial displacement for panel %d"%(tag, panel_id))
      plt.xlabel("Radial displacement (mm)")
      plt.ylabel("Intensity")
      ax = plt.colorbar()
      ax.set_label("Counts")

      # Plot delta 2theta vs. deltapsi
      a = reflections['delpsical.rad']*180/math.pi
      b = reflections['two_theta_obs'] - reflections['two_theta_cal']
      fig = plt.figure()
      sel = (a > -0.2) & (a < 0.2) & (b > -0.002) & (b < 0.002)
      plt.hist2d(a.select(sel), b.select(sel), bins=50)
      cb = plt.colorbar()
      cb.set_label("N reflections")
      plt.title(r'%s$\Delta2\Theta$ vs. $\Delta\Psi$. Showing %d of %d refls'%(tag,len(a.select(sel)),len(a)))
      plt.xlabel(r'$\Delta\Psi \circ$')
      plt.ylabel(r'$\Delta2\Theta \circ$')

      # Plot delta 2theta vs. 2theta
      a = reflections['two_theta_obs']
      b = reflections['two_theta_obs'] - reflections['two_theta_cal']
      fig = plt.figure()
      sel = (b > -0.002) & (b < 0.002)
      plt.hist2d(a.select(sel), b.select(sel), bins=100)
      cb = plt.colorbar()
      cb.set_label("N reflections")
      plt.title(r'%s$\Delta2\Theta$ vs. 2$\Theta$. Showing %d of %d refls'%(tag,len(a.select(sel)),len(a)))
      plt.xlabel(r'2$\Theta \circ$')
      plt.ylabel(r'$\Delta2\Theta \circ$')

      # Plots with single values per panel
      self.detector_plot_dict(detector, refl_counts, u"%s N reflections"%t, u"%6d", show=False)
      self.detector_plot_dict(detector, rmsds, "%s Positional RMSDs (microns)"%t, u"%4.1f", show=False)
      self.detector_plot_dict(detector, radial_rmsds, "%s Radial RMSDs (microns)"%t, u"%4.1f", show=False)
      self.detector_plot_dict(detector, transverse_rmsds, "%s Transverse RMSDs (microns)"%t, u"%4.1f", show=False)
      self.detector_plot_dict(detector, ttdpcorr, r"%s $\Delta2\Theta$ vs. $\Delta\Psi$ CC"%t, u"%5.3f", show=False)

      self.histogram(reflections, '%sDifference vector norms (mm)'%tag)
      self.plot_unitcells(experiments)

      plt.show()

  def histogram(self, reflections, title):
    data = reflections['difference_vector_norms']
    n_slots = 100
    if self.params.residuals.histogram_max is None:
      h = flex.histogram(data, n_slots=n_slots)
    else:
      h = flex.histogram(data.select(data <= self.params.residuals.histogram_max), n_slots=n_slots)

    n = len(reflections)
    rmsd = math.sqrt((reflections['xyzcal.mm']-reflections['xyzobs.mm.value']).sum_sq()/n)
    sigma = mode = h.slot_centers()[list(h.slots()).index(flex.max(h.slots()))]
    mean = flex.mean(data)
    median = flex.median(data)
    print "RMSD", rmsd * 1000
    print "Histogram mode:", mode * 1000
    print "Overall mean:", mean * 1000
    print "Overall median:", median * 1000
    mean2 = math.sqrt(math.pi/2)*sigma
    rmsd2 = math.sqrt(2)*sigma
    print "Rayleigh Mean", mean2 * 1000
    print "Rayleigh RMSD", rmsd2 * 1000

    r = reflections['radial_displacements']
    t = reflections['transverse_displacements']
    print "Overall radial RMSD", math.sqrt(flex.sum_sq(r)/len(r)) * 1000
    print "Overall transverse RMSD", math.sqrt(flex.sum_sq(t)/len(t)) * 1000

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(h.slot_centers().as_numpy_array(), h.slots().as_numpy_array(), '-')

    vmax = self.params.residuals.plot_max
    if self.params.residuals.histogram_xmax is not None:
      ax.set_xlim((0,self.params.residuals.histogram_xmax))
    if self.params.residuals.histogram_ymax is not None:
      ax.set_ylim((0,self.params.residuals.histogram_ymax))
    plt.title(title)


    ax.plot((mean, mean), (0, flex.max(h.slots())), 'g-')
    ax.plot((mean2, mean2), (0, flex.max(h.slots())), 'g--')
    ax.plot((mode, mode), (0, flex.max(h.slots())), 'r-')
    ax.plot((rmsd, rmsd), (0, flex.max(h.slots())), 'b-')
    ax.plot((rmsd2, rmsd2), (0, flex.max(h.slots())), 'b--')

    ax.legend([r"$\Delta$XY", "MeanObs", "MeanRayl", "Mode", "RMSDObs", "RMSDRayl"])

  def detector_plot_refls(self, detector, reflections, data, title, show=True, plot_callback=None, colorbar_units=None):
    """
    Use matplotlib to plot a detector, color coding panels according to data
    @param detector detector reference detector object
    @param data python dictionary of panel names as keys and numbers as values
    @param title title string for plot
    @param units_str string with a formatting statment for units on each panel
    """
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
      bounds = (p0[0:2],p1[0:2],p2[0:2],p3[0:2])

      v1 = p1-p0
      v2 = p3-p0
      vcen = ((v2/2) + (v1/2)) + p0

     # add the panel to the plot
      ax.add_patch(Polygon(bounds, closed=True, fill=False))
      if self.params.panel_numbers:
        ax.annotate("%d"%(panel_id), vcen[0:2], ha='center')

      # find the plot maximum dimensions
      for p in [p0, p1, p2, p3]:
        for c in p[0:2]:
          if abs(c) > max_dim:
            max_dim = abs(c)

      panel_refls = reflections.select(reflections['panel'] == panel_id)
      if len(panel_refls) == 0:
        sm = color_vals = None
        continue

      if plot_callback is None:
        sm = color_vals = None
      else:
        result = plot_callback(panel_refls, panel, ax, bounds)
        if result is not None:
          sm, color_vals = result
        else:
          sm = color_vals = None

    # plot the results
    ax.set_xlim((-max_dim,max_dim))
    ax.set_ylim((-max_dim,max_dim))
    ax.set_xlabel("mm")
    ax.set_ylabel("mm")
    if sm is not None and color_vals is not None:
      if colorbar_units is None:
        colorbar_units = "mm"
      cb = ax.figure.colorbar(sm, ticks=color_vals)
      cb.ax.set_yticklabels(["%3.2f %s"%(i,colorbar_units) for i in color_vals])

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