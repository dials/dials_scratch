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

help_message = '''

This program is used to calculate statisical measurements of consistency
between observed and predicted reflections

Example:

  libtbx.python detector_residuals.py experiment.json reflections.pickle
'''

# Create the phil parameters
phil_scope = parse('''
tag = None
  .type = str
  .help = Used in the plot titles
colormap = RdYlGn_r
  .type = str
  .help = matplotlib color map. See e.g.: \
          http://matplotlib.org/examples/color/colormaps_reference.html
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
''')

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
    if len(detectors) != 1 or len(reflections) != 1:
      print "Please provide a set of reflections and an experiment list with one detector model"
      return

    detector = detectors[0]
    reflections = reflections[0]

    from dials.algorithms.refinement.prediction import ExperimentsPredictor
    ref_predictor = ExperimentsPredictor(experiments, force_stills=experiments.all_stills())

    reflections = ref_predictor.predict(reflections)
    reflections['residuals'] = (reflections['xyzcal.mm']-reflections['xyzobs.mm.value']).norms()

    if params.tag is None:
      title = 'XY residuals (mm)'
    else:
      title = '%s XY residuals (mm)'%params.tag
    self.detector_plot(detector, reflections, 'residuals', title)

  def detector_plot(self, detector, reflections, data_key, title):
    """
    Use matplotlib to plot a detector, color coding panels according to data
    @param detector detector reference detector object
    @param data python dictionary of panel names as keys and numbers as values
    @param title title string for plot
    @param units_str string with a formatting statment for units on each panel
    """
    if self.params.residuals.histogram_max is None:
      h = flex.histogram(reflections[data_key])
    else:
      h = flex.histogram(reflections[data_key].select(reflections[data_key] <= self.params.residuals.histogram_max))
    print "Histogram mode:", h.slot_centers()[list(h.slots()).index(flex.max(h.slots()))]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(h.slot_centers().as_numpy_array(), h.slots().as_numpy_array(), '-')

    vmax = self.params.residuals.plot_max
    if vmax is not None:
      ax.plot((vmax, vmax), (0, flex.max(h.slots())), 'r-')
    if self.params.residuals.histogram_xmax is not None:
      ax.set_xlim((0,self.params.residuals.histogram_xmax))
    if self.params.residuals.histogram_ymax is not None:
      ax.set_ylim((0,self.params.residuals.histogram_ymax))
    plt.title(title)

    # initialize the color map
    if vmax is None:
      vmax = flex.max(values)
    else:
      sel = reflections[data_key] <= vmax
      reflections = reflections.select(sel)
    values = reflections[data_key]

    norm = Normalize(vmin=flex.min(values), vmax=vmax)
    cmap = plt.cm.get_cmap(self.params.colormap)
    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    color_vals = np.linspace(flex.min(values), vmax, 11)
    sm.set_array(color_vals) # needed for colorbar

    mm_panel_coords = flex.vec2_double(reflections['xyzobs.mm.value'].parts()[0], reflections['xyzobs.mm.value'].parts()[1])

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
      ax.add_patch(Polygon((p0[0:2],p1[0:2],p2[0:2],p3[0:2]), closed=True, fill=False))
      if self.params.panel_numbers:
        ax.annotate("%d"%(panel_id), vcen[0:2], ha='center')

      lab_coords = panel.get_lab_coord(mm_panel_coords.select(reflections['panel'] == panel_id))
      panel_values = values.select(reflections['panel'] == panel_id)

      ax.scatter(lab_coords.parts()[0], lab_coords.parts()[1], c = panel_values, norm=norm, cmap = cmap, linewidths=0, s=self.params.dot_size)

      # find the plot maximum dimensions
      for p in [p0, p1, p2, p3]:
        for c in p[0:2]:
          if abs(c) > max_dim:
            max_dim = abs(c)

    # plot the results, full size
    ax.set_xlim((-max_dim,max_dim))
    ax.set_ylim((-max_dim,max_dim))
    # inset
    #ax.set_xlim((-4,41))
    #ax.set_ylim((4,50))
    ax.set_xlabel("mm")
    ax.set_ylabel("mm")
    cb = fig.colorbar(sm, ticks=color_vals)
    cb.ax.set_yticklabels(["%3.2f mm"%i for i in color_vals])
    plt.title(title)

    plt.show()

if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
