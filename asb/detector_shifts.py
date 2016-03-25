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
from libtbx.phil import parse

help_message = '''

This program is used to show differences between a reference and a moving set of detectors
Example:

  libtbx.python detector_shifts.py experiment1.json experiment2.json reflections1.pickle reflections2.pickle
'''

# Create the phil parameters
phil_scope = parse('''
  include scope dials_scratch.asb.detector_congruence.phil_scope
''', process_includes=True)

from dials_scratch.asb.detector_congruence import iterate_detector_at_level, iterate_panels, id_from_name
from dials_scratch.asb.detector_congruence import Script as ParentScript

class Script(ParentScript):
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
      print "Please provide a reference and a moving set of experiments and or datablocks"
      return

    detector = detectors[1]
    reflections = reflections[1]
    reference_root = detectors[0].hierarchy()
    moving_root = detector.hierarchy()
    rori = col(reference_root.get_origin())
    r_norm = col(reference_root.get_normal())
    s0 = col(flex.vec3_double([col(b.get_s0()) for b in experiments.beams()]).mean())

    table_header = ["Hierarchy", "X Offsets", "X Offsets", "Y Offsets", "Y Offsets", "Delta XY", "Delta XY", "Z Offsets", "Z Offsets" ]
    table_header2 = ["Level","","Sigma","","Sigma","","Sigma", "", "Sigma"]
    table_header3 = ["","(microns)","(microns)","(microns)","(microns)","(microns)","(microns)", "(microns)", "(microns)"]
    table_data = []
    table_data.append(table_header)
    table_data.append(table_header2)
    table_data.append(table_header3)

    for level in xrange(4):
      x_offsets = flex.double()
      y_offsets = flex.double()
      delta_xy = flex.double()
      z_offsets = flex.double()
      weights = flex.double()

      for pg_id, (pg1, pg2) in enumerate(zip(iterate_detector_at_level(reference_root, 0, level),
                                             iterate_detector_at_level(moving_root, 0, level))):
        weight = 0
        for panel_id, p in enumerate(iterate_panels(pg2)):
          weight += len(reflections.select(reflections['panel'] == id_from_name(detector, p.get_name())))
        weights.append(weight)

        x_dists = []
        y_dists = []
        z_dists = []
        ori_xy = []
        for pg in [pg1,pg2]:
          ori = pg.get_local_origin()
          x_dists.append(ori[0]*1000)
          y_dists.append(ori[1]*1000)
          ori_xy.append(col((ori[0], ori[1])))
          z_dists.append(ori[2]*1000)
        x_offsets.append(x_dists[1]-x_dists[0])
        y_offsets.append(y_dists[1]-y_dists[0])
        delta_xy.append((ori_xy[1]-ori_xy[0]).length()*1000)
        z_offsets.append(z_dists[1]-z_dists[0])

      row = ["%d"%level]
      if len(z_offsets) == 0:
        row.extend(["%6.1f"%0]*6)
      elif len(z_offsets) == 1:
        for data in [x_offsets, y_offsets, delta_xy, z_offsets]:
          row.append("%6.1f"%data[0])
          row.append("%6.1f"%0)
      else:
        for data in [x_offsets, y_offsets, delta_xy, z_offsets]:
          stats = flex.mean_and_variance(data, weights)
          row.append("%6.1f"%stats.mean())
          row.append("%6.1f"%stats.gsl_stats_wsd())
      table_data.append(row)

    from libtbx import table_utils
    print "Detector shifts"
    print table_utils.format(table_data,has_header=3,justify='center',delim=" ")

if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
