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

    for level in xrange(4):
      print "Z offsets at level", level
      data = flex.double()
      weights = flex.double()

      for pg_id, (pg1, pg2) in enumerate(zip(iterate_detector_at_level(reference_root, 0, level),
                                             iterate_detector_at_level(moving_root, 0, level))):
        weight = 0
        for panel_id, p in enumerate(iterate_panels(pg2)):
          weight += len(reflections.select(reflections['panel'] == id_from_name(detector, p.get_name())))
        weights.append(weight)

        dists = []
        for pg in [pg1,pg2]:
          ori = pg.get_local_origin()
          dists.append(ori[2]*1000)
        data.append(dists[1]-dists[0])

      if len(data) > 1:
        stats = flex.mean_and_variance(data, weights)
        print r"Weighted mean %5.1f microns"%stats.mean()
        print r"Weighted standard deviation %5.1f microns"%stats.gsl_stats_wsd()
      else:
        print r"Single value %5.1f microns"%data[0]

if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
