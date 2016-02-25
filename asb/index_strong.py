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
from libtbx.phil import parse
from libtbx import easy_pickle

help_message = '''
  Script to index strong reflections without outlier rejection
'''

# Create the phil parameters
phil_scope = parse('''
  d_min = None
    .type=float
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
      read_reflections=True,
      epilog=help_message)

  def run(self):
    ''' Parse the options. '''
    from dials.util.options import flatten_experiments, flatten_reflections
    # Parse the command line arguments
    params, options = self.parser.parse_args(show_diff_phil=True)
    self.params = params
    experiments = flatten_experiments(params.input.experiments)
    reflections = flatten_reflections(params.input.reflections)

    assert len(reflections) == len(experiments) == 1
    reflections = reflections[0]
    exp = experiments[0]

    from dials.algorithms.indexing import index_reflections
    from dials.algorithms.indexing.indexer import indexer_base

    reflections['id'] = flex.int(len(reflections), -1)
    reflections['imageset_id'] = flex.int(len(reflections), 0)
    reflections = indexer_base.map_spots_pixel_to_mm_rad(reflections, exp.detector, exp.scan)

    indexer_base.map_centroids_to_reciprocal_space(
      reflections, exp.detector, exp.beam, exp.goniometer,)

    index_reflections(reflections,
                      experiments, params.d_min,
                      tolerance=0.3)
    indexed_reflections = reflections.select(reflections['miller_index'] != (0,0,0))
    print "Indexed %d reflections out of %d"%(len(indexed_reflections), len(reflections))
    easy_pickle.dump("indexedstrong.pickle", indexed_reflections)

if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
