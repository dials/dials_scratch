#!/usr/bin/env python
#
# dials.custard.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import absolute_import, division
import libtbx.load_env
import logging
logger = logging.getLogger(libtbx.env.dispatcher_name)
# DIALS_ENABLE_COMMAND_LINE_COMPLETION

help_message = '''

'''

# Create the phil scope
from libtbx.phil import parse
phil_scope = parse(
'''

  output {
    experiments = 'integrated_experiments.json'
      .type = str
      .help = "The experiments output filename"

    reflections = 'integrated.pickle'
      .type = str
      .help = "The integrated output filename"

  }

  include scope dials_scratch.jmp.stills.custard.phil_scope

''', process_includes=True)


class Script(object):
  ''' The integration program. '''

  def __init__(self):
    '''Initialise the script.'''
    from dials.util.options import OptionParser
    import libtbx.load_env

    # The script usage
    usage  = "usage: %s [options] experiment.json" % libtbx.env.dispatcher_name

    # Create the parser
    self.parser = OptionParser(
      usage=usage,
      phil=phil_scope,
      epilog=help_message,
      read_experiments=True)

  def run(self):
    ''' Perform the integration. '''

    from dials.util.options import flatten_reflections, flatten_experiments
    from dxtbx.model.experiment_list import ExperimentListDumper
    from dials_scratch.jmp.stills.custard import Integrator
    from dials.array_family import flex
    import sys

    params, options = self.parser.parse_args(show_diff_phil=False)
    experiments = flatten_experiments(params.input.experiments)

    assert len(experiments) == 1

    integrator = Integrator(experiments[0], params)
    integrator.process()

    reflections = integrator.reflections
    experiments[0] = integrator.experiment

    selection = reflections.get_flags(reflections.flags.integrated_sum)
    partiality = reflections['partiality'].select(selection)
    min_partiality = flex.min(partiality)
    max_partiality = flex.max(partiality)

    print ""
    print "Mosaicity: %f" % integrator.mosaicity
    print "Min partiality: %f, Max partiality: %f" % (
      min_partiality, max_partiality)
    print ""

    from matplotlib import pylab
    pylab.hist(partiality, bins=100)
    pylab.show()

    reflections.as_pickle("integrated.pickle")
    ExperimentListDumper(experiments).as_json("integrated_experiments.json")

if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
