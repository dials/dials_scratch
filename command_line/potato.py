#!/usr/bin/env python
#
# dials.potato.py
#
#  Copyright (C) 2018 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

# LIBTBX_SET_DISPATCHER_NAME dials.potato

from __future__ import absolute_import, division
import libtbx.load_env
from libtbx.utils import Sorry
from libtbx.phil import parse
from dials.util import log
from dials.util.options import OptionParser
from dials.util.options import flatten_reflections, flatten_experiments
from dials_scratch.jmp.potato.potato import Integrator
from dxtbx.model.experiment_list import ExperimentListDumper
import logging


logger = logging.getLogger(libtbx.env.dispatcher_name)

help_message = '''

This script does profile modelling for stills

'''

# Create the phil scope
phil_scope = parse(
'''

''', process_includes=True)


class Script(object):
  ''' 
  The integration program. 
  
  '''

  def __init__(self):
    '''
    Initialise the script.
    
    '''

    # The script usage
    usage  = "usage: %s [options] experiment.json" % libtbx.env.dispatcher_name
    
    # Create the parser
    self.parser = OptionParser(
      usage=usage,
      phil=phil_scope,
      epilog=help_message,
      read_experiments=True,
      read_reflections=True)

  def run(self):
    '''
    Perform the integration. 
    
    '''
    from time import time

    # Check the number of arguments is correct
    start_time = time()

    # Parse the command line
    params, options = self.parser.parse_args(show_diff_phil=False)
    reflections = flatten_reflections(params.input.reflections)
    experiments = flatten_experiments(params.input.experiments)
    if len(reflections) == 0 or len(experiments) == 0:
      self.parser.print_help()
      return
    elif len(reflections) != 1:
      raise Sorry('more than 1 reflection file was given')
    elif len(experiments) == 0:
      raise Sorry('no experiment list was specified')
    reflections = reflections[0]

    # Contruct the integrator
    integrator = Integrator(experiments, reflections)

    # Do the integration
    integrator.initial_integration()
    integrator.refine()
    integrator.predict()
    integrator.integrate()
    
    # Get the reflections
    reflections = integrator.reflections
    experiments = integrator.experiments

    # Save the reflections
    reflections.as_pickle("integrated.pickle")

    dump = ExperimentListDumper(experiments)
    with open("integrated_experiments.json", "w") as outfile:
      outfile.write(dump.as_json())


if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
