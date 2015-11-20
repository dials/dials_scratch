#!/usr/bin/env python
#
# dials.idials.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division
from cmd import Cmd
import sys


class ActionError(RuntimeError):
  '''
  Class to represent exception for when an action can't be performed

  '''

  def __init__(self, action, parent_action):
    '''
    :param action: The action we want to do
    :param parent_action: The parent's action

    '''
    text = 'Unable to perform "%s" after "%s"' % (action, parent_action)
    super(ActionError, self).__init__(text)


class ExternalCommand(object):
  '''
  Class to run an external command

  '''

  def __init__(self, command, filename, output=sys.stdout, wait_time=0.1):
    '''
    Run the command

    :param command: The command to run
    :param filename: The filename for stdout and stderr
    :param output: Write stdout and stderr to additional output
    :param wait_time: Wait time for polling file output

    '''
    import subprocess
    import time

    # Create the command string
    if not isinstance(command, str):
      command = subprocess.list2cmdline(command)

    # Open the file for output and input
    with open(filename, "w") as outfile:
      with open(filename, "r") as infile:

        # Start the process
        process = subprocess.Popen(
          command,
          stdout=outfile,
          stderr=outfile,
          shell=True,
          universal_newlines=True,
          bufsize=-1)

        # Write the lines of the file to the output
        if output is not None:
          while True:
            line = infile.readline()
            if not line:
              if process.poll() is not None:
                break
              time.sleep(wait_time)
            else:
              output.write(line)

        # Get the result
        self.result = process.wait()


def run_external_command(command, filename, output=sys.stdout, wait_time=0.1):
  '''
  Helper function to run command

  :param command: The command to run
  :param filename: The filename to output
  :param output: The buffer to write to
  :param wait_time: The polling timeout

  '''
  command = ExternalCommand(command, filename, output, wait_time)
  if command.result != 0:
    raise RuntimeError('Error: external command failed')


class ParameterManager(object):
  '''
  A class to manage the current set of parameters.

  '''

  def __init__(self, phil_scope):
    '''
    Create the master phil and set working phil to default parameters

    '''
    self.master_phil = phil_scope
    self.reset()

  def reset(self):
    '''
    Reset the working phil to the default parameters

    '''
    from libtbx.phil import parse
    self.working_phil = self.master_phil.fetch(source=parse(''))

  def set(self, parameters, short_syntax=False):
    '''
    Set a parameter and update the working phil
    :param parameter: The text string of parameters
    :param short_syntax: True/False treat as command line parameter

    '''
    from libtbx.phil import parse
    from libtbx.utils import Sorry
    import shlex
    if short_syntax == True:
      for parameter in shlex.split(parameters):
        interpretor = self.master_phil.command_line_argument_interpreter()
        self.working_phil = self.working_phil.fetch(
          interpretor.process_arg(parameter))
    else:
      working_phil, unused = self.working_phil.fetch(
        source=parse(parameters),
        track_unused_definitions=True)
      if len(unused) > 0:
        msg = [item.object.as_str().strip() for item in unused]
        msg = '\n'.join(['  %s' % line for line in msg])
        raise Sorry('The following definitions were not recognised\n%s' % msg)
      self.working_phil = working_phil

  def get(self, diff=True):
    '''
    Get the phil parameters

    :param diff: Get the diff phil

    '''
    if diff == False:
      result = self.working_phil
    else:
      result = self.master_phil.fetch_diff(source=self.working_phil)
    return result


class ImportParameterManager(ParameterManager):
  '''
  Specialization for import parameters

  '''

  def __init__(self):
    '''
    Import phil scope and set up

    '''
    from libtbx.phil import parse
    phil_scope = parse('''
      include scope dials.command_line.import.phil_scope
    ''', process_includes=True)
    super(ImportParameterManager, self).__init__(phil_scope)


class FindSpotsParameterManager(ParameterManager):
  '''
  Specialization for find spots parameters

  '''

  def __init__(self):
    '''
    Import phil scope and set up

    '''
    from libtbx.phil import parse
    phil_scope = parse('''
      input {
        datablock = None
          .type = str
      }
      include scope dials.command_line.find_spots.phil_scope
    ''', process_includes=True)
    super(FindSpotsParameterManager, self).__init__(phil_scope)


class IndexParameterManager(ParameterManager):
  '''
  Specialization for index parameters

  '''

  def __init__(self):
    '''
    Import phil scope and set up

    '''
    from libtbx.phil import parse
    phil_scope = parse('''
      input {
        datablock = None
          .type = str
        reflections = None
          .type = str
      }
      include scope dials.command_line.index.phil_scope
    ''', process_includes=True)
    super(IndexParameterManager, self).__init__(phil_scope)


class RefineBSParameterManager(ParameterManager):
  '''
  Specialization for refine_bs parameters

  '''

  def __init__(self):
    '''
    Import phil scope and set up

    '''
    from libtbx.phil import parse
    phil_scope = parse('''
      input {
        experiments = None
          .type = str
        reflections = None
          .type = str
      }
      include scope dials.command_line.refine_bravais_settings.phil_scope
    ''', process_includes=True)
    super(RefineBSParameterManager, self).__init__(phil_scope)


class ReIndexParameterManager(ParameterManager):
  '''
  Specialization for reindex parameters

  '''

  def __init__(self):
    '''
    Import phil scope and set up

    '''
    from libtbx.phil import parse
    phil_scope = parse('''
      solution = None
        .type = int
      input {
        reflections = None
          .type = str
      }
      include scope dials.command_line.reindex.phil_scope
    ''', process_includes=True)
    super(ReIndexParameterManager, self).__init__(phil_scope)


class RefineParameterManager(ParameterManager):
  '''
  Specialization for refine parameters

  '''

  def __init__(self):
    '''
    Import phil scope and set up

    '''
    from libtbx.phil import parse
    phil_scope = parse('''
      input {
        experiments = None
          .type = str
        reflections = None
          .type = str
      }
      include scope dials.command_line.refine.phil_scope
    ''', process_includes=True)
    super(RefineParameterManager, self).__init__(phil_scope)


class IntegrateParameterManager(ParameterManager):
  '''
  Specialization for integrate parameters

  '''

  def __init__(self):
    '''
    Import phil scope and set up

    '''
    from libtbx.phil import parse
    phil_scope = parse('''
      input {
        experiments = None
          .type = str
        reflections = None
          .type = str
      }
      include scope dials.command_line.integrate.phil_scope
    ''', process_includes=True)
    super(IntegrateParameterManager, self).__init__(phil_scope)


class ExportParameterManager(ParameterManager):
  '''
  Specialization for export parameters

  '''

  def __init__(self):
    '''
    Import phil scope and set up

    '''
    from libtbx.phil import parse
    phil_scope = parse('''
      input {
        experiments = None
          .type = str
        reflections = None
          .type = str
      }
      include scope dials.command_line.export_mtz.phil_scope
    ''', process_includes=True)
    super(ExportParameterManager, self).__init__(phil_scope)


class GlobalParameterManager(dict):
  '''
  Class to hold all parameter managers

  '''

  def __init__(self):
    '''
    Init everything

    '''
    super(GlobalParameterManager, self).__init__()
    self.update({
      'import'     : ImportParameterManager(),
      'find_spots' : FindSpotsParameterManager(),
      'index'      : IndexParameterManager(),
      'refine_bs'  : RefineBSParameterManager(),
      'reindex'    : ReIndexParameterManager(),
      'refine'     : RefineParameterManager(),
      'integrate'  : IntegrateParameterManager(),
      'export'     : ExportParameterManager(),
    })


class Counter(object):
  '''
  A counter class to update command indices

  '''

  def __init__(self):
    '''
    Counter begins at zero

    '''
    self.count = 0

  def current(self):
    '''
    :return: The current counter

    '''
    return self.count

  def next(self):
    '''
    Update the counter value

    :return: The new counter value

    '''
    result = self.count
    self.count += 1
    return result


class CommandNode(object):
  '''
  A class to represent the commands

  '''

  parent_actions = []

  def __init__(self, parent=None, action='', parameters=None, directory=None):
    '''
    Initialise the tree parent and children and set the index

    :param parent: The command parent

    '''
    from os.path import join
    import copy

    # Check the parent is OK
    if parent is None:
      parent_action = None
    else:
      parent_action = parent.action
    if parent_action not in self.parent_actions:
      raise ActionError(action, parent_action)

    # Raise exception if trying to job after failure
    if parent is not None and parent.success == False:
      raise RuntimeError('Error: parent job %d failed' % parent.index)

    # Set the parent and counter
    self.parent = parent
    if self.parent is None:
      self.counter = Counter()
    else:
      self.counter = self.parent.counter

    # Set the index and some tree stuff
    self.index = self.counter.next()
    self.children = []
    if self.parent is not None:
      self.parent.children.append(self)

    # Init the result
    self.success = False

    # Save the info
    self.action = action
    self.parameters = copy.deepcopy(parameters)
    if directory is not None:
      self.directory = join(directory, "%d_%s" % (self.index, self.action))
    else:
      self.directory = None

  def __iter__(self):
    '''
    Iterate through the children and their children

    '''
    yield self, 0
    for child in self.children:
      for node, depth in child:
        yield node, depth+1

  def apply(self):
    '''
    Apply the command

    :return: True/False success or failure

    '''
    from os.path import exists, join
    from os import makedirs

    # Check the output path does not exist already
    if exists(self.directory):
      raise RuntimeError('Output directory %s already exists' % self.directory)

    # Initialise running the command
    self.initialize()

    # Make the directory to store output
    makedirs(self.directory)

    # Set the parameter filename and write to file
    parameters = join(self.directory, "parameters.phil")
    with open(parameters, "w") as outfile:
      outfile.write(self.parameters.get(diff=True).as_str())
    outfile.close()

    # Set the output filename
    output = join(self.directory, "output.txt")

    # Run the command (override this method)
    self.run(parameters, output)

    # Grab the result (override this method)
    self.finalize()

    # Set success
    self.success = True


class CommandTree(object):
  '''
  A class to provide to helpful tree functions

  '''

  def __init__(self, root):
    '''
    :param root: The tree root

    '''
    self.root = root

  def goto(self, index):
    '''
    Go to the desired node in the tree

    :param index: the index of the node to go to
    :return: The node at index

    '''
    for node, level in self.iternodes():
      if node.index == index:
        return node
    raise IndexError('Node %d not found' % index)

  def iternodes(self):
    '''
    Iterate through the tree nodes depth first

    '''
    for node, depth in self.root:
      yield node, depth

  def string(self, current=None):
    '''
    :return: The tree as a string

    '''
    size = len(str(self.root.counter.current()))
    def draw_tree(node, prefix):
      from cStringIO import StringIO
      buf = StringIO()
      if prefix:
        buf.write(('%%%dd' % size) % node.index)
        buf.write(' %s' % ('S' if node.success else 'F'))
        buf.write(prefix[:-3])
        buf.write('  +--')
      buf.write(node.action)
      if current is not None and node.index == current:
        buf.write(" (current)")
      buf.write('\n')
      for index, child in enumerate(node.children):
        if index+1 == len(node.children):
          sub_prefix = prefix + '   '
        else:
          sub_prefix = prefix + '  |'
        buf.write(draw_tree(child, sub_prefix))
      return buf.getvalue()
    return draw_tree(self.root, '')


class InitialState(CommandNode):
  '''
  A class to represent the initial clean state

  '''

  parent_actions = [None]

  def __init__(self):
    '''
    Initialise the command

    '''
    super(InitialState, self).__init__(None, 'clean')

    # Set success to True
    self.success = True

  def apply(self):
    '''
    Override apply since we have nothing to do here

    '''
    raise RuntimeError("Programming error: nothing to do")


class ImportCommand(CommandNode):
  '''
  A command to perform an import operation

  '''

  parent_actions = ['clean']

  def __init__(self, parent, parameters, directory):
    '''
    Initialise the command node

    :param parent: The parent command
    :param parameters: The parameters to use
    :param directory: The output directory

    '''
    super(ImportCommand, self).__init__(
      parent, 'import', parameters, directory)

  def initialize(self):
    '''
    Initialise the processing

    '''
    from os.path import join
    self.filenames = {
      'output.datablock' : join(self.directory, "datablock.json"),
      'output.log'       : join(self.directory, "info.log"),
      'output.debug_log' : join(self.directory, "debug.log")
    }
    for name, value in self.filenames.iteritems():
      self.parameters.set('%s=%s' % (name, value))

  def run(self, parameters, output):
    '''
    Run the import command

    :param parameters: The input parameter filename
    :param output: The stdout file

    '''
    print "Running import: for output see %s" % output
    run_external_command(['dials.import', parameters], output)

  def finalize(self):
    '''
    Finalize the processing

    '''
    from os.path import exists
    for name, value in self.filenames.iteritems():
      if not exists(value):
        raise RuntimeError("File %s could not be found" % value)


class FindSpotsCommand(CommandNode):
  '''
  A command to perform an find_spots operation

  '''

  parent_actions = ['import']

  def __init__(self, parent, parameters, directory):
    '''
    Initialise the command node

    :param parent: The parent command
    :param parameters: The parameters to use
    :param directory: The output directory

    '''
    super(FindSpotsCommand, self).__init__(
      parent, 'find_spots', parameters, directory)

  def initialize(self):
    '''
    Initialise the processing

    '''
    from os.path import join
    self.filenames = {
      'input.datablock'    : self.parent.filenames['output.datablock'],
      'output.datablock'   : join(self.directory, "datablock.json"),
      'output.reflections' : join(self.directory, "reflections.pickle"),
      'output.log'         : join(self.directory, "info.log"),
      'output.debug_log'   : join(self.directory, "debug.log"),
    }
    for name, value in self.filenames.iteritems():
      self.parameters.set('%s=%s' % (name, value))

  def run(self, parameters, output):
    '''
    Run the find_spots command

    :param parameters: The input parameter filename
    :param output: The stdout file

    '''
    print "Running find_spots: for output see %s" % output
    run_external_command(['dials.find_spots', parameters], output)

  def finalize(self):
    '''
    Finalize the processing

    '''
    from os.path import exists
    for name, value in self.filenames.iteritems():
      if not exists(value):
        raise RuntimeError("File %s could not be found" % value)


class IndexCommand(CommandNode):
  '''
  A command to perform an index operation

  '''

  parent_actions = ['find_spots']

  def __init__(self, parent, parameters, directory):
    '''
    Initialise the command node

    :param parent: The parent command
    :param parameters: The parameters to use
    :param directory: The output directory

    '''
    super(IndexCommand, self).__init__(
      parent, 'index', parameters, directory)

  def initialize(self):
    '''
    Initialise the processing

    '''
    from os.path import join
    self.filenames = {
      'input.datablock'    : self.parent.filenames['output.datablock'],
      'input.reflections'  : self.parent.filenames['output.reflections'],
      'output.reflections' : join(self.directory, "reflections.pickle"),
      'output.experiments' : join(self.directory, "experiments.json"),
      'output.log'         : join(self.directory, "info.log"),
      'output.debug_log'   : join(self.directory, "debug.log"),
    }
    for name, value in self.filenames.iteritems():
      self.parameters.set('%s=%s' % (name, value))

  def run(self, parameters, output):
    '''
    Run the index command

    :param parameters: The input parameter filename
    :param output: The stdout file

    '''
    print "Running index: for output see %s" % output
    run_external_command(['dials.index', parameters], output)

  def finalize(self):
    '''
    Finalize the processing

    '''
    from os.path import exists
    for name, value in self.filenames.iteritems():
      if not exists(value):
        raise RuntimeError("File %s could not be found" % value)


class RefineBSCommand(CommandNode):
  '''
  A command to perform an refine_bs operation

  '''

  parent_actions = ['index']

  def __init__(self, parent, parameters, directory):
    '''
    Initialise the command node

    :param parent: The parent command
    :param parameters: The parameters to use
    :param directory: The output directory

    '''
    super(RefineBSCommand, self).__init__(
      parent, 'refine_bs', parameters, directory)

  def initialize(self):
    '''
    Initialise the processing

    '''
    from os.path import join
    self.filenames = {
      'input.experiments'  : self.parent.filenames['output.experiments'],
      'input.reflections'  : self.parent.filenames['output.reflections'],
      'output.log'         : join(self.directory, "info.log"),
      'output.debug_log'   : join(self.directory, "debug.log"),
      'output.directory'   : self.directory,
    }
    for name, value in self.filenames.iteritems():
      self.parameters.set('%s=%s' % (name, value))

  def run(self, parameters, output):
    '''
    Run the refine_bravais_settings command

    :param parameters: The input parameter filename
    :param output: The stdout file

    '''
    print "Running refine_bs: for output see %s" % output
    run_external_command(['dials.refine_bravais_settings', parameters], output)

  def finalize(self):
    '''
    Finalize the processing

    '''
    from os.path import exists, join
    import json

    # Add expected output files to check
    self.filenames.update({
      'output.reflections' : self.parent.filenames['output.reflections'],
      'output.summary'     : join(self.directory, 'bravais_summary.json'),
    })

    # Check the output files
    for name, value in self.filenames.iteritems():
      if not exists(value):
        raise RuntimeError("File %s could not be found" % value)

    # Read the summary and check all json files exist
    with open(self.filenames['output.summary']) as summary_file:
      self.summary = json.load(summary_file)
      self.bs_filenames = {}
      for name, value in self.summary.iteritems():
        self.bs_filenames[name] = join(self.directory, 'bravais_setting_%s.json' % name)
      for name, value in self.bs_filenames.iteritems():
        if not exists(value):
          raise RuntimeError("File %s could not be found" % value)


class ReIndexCommand(CommandNode):
  '''
  A command to perform an reindex operation

  '''

  parent_actions = ['refine_bs']

  def __init__(self, parent, parameters, directory):
    '''
    Initialise the command node

    :param parent: The parent command
    :param parameters: The parameters to use
    :param directory: The output directory

    '''
    super(ReIndexCommand, self).__init__(
      parent, 'reindex', parameters, directory)

  def initialize(self):
    '''
    Initialise the processing

    '''
    from os.path import join

    # The files which can be set as parameters
    self.filenames = {
      'input.reflections'  : self.parent.filenames['output.reflections'],
      'output.reflections' : join(self.directory, "reflections.pickle"),
    }
    for name, value in self.filenames.iteritems():
      self.parameters.set('%s=%s' % (name, value))

    # Get the solution we want and convert to the change_of_basis_op
    solution = self.parameters.get(diff=False).extract().solution
    if solution is None:
      raise RuntimeError("No solution selected")
    change_of_basis_op = self.parent.summary[str(solution)]['cb_op']

    # Set the solution parameter to None and set the cb_op
    self.parameters.set("solution=None")
    self.parameters.set("change_of_basis_op=%s" % change_of_basis_op)

    # Set the output experiments to the bravais settings file
    self.filenames.update({
      'output.experiments' : self.parent.bs_filenames[str(solution)]
    })

  def run(self, parameters, output):
    '''
    Run the index command

    :param parameters: The input parameter filename
    :param output: The stdout file

    '''
    print "Running index: for output see %s" % output
    run_external_command(['dials.reindex', parameters], output)

  def finalize(self):
    '''
    Finalize the processing

    '''
    from os.path import exists
    for name, value in self.filenames.iteritems():
      if not exists(value):
        raise RuntimeError("File %s could not be found" % value)


class RefineCommand(CommandNode):
  '''
  A command to perform an refine operation

  '''

  parent_actions = ['index', 'reindex', 'refine', 'integrate']

  def __init__(self, parent, parameters, directory):
    '''
    Initialise the command node

    :param parent: The parent command
    :param parameters: The parameters to use
    :param directory: The output directory

    '''
    super(RefineCommand, self).__init__(
      parent, 'refine', parameters, directory)

  def initialize(self):
    '''
    Initialise the processing

    '''
    from os.path import join
    self.filenames = {
      'input.experiments'  : self.parent.filenames['output.experiments'],
      'input.reflections'  : self.parent.filenames['output.reflections'],
      'output.reflections' : join(self.directory, "reflections.pickle"),
      'output.experiments' : join(self.directory, "experiments.json"),
      'output.log'         : join(self.directory, "info.log"),
      'output.debug_log'   : join(self.directory, "debug.log"),
      'output.matches'     : join(self.directory, "matches.pickle"),
      'output.centroids'   : join(self.directory, "centroids.txt"),
      'output.history'     : join(self.directory, "history.txt"),
    }
    for name, value in self.filenames.iteritems():
      self.parameters.set('%s=%s' % (name, value))

  def run(self, parameters, output):
    '''
    Run the refine command

    :param parameters: The input parameter filename
    :param output: The stdout file

    '''
    print "Running refine: for output see %s" % output
    run_external_command(['dials.refine', parameters], output)

  def finalize(self):
    '''
    Finalize the processing

    '''
    from os.path import exists
    for name, value in self.filenames.iteritems():
      if not exists(value):
        raise RuntimeError("File %s could not be found" % value)


class IntegrateCommand(CommandNode):
  '''
  A command to perform an integrate operation

  '''

  parent_actions = ['index', 'reindex', 'refine', 'integrate']

  def __init__(self, parent, parameters, directory):
    '''
    Initialise the command node

    :param parent: The parent command
    :param parameters: The parameters to use
    :param directory: The output directory

    '''
    super(IntegrateCommand, self).__init__(
      parent, 'integrate', parameters, directory)

  def initialize(self):
    '''
    Initialise the processing

    '''
    from os.path import join
    self.filenames = {
      'input.experiments'  : self.parent.filenames['output.experiments'],
      'input.reflections'  : self.parent.filenames['output.reflections'],
      'output.reflections' : join(self.directory, "reflections.pickle"),
      'output.experiments' : join(self.directory, "experiments.json"),
      'output.log'         : join(self.directory, "info.log"),
      'output.debug_log'   : join(self.directory, "debug.log"),
      'output.report'      : join(self.directory, "report.json"),
      'output.phil'        : 'None'
    }
    for name, value in self.filenames.iteritems():
      self.parameters.set('%s=%s' % (name, value))
  def run(self, parameters, output):
    '''
    Run the integrate command

    :param parameters: The input parameter filename
    :param output: The stdout file

    '''
    print "Running integrate: for output see %s" % output
    run_external_command(['dials.integrate', parameters], output)

  def finalize(self):
    '''
    Finalize the processing

    '''
    from os.path import exists
    for name, value in self.filenames.iteritems():
      if value is not 'None' and not exists(value):
        raise RuntimeError("File %s could not be found" % value)


class ExportCommand(CommandNode):
  '''
  A command to perform an export operation

  '''

  parent_actions = ['integrate']

  def __init__(self, parent, parameters, directory):
    '''
    Initialise the command node

    :param parent: The parent command
    :param parameters: The parameters to use
    :param directory: The output directory

    '''
    super(ExportCommand, self).__init__(
      parent, 'export', parameters, directory)

  def initialize(self):
    '''
    Initialise the processing

    '''
    from os.path import join
    self.filenames = {
      'input.experiments'  : self.parent.filenames['output.experiments'],
      'input.reflections'  : self.parent.filenames['output.reflections'],
      'hklout'             : join(self.directory, "reflections.mtz"),
      'output.log'         : join(self.directory, "info.log"),
      'output.debug_log'   : join(self.directory, "debug.log"),
    }
    for name, value in self.filenames.iteritems():
      self.parameters.set('%s=%s' % (name, value))

  def run(self, parameters, output):
    '''
    Run the export command

    :param parameters: The input parameter filename
    :param output: The stdout file

    '''
    print "Running export: for output see %s" % output
    run_external_command(['dials.export_mtz', parameters], output)

  def finalize(self):
    '''
    Finalize the processing

    '''
    from os.path import exists
    import shutil
    for name, value in self.filenames.iteritems():
      print name, value
      if not exists(value):
        raise RuntimeError("File %s could not be found" % value)

    # Copy the resulting mtz file to the working directory
    result_filename = "%d_integrated.mtz" % self.index
    shutil.copy2(self.filenames['hklout'], result_filename)


class ApplicationState(object):
  '''
  A class to hold all the application state

  '''

  def __init__(self, directory):
    '''
    Initialise the state

    :param directory: The output directory

    '''
    # Create the parameters
    self.parameters = GlobalParameterManager()

    # Set the initial state to current
    self.current = InitialState()

    # Create the command tree
    self.command_tree = CommandTree(self.current)

    # Save the parameters and directory
    self.directory = directory

    # Initialise the mode
    self.mode = 'import'

  def run(self):
    '''
    Run the command for the given mode

    '''

    # The command classes
    CommandClass = {
      'import'                  : ImportCommand,
      'find_spots'              : FindSpotsCommand,
      'index'                   : IndexCommand,
      'refine_bs' : RefineBSCommand,
      'reindex'                 : ReIndexCommand,
      'refine'                  : RefineCommand,
      'integrate'               : IntegrateCommand,
      'export'                  : ExportCommand
    }

    # Create the command
    command = CommandClass[self.mode](
      self.current,
      self.parameters[self.mode],
      self.directory)

    # Apply the command
    command.apply()

    # If successful update current
    self.current = command

  def goto(self, index):
    '''
    Goto a specific command

    :param index: The command index

    '''
    self.current = self.command_tree.goto(index)

  def history(self):
    '''
    Get the command history

    :return The command history

    '''
    return self.command_tree.string(current=self.current.index)

  def dump(self, filename):
    '''
    Dump the state to file

    :param filename: The filename

    '''
    import cPickle as pickle
    with open(filename, "w") as outfile:
      pickle.dump(self, outfile)

  @classmethod
  def load(Class, filename):
    '''
    Load the state from file

    :param filename: The filename
    :return: The state object

    '''
    import cPickle as pickle
    with open(filename) as infile:
      return pickle.load(infile)


class Controller(object):
  '''
  The controller class.

  This defines the interface the DIALS GUI and CLI programs can use to interact
  with the DIALS programs in a standard way.

  '''

  # The list of program modes
  mode_list = [
    'import',
    'find_spots',
    'index',
    'refine_bs',
    'reindex',
    'refine',
    'integrate',
    'export']

  def __init__(self,
               directory=".",
               state_filename="dials.state",
               recover=True):
    '''
    Initialise the controller

    :param directory: The output directory
    :param state_filename: The filename to save the state to
    :param recover: Recover the state if available

    '''
    from os.path import exists, abspath, join

    # Set some stuff
    self.state_filename = join(directory, state_filename)

    # Read state if available
    if recover == True and exists(state_filename):
      self.state = ApplicationState.load(state_filename)
      print "Recovered state from %s" % state_filename
      print self.get_history()
    else:
      def find_directory(working_directory):
        counter = 1
        while True:
          directory = join(working_directory, "dials-%d" % counter)
          if not exists(directory):
            return directory
          counter += 1
      self.state = ApplicationState(find_directory(abspath(directory)))

  def set_mode(self, mode):
    '''
    Set the current mode.

    :param mode: The mode to set

    '''
    # Is mode available?
    if mode not in self.mode_list:
      raise RuntimeError('Unknown mode: %s' % mode)

    # Set the mode
    self.state.mode = mode
    self.state.dump(self.state_filename)

  def get_mode(self):
    '''
    Get the current mode

    :return: The current mode

    '''
    return self.state.mode

  def set_parameters(self, parameters, short_syntax=False):
    '''
    Set the parameters.

    :param parameters: The parameters as a string
    :param show_syntax: Use command line string

    '''
    from libtbx.utils import Sorry
    self.state.parameters[self.get_mode()].set(parameters, short_syntax=short_syntax)
    self.state.dump(self.state_filename)

  def reset_parameters(self):
    '''
    Reset the parameters to the default values

    '''
    self.state.parameters[self.get_mode()].reset()
    self.state.dump(self.state_filename)

  def get_parameters(self, diff=True):
    '''
    Get the current parameters

    :param diff: Show only the modified parameters

    '''
    return self.state.parameters[self.get_mode()].get(diff=diff)

  def get_history(self):
    '''
    Get the history as a string

    :return: The history string

    '''
    return self.state.history()

  def goto(self, index):
    '''
    Change state to a different index

    :param index: The index to go to

    '''
    self.state.goto(index)
    self.state.dump(self.state_filename)

  def run(self):
    '''
    Run a program

    '''
    try:
      self.state.run()
      self.state.dump(self.state_filename)
    except Exception:
      self.state.dump(self.state_filename)
      raise



def print_error(exception):
  '''
  Print out the error message

  '''
  print ''
  print '*' * 80
  print 'USER ERROR: PLEASE REPLACE USER'
  print ''
  print exception
  print '*' * 80
  print ''


class Console(Cmd):
  '''
  A class to implement an interactive dials console

  '''

  # The default prompt
  prompt = ">> "

  def __init__(self):
    '''
    Initialise the console

    '''

    # Initialise the console base
    Cmd.__init__(self)

    # Create the controller object
    self.controller = Controller()

    # Set the prompt to show the current mode
    self.prompt = "%s >> " % self.controller.get_mode()

  def emptyline(self):
    ''' Do nothing on empty line '''
    pass

  def default(self, line):
    ''' The default line handler '''
    try:
      self.controller.set_parameters(line, short_syntax=True)
    except Exception:
      return Cmd.default(self, line)

  def do_mode(self, mode):
    ''' Set the program mode '''
    try:
      self.controller.set_mode(mode)
      self.prompt = "%s >> " % self.controller.get_mode()
    except Exception, e:
      print_error(e)

  def do_set(self, parameter):
    ''' Set a phil parameter '''
    try:
      self.controller.set_parameters(parameter, short_syntax=True)
    except Exception, e:
      print_error(e)

  def do_reset(self, line):
    ''' Reset parameters to default. '''
    try:
      self.controller.reset_parameters()
    except Exception, e:
      print_error(e)

  def do_load(self, filename):
    ''' Load a phil parameter file '''
    try:
      with open(filename) as infile:
        self.controller.set_parameters(infile.read())
    except Exception, e:
      print_error(e)

  def do_run(self, line):
    ''' Run a program '''
    try:
      self.controller.run()
      self.print_history()
    except Exception, e:
      print_error(e)

  def do_goto(self, line):
    ''' Goto a particular history state '''
    try:
      self.controller.goto(int(line))
      self.print_history()
    except Exception, e:
      print_error(e)

  def do_get(self, line):
    ''' Show all the possible parameters '''
    print self.controller.get_parameters(diff=True).as_str()

  def do_all(self, line):
    ''' Show all the possible parameters '''
    print self.controller.get_parameters(diff=False).as_str()

  def do_history(self, line):
    ''' Show the history. '''
    self.print_history()

  def do_import(self, params):
    ''' Imperative import command '''
    self.run_as_imperative("import", params)

  def do_find_spots(self, params):
    ''' Imperative find_spots command '''
    self.run_as_imperative("find_spots", params)

  def do_index(self, params):
    ''' Imperative index command '''
    self.run_as_imperative("index", params)

  def do_refine_bs(self, params):
    ''' Imperative refine_bs command '''
    self.run_as_imperative("refine_bs", params)

  def do_reindex(self, params):
    ''' Imperative reindex command '''
    self.run_as_imperative("reindex", params)

  def do_refine(self, params):
    ''' Imperative refine command '''
    self.run_as_imperative("refine", params)

  def do_integrate(self, params):
    ''' Imperative integrate command '''
    self.run_as_imperative("integrate", params)

  def do_export(self, params):
    ''' Imperative export command '''
    self.run_as_imperative("export", params)

  def do_exit(self, line):
    ''' Exit the console '''
    return True

  def do_EOF(self, line):
    ''' Exit the console '''
    print ''
    return True

  def run_as_imperative(self, mode, parameters):
    '''
    Helper for imperative mode. Change mode, set parameters and run the job

    '''
    try:
      self.controller.set_mode(mode)
      self.prompt = "%s >> " % self.controller.get_mode()
      self.controller.set_parameters(parameters, short_syntax=True)
      self.controller.run()
      self.print_history()
    except Exception, e:
      print_error(e)

  def print_history(self):
    '''
    Print the history

    '''
    print ''
    print 'History'
    print self.controller.get_history()

  def complete_mode(self, text, line, begidx, endidx):
    '''
    Offer tab completion options for changing mode.

    '''
    return [i for i in self.controller.mode_list if i.startswith(text)]


# The intro string for the console
CONSOLE_INTRO = '''
DIALS interactive mode
Type "help" for more information
'''


if __name__ == '__main__':

  # Print the console intro
  print CONSOLE_INTRO

  # Create the console
  console = Console()

  # Enter the command loop
  console.cmdloop()
