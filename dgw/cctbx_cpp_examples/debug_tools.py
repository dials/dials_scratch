#!/usr/bin/env cctbx.python

"""Examples of debugging tools"""

import textwrap

run_number = 1
def runner(fn):
  global run_number
  title = 'Example {0}'.format(run_number)
  print title
  print '-' * len(title)
  try:
    desc = textwrap.fill(' '.join(fn.__doc__.split()), width=50,
                            initial_indent='  ', subsequent_indent='  ')
  except AttributeError:
    desc = 'No description'
  print desc
  print
  fn()
  print
  run_number += 1

def demo_print_array():
  '''Demonstrate simple array printing from C++ using simple_io.h.'''

  from dials_scratch_cctbx_cpp_examples_ext import print_array
  from cctbx.array_family import flex

  print 'Create a flex.double array'
  v = flex.double([0,1,2,3,4,5,6,7,8,9])

  print 'Print from Python:'
  print list(v)

  print 'Print from C++:'
  print_array(v)

  return

def demo_print_array_head():
  '''Demonstrate simple array slicing from C++'''

  from dials_scratch_cctbx_cpp_examples_ext import print_array_head
  from cctbx.array_family import flex
  v = flex.double(range(100))

  print "Print array head, default head length of 10:"
  print_array_head(v, 10)

  print "Print array head, set head length to 5"
  print_array_head(v, 5)

if __name__ == '__main__':

  runner(demo_print_array)
  runner(demo_print_array_head)


