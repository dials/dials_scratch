#!/usr/bin/env cctbx.python

"""Examples of using and abusing scitbx::array_family types when writing
C++ extension functions and classes for use in Python"""

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

def demo_no_converter_for_const_ref():
  '''Demonstrate that a C++ function with return value scitbx::af::const_ref<>
  is not automatically wrapped to return a flex array in Python.'''

  from dials_scratch_cctbx_cpp_examples_ext import BadBucket
  from cctbx.array_family import flex

  print 'Create a C++ class that contains an array of Miller indices.'
  hkl = flex.miller_index([(0,0,i) for i in range(10)])
  bb = BadBucket(hkl)

  # cannot return an af::const_ref to Python
  print 'Attempt to access the contained array as an af::const_ref.'
  try:
    bb.get_const_ref_hkl()
  except TypeError as e:
    print "This fails with a TypeError and this message:"
    print e.message

  # however, can return an af::shared copy of the stored data
  print 'Now access the contained array as an af::shared.'
  hkl2 = bb.get_shared_hkl()
  print 'This works, and the result is a {0}'.format(type(hkl2))

  # check the initial and returned arrays are the same
  assert (hkl == hkl2).all_eq(True)

def demo_data_loss_with_const_ref_storage():
  '''Demonstrate that a C++ object that stores data as a
  scitbx::af::const_ref<> without taking a copy of the data cannot guarantee
  veracity of that data.'''

  from dials_scratch_cctbx_cpp_examples_ext import BadBucket
  from cctbx.array_family import flex

  print 'Create a C++ class that contains an array of Miller indices.'
  hkl = flex.miller_index([(0,0,i) for i in range(10)])
  bb = BadBucket(hkl)

  print 'Access the data from this object.'
  hkl2 = bb.get_shared_hkl()

  print 'Check all values are as expected. So far, so good.'
  assert (hkl2 == hkl).all_eq(True)

  print 'Now alter the original data.'
  hkl.fill((0,0,0))

  print 'Access data from the object again.'
  hkl3 = bb.get_shared_hkl()

  print ("See that altering the original array of Miller indices affects "
         "the data accessible by the class. Worse, if the original data goes "
         "out of scope or is deleted, the values returned from the C++ object "
         "will become garbage.")
  assert (hkl3 == flex.miller_index(10, (0,0,0))).all_eq(True)

  return

def demo_cannot_pass_versa_from_python():
  '''Demonstrate that multidimensional arrays using scitbx::af::versa type
  cannot be passed from Python by a setter taking that type. They can be
  passed as a scitbx::af::const_ref though.'''

  from dials_scratch_cctbx_cpp_examples_ext import TwoDimensionalArrayBroken
  from dials_scratch_cctbx_cpp_examples_ext import TwoDimensionalArrayFixed
  from scitbx.array_family import flex

  print 'Create a 2D flex.double array'
  nrow = 10
  ncol = 2
  array = flex.double(range(nrow * ncol))
  array.reshape(flex.grid(nrow, ncol))

  broken = TwoDimensionalArrayBroken()
  fixed = TwoDimensionalArrayFixed()

  # Cannot pass array from Python as versa
  print 'Attempt to pass the array to a C++ extension as an af::versa.'
  try:
    broken.set_array_data(array)
  except TypeError as e:
    print "This fails with a Boost.Python.ArgumentError and this message:"
    print e.message

  # Can pass the array from Python as const_ref
  print 'Now attempt to pass the array to a C++ extension as an af::const_ref.'
  fixed.set_array_data(array)

  # Getter as versa works fine
  array2 = fixed.get_array_data()

  assert (array2 == array).all_eq(True)
  print 'This works, and the data can be passed back to Python via af::versa'


if __name__ == '__main__':

  runner(demo_no_converter_for_const_ref)
  runner(demo_data_loss_with_const_ref_storage)
  runner(demo_cannot_pass_versa_from_python)
