/*
 * debug_tools.cc
 *
 *  Copyright (C) (2016) STFC Rutherford Appleton Laboratory, UK.
 *
 *  Author: David Waterman
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dials_scratch/dgw/cctbx_cpp_examples/debug_tools.h>
#include <scitbx/array_family/ref.h>
#include <scitbx/array_family/shared.h>

namespace dials_scratch {
namespace examples {
namespace boost_python {

using namespace boost::python;

void export_print_array() { def("print_array", &print_array, (arg("input"))); }

// http://www.boost.org/doc/libs/1_62_0/libs/python/doc/html/tutorial/tutorial/functions.html#tutorial.functions.default_arguments
// The default argument value defined in the C++ function is not automatically
// picked up, but this macro will produce the right thin wrappers to get
// around this. The values 1 and 2 refer to the minimum and maximum number
// of arguments to the function. The function print_array_head_overloads
// created by the macro is then passed to def, below.
BOOST_PYTHON_FUNCTION_OVERLOADS(print_array_head_overloads, print_array_head, 1,
                                2)
void export_print_array_head() {

  def("print_array_head", print_array_head, print_array_head_overloads());
}

} // namespace boost_python
} // namespace examples
} // namespace dials_scratch
