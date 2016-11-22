/*
 * debug_tools.h
 *
 *  Copyright (C) (2016) STFC Rutherford Appleton Laboratory, UK.
 *
 *  Author: David Waterman
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */

#ifndef DIALS_SCRATCH_EXAMPLES_DEBUG_TOOLS_H
#define DIALS_SCRATCH_EXAMPLES_DEBUG_TOOLS_H

#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/ref.h>
#include <dials/error.h>
#include <scitbx/array_family/simple_io.h> // for std::cout overloads
#include <algorithm> // for std::min

namespace dials_scratch { namespace examples {

  void print_array(scitbx::af::shared< double > input) {
    std::cout << input.const_ref() << "\n";
  }

  // Demonstrate 'slicing' in C++ and default argument values
  void print_array_head(scitbx::af::shared< double > input,
                        std::size_t head_len = 10) {
    // The default value of 10 for head_len is not automatically picked up
    // by the Boost.Python wrapper. See debug_tools.cc for use of the macro
    // BOOST_PYTHON_FUNCTION_OVERLOADS to create the right wrappers.

    // ensure the head length does not exceed the array size
    head_len = std::min(head_len, input.size());

    // array 'slicing' using the iterator constructor (JMP)
    scitbx::af::shared<double> head(
      &input[0],
      &input[0] + head_len
    );

    // now the same std::cout overload as print_array
    std::cout << head.const_ref() << "\n";
  }

}} // namespace dials_scratch::examples

#endif // DIALS_SCRATCH_EXAMPLES_DEBUG_TOOLS_H
