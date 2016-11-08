/*
 * array_gotchas.h
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

namespace dials_scratch { namespace examples {

  void print_array(scitbx::af::shared< double > input) {
    std::cout << input.const_ref() << "\n";
  }

}} // namespace dials_scratch::examples

#endif // DIALS_SCRATCH_EXAMPLES_DEBUG_TOOLS_H
