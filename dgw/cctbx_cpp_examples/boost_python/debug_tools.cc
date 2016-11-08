/*
 * array_gotchas.cc
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
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/ref.h>
#include <dials_scratch/dgw/cctbx_cpp_examples/debug_tools.h>

namespace dials_scratch { namespace examples { namespace boost_python {

  using namespace boost::python;

  void export_print_array()
  {
    def("print_array", &print_array, (
      arg("input")));
  }

}}} // namespace = dials_scratch::examples::boost_python

