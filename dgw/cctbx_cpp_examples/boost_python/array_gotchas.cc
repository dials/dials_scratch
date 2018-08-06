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
#include <cctbx/miller.h>
#include <dials_scratch/dgw/cctbx_cpp_examples/array_gotchas.h>

namespace dials_scratch { namespace examples { namespace boost_python {

  using namespace boost::python;

  typedef cctbx::miller::index<> miller_index;

  void export_bad_bucket()
  {
    // Create and return the wrapper
    class_ <BadBucket> ("BadBucket", no_init)
      .def(init< scitbx::af::const_ref< miller_index > >((
        arg("hkl"))))
      .def("get_const_ref_hkl", &BadBucket::get_const_ref_hkl)
      .def("get_shared_hkl", &BadBucket::get_shared_hkl)
      ;
  }

  void export_two_dimensional_array()
  {
    // Broken version
    class_ <TwoDimensionalArray> ("TwoDimensionalArrayBroken")
      .def("set_array_data", &TwoDimensionalArray::set_array_data_from_versa)
      .def("get_array_data", &TwoDimensionalArray::get_array_data)
      ;

    // Fixed version
    class_ <TwoDimensionalArray> ("TwoDimensionalArrayFixed")
      .def("set_array_data", &TwoDimensionalArray::set_array_data_from_const_ref)
      .def("get_array_data", &TwoDimensionalArray::get_array_data)
      ;
  }

}}} // namespace = dials_scratch::examples::boost_python

