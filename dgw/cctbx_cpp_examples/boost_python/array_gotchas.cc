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
#include <cctbx/miller.h>
#include <dials_scratch/dgw/cctbx_cpp_examples/array_gotchas.h>
#include <scitbx/array_family/ref.h>
#include <scitbx/array_family/shared.h>

namespace dials_scratch {
namespace examples {
namespace boost_python {

using namespace boost::python;

typedef cctbx::miller::index<> miller_index;

void export_bad_bucket() {
  // Create and return the wrapper
  class_<BadBucket>("BadBucket", no_init)
      .def(init<scitbx::af::const_ref<miller_index>>((arg("hkl"))))
      .def("get_const_ref_hkl", &BadBucket::get_const_ref_hkl)
      .def("get_shared_hkl", &BadBucket::get_shared_hkl);
}

void export_two_dimensional_array() {
  // Broken version
  class_<TwoDimensionalArray>("TwoDimensionalArrayBroken")
      .def("set_array_data", &TwoDimensionalArray::set_array_data_from_versa)
      .def("get_array_data", &TwoDimensionalArray::get_array_data);

  // Fixed version
  class_<TwoDimensionalArray>("TwoDimensionalArrayFixed")
      .def("set_array_data",
           &TwoDimensionalArray::set_array_data_from_const_ref)
      .def("get_array_data", &TwoDimensionalArray::get_array_data);
}

void export_vector_of_arrays() {
  // Broken version (segfaults)
  class_<BadVectorOfArrays>("BadVectorOfArrays")
      .def("add_array_to_vector", &BadVectorOfArrays::add_array_to_vector)
      .def("get_sum", &BadVectorOfArrays::get_sum);

  // Working version
  class_<VectorOfArrays>("VectorOfArrays")
      .def("add_array_to_vector", &VectorOfArrays::add_array_to_vector)
      .def("get_sum", &VectorOfArrays::get_sum);
}

} // namespace boost_python
} // namespace examples
} // namespace dials_scratch
