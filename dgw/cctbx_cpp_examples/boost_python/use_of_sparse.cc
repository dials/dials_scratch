/*
 * use_of_sparse.cc
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
#include <dials_scratch/dgw/cctbx_cpp_examples/use_of_sparse.h>

namespace dials_scratch {
namespace examples {
namespace boost_python {

using namespace boost::python;

using scitbx::sparse::matrix;

void export_create_sparse_matrix() {
  // note this needs "from scitbx import sparse" before use in Python
  def("create_sparse_matrix", &create_sparse_matrix);
}

void export_mat_sum() { def("mat_sum", &mat_sum); }

} // namespace boost_python
} // namespace examples
} // namespace dials_scratch
