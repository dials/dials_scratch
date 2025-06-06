#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dials/error.h>

using namespace boost::python;

namespace dials_scratch {
namespace examples {
namespace boost_python {

void export_bad_bucket();
void export_two_dimensional_array();
void export_print_array();
void export_print_array_head();
void export_create_sparse_matrix();
void export_mat_sum();
void export_vector_of_arrays();

BOOST_PYTHON_MODULE(dials_scratch_cctbx_cpp_examples_ext) {
  export_bad_bucket();
  export_two_dimensional_array();
  export_print_array();
  export_print_array_head();
  export_create_sparse_matrix();
  export_mat_sum();
  export_vector_of_arrays();
}

} // namespace boost_python
} // namespace examples
} // namespace dials_scratch
