#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dials/error.h>

using namespace boost::python;

namespace dials_scratch { namespace examples { namespace boost_python {

  void export_bad_bucket();

  BOOST_PYTHON_MODULE(dials_scratch_cctbx_cpp_examples_ext)
  {
    export_bad_bucket();
  }

}}} // namespace dials_scratch::examples::boost_python
