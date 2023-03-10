#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dials_scratch/dgw/gemmi_mtz/gemmi_mtz.h>

namespace dials_scratch { namespace gemmi_mtz { namespace boost_python {

  using namespace boost::python;

  void export_create_mtz()
  {
    def("create_mtz", &create_mtz);
  }

}}} // namespace = dials_scratch::gemmi_mtz::boost_python

