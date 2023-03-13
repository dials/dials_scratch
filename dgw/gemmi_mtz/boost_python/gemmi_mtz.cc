#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dials_scratch/dgw/gemmi_mtz/gemmi_mtz.h>

using namespace boost::python;
namespace dials {
namespace gemmi_mtz {
namespace boost_python {

void export_create_mtz() {
  def("create_mtz", &create_mtz, (arg("title"), arg("reflections")));
}

BOOST_PYTHON_MODULE(dials_gemmi_mtz_ext) { export_create_mtz(); }

} // namespace boost_python
} // namespace gemmi_mtz
} // namespace dials
