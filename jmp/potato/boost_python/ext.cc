#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <scitbx/constants.h>
#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/mat2.h>
#include <scitbx/mat3.h>
#include <scitbx/math/r3_rotation.h>
#include <scitbx/matrix/multiply.h>
#include <dials/error.h>

using namespace boost::python;

namespace dials { namespace algorithms { namespace boost_python {

  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::mat2;
  using scitbx::mat3;

  /**
   * Function to given chisq quantile value
   * @param k The degrees of freedom
   * @param p The probability
   */
  double chisq_quantile(int k, double p) {
    DIALS_ASSERT(k > 0);
    DIALS_ASSERT(p >= 0 && p <= 1);
    boost::math::chi_squared_distribution<> dist(k);
    return boost::math::quantile(dist, p);
  }


  BOOST_PYTHON_MODULE(dials_scratch_jmp_potato_ext)
  {
    def("chisq_quantile", &chisq_quantile);
  }

}}}
