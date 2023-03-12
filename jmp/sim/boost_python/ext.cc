#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dials/algorithms/profile_model/gaussian_rs/coordinate_system.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error.h>
#include <dxtbx/model/detector.h>
#include <dxtbx/model/scan.h>
#include <random>

using namespace boost::python;

namespace dials {
namespace algorithms {
namespace boost_python {

using dials::algorithms::profile_model::gaussian_rs::CoordinateSystem;
using dxtbx::model::Detector;
using dxtbx::model::Panel;
using dxtbx::model::Scan;
using scitbx::vec2;
using scitbx::vec3;
using scitbx::af::int6;

af::versa<double, af::c_grid<3>>
compute_profile_internal(af::const_ref<double, af::c_grid<3>> grid, int6 bbox,
                         std::size_t zs, std::size_t ys, std::size_t xs,
                         std::size_t N, double delta_d, double delta_m,
                         const Detector &detector, const Scan &scan,
                         const CoordinateSystem &cs) {

  std::default_random_engine generator;
  std::uniform_real_distribution<double> uniform(0.0, 1.0);

  af::versa<double, af::c_grid<3>> profile(af::c_grid<3>(zs, ys, xs));

  double xoff = -delta_d;
  double yoff = -delta_d;
  double zoff = -delta_m;

  double xstep = (2 * delta_d) / grid.accessor()[2];
  double ystep = (2 * delta_d) / grid.accessor()[1];
  double zstep = (2 * delta_m) / grid.accessor()[0];

  for (std::size_t k = 0; k < grid.accessor()[0]; ++k) {
    for (std::size_t j = 0; j < grid.accessor()[1]; ++j) {
      for (std::size_t i = 0; i < grid.accessor()[2]; ++i) {
        double counts = grid(k, j, i);
        if (counts > 0) {
          double fraction = counts / double(N);
          for (std::size_t l = 0; l < N; ++l) {

            double gz = uniform(generator) + k;
            double gy = uniform(generator) + j;
            double gx = uniform(generator) + i;

            double e1 = xoff + gx * xstep;
            double e2 = yoff + gy * ystep;
            double e3 = zoff + gz * zstep;

            vec3<double> s1p = cs.to_beam_vector(vec2<double>(e1, e2));
            vec2<double> pxy = detector[0].get_ray_intersection_px(s1p);
            double px = pxy[0] - bbox[0];
            double py = pxy[1] - bbox[2];

            double phip = cs.to_rotation_angle_fast(e3);
            double pz = scan.get_array_index_from_angle(phip) - bbox[4];

            int kk = (int)std::floor(pz);
            int jj = (int)std::floor(py);
            int ii = (int)std::floor(px);

            /* std::cout << zs << ", " << cs.zeta() << ", " << cs.phi() << ", "
             * << phip << ", " << kk << std::endl; */

            /* std::cout << zs << ", " << ys << "," << xs << "->" << px << ", "
             * << py << " " << pz << std::endl; */

            if (kk >= 0 && kk < zs && jj >= 0 && jj < ys && ii >= 0 &&
                ii < xs) {
              profile(kk, jj, ii) += fraction;
            }
          }
        }
      }
    }
  }

  return profile;
}

BOOST_PYTHON_MODULE(dials_scratch_jmp_sim_ext) {
  def("compute_profile_internal", &compute_profile_internal);
}

} // namespace boost_python
} // namespace algorithms
} // namespace dials
