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
#include <dxtbx/model/experiment.h>
#include <dials/model/data/shoebox.h>
#include <dials/algorithms/profile_model/gaussian_rs/coordinate_system.h>
#include <dials/algorithms/shoebox/mask_code.h>
#include <dials/array_family/reflection_table.h>
#include <dials/error.h>

using namespace boost::python;

namespace dials { namespace algorithms { namespace boost_python {

  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::mat2;
  using scitbx::mat3;
  using scitbx::matrix::transpose_multiply;
  using scitbx::matrix::multiply_transpose;
  using scitbx::af::int6;
  using dxtbx::model::Experiment;
  using dxtbx::model::BeamBase;
  using dxtbx::model::Panel;
  using dxtbx::model::Detector;
  using dials::model::Shoebox;
  using dials::model::Valid;
  using dials::model::Foreground;
  using dials::model::Background;
  using dials::model::Overlapped;
  using dials::algorithms::profile_model::gaussian_rs::CoordinateSystem2d;

  namespace detail {

    double AT_B_A(vec2<double> A, mat2<double> B) {
      vec2<double> ATB = A * B;
      return ATB*A;
    }

  }


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

  mat3<double> compute_change_of_basis_operation(vec3<double> s0, vec3<double> s2) {
    const double TINY = 1e-7;
    DIALS_ASSERT((s2 - s0).length() > TINY);
    vec3<double> e1 = s2.cross(s0).normalize();
    vec3<double> e2 = s2.cross(e1).normalize();
    vec3<double> e3 = s2.normalize();
    mat3<double> R(
        e1[0], e1[1], e1[2],
        e2[0], e2[1], e2[2],
        e3[0], e3[1], e3[2]);
    return R;
  }

  class MaskCalculator {
  public:

    MaskCalculator(Experiment experiment, mat3<double> sigma)
      : experiment_(experiment),
        sigma_(sigma) {
    }

    void compute(af::reflection_table reflections) const {

      // Get some array from the reflection table
      af::const_ref< vec3<double> > s1 = reflections["s1"];
      af::const_ref< vec3<double> > s2 = reflections["s2"];
      af::ref< Shoebox<> > sbox = reflections["shoebox"];

      // Compute quantile
      double D = chisq_quantile(2, 0.997);

      // Compute mask for all reflections
      for (std::size_t i = 0; i < reflections.size(); ++i) {
        compute_single(s1[i], s2[i], sbox[i], D);
      }
    }

  protected:

    void compute_single(
            vec3<double> s1,
            vec3<double> s2,
            Shoebox<> sbox,
            double D) const {

      const double TINY = 1e-7;
      DIALS_ASSERT(s1.length() > 0);
      DIALS_ASSERT(s2.length() > 0);
      DIALS_ASSERT(sbox.is_consistent());
      DIALS_ASSERT(D > 0);

      // Get the beam and detector
      DIALS_ASSERT(experiment_.get_beam() != NULL);
      DIALS_ASSERT(experiment_.get_detector() != NULL);
      Detector detector = *experiment_.get_detector();

      // The the indicident beam vector
      vec3<double> s0 = experiment_.get_beam()->get_s0();
      double s0_length = s0.length();
      DIALS_ASSERT(std::abs(s0_length - s1.length()) < TINY);

      // Compute the change of basis for the reflection
      mat3<double> R = compute_change_of_basis_operation(s0, s2);

      // Rotate the covariance matrix and s2 vector
      mat3<double> S = R*sigma_*R.transpose();
      vec3<double> mu = R*s2;
      vec3<double> zaxis(0,0,1);
      DIALS_ASSERT(std::abs((mu.normalize() * zaxis) - 1) < TINY);

      // Partition the covariance matrix
      mat2<double> S11(S[0], S[1], S[3], S[4]);
      vec2<double> S12(S[2], S[5]);
      vec2<double> S21(S[6], S[7]);
      double S22 = S[8];

      // Partition the mean vector
      vec2<double> mu1(mu[0], mu[1]);
      double mu2 = mu[2];

      // Compute epsilon the distance to the Ewald sphere
      double epsilon = s0.length() - mu2;

      // Compute the mean of the conditional distribution
      DIALS_ASSERT(S22 > 0);
      double S22_inv = 1.0 / S22;
      vec2<double> mubar = mu1 + S12*S22_inv*epsilon;

      // Compute the covariance of the conditional distribution
      mat2<double> S12_S21;
      multiply_transpose(&S12[0], &S21[0], 2, 1, 2, &S12_S21[0]);
      mat2<double> Sbar = S11 - S12_S21 * S22_inv;
      mat2<double> Sbar_inv = Sbar.inverse();

      // Get the mask array
      af::ref< int, af::c_grid<3> > mask = sbox.mask.ref();

      // Get the bounding box
      int x0 = sbox.bbox[0];
      int x1 = sbox.bbox[1];
      int y0 = sbox.bbox[2];
      int y1 = sbox.bbox[3];
      int z0 = sbox.bbox[4];
      int z1 = sbox.bbox[5];
      DIALS_ASSERT(x0 < x1);
      DIALS_ASSERT(y0 < y1);
      DIALS_ASSERT(z0 < z1);
      DIALS_ASSERT(z1 - z0 == 1);

      // Create the coordinate system
      CoordinateSystem2d cs(s0, s2);

      // Get the panel model
      Panel panel = detector[sbox.panel];

      // Set the mask value for each pixel
      DIALS_ASSERT(mask.accessor()[0] == 1);
      for (std::size_t j = 0; j < mask.accessor()[1]; ++j) {
        for (std::size_t i = 0; i < mask.accessor()[2]; ++i) {
          int ii = x0 + ((int)i);
          int jj = y0 + ((int)j);

          // The pixel coordinates of the corners
          vec2<double> p1(ii,jj);
          vec2<double> p2(ii+1,jj);
          vec2<double> p3(ii,jj+1);
          vec2<double> p4(ii+1,jj+1);

          // The lab coordinates of the pixel corners
          vec3<double> sp1 = panel.get_pixel_lab_coord(p1).normalize() * s0_length;
          vec3<double> sp2 = panel.get_pixel_lab_coord(p2).normalize() * s0_length;
          vec3<double> sp3 = panel.get_pixel_lab_coord(p3).normalize() * s0_length;
          vec3<double> sp4 = panel.get_pixel_lab_coord(p4).normalize() * s0_length;

          // The coordinates in kabsch space
          vec2<double> x1 = cs.from_beam_vector(sp1);
          vec2<double> x2 = cs.from_beam_vector(sp2);
          vec2<double> x3 = cs.from_beam_vector(sp3);
          vec2<double> x4 = cs.from_beam_vector(sp4);

          // The distance from the mean
          double d1 = detail::AT_B_A(x1-mubar, Sbar_inv);
          double d2 = detail::AT_B_A(x2-mubar, Sbar_inv);
          double d3 = detail::AT_B_A(x3-mubar, Sbar_inv);
          double d4 = detail::AT_B_A(x4-mubar, Sbar_inv);

          // The minimum distance
          if (std::min(std::min(d1, d2), std::min(d3, d4)) < D) {
            mask(0,j,i) |= Foreground;
          } else {
            mask(0,j,i) |= Background;
          }
        }
      }
    }

    Experiment experiment_;
    mat3<double> sigma_;
  };

  BOOST_PYTHON_MODULE(dials_scratch_jmp_potato_ext)
  {
    def("chisq_quantile", &chisq_quantile);


    class_<MaskCalculator>("MaskCalculator", no_init)
      .def(init<Experiment, mat3<double> >())
      .def("compute", &MaskCalculator::compute)
      ;
  }

}}}
