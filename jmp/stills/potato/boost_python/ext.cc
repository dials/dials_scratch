#include <boost/python.hpp>
#include <boost/python/def.hpp>
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

  namespace detail {

    /**
     * The conditional distribution of a multivariate normal P(x,y|z=z0)
     */
    class ConditionalDistributionAtZ {
    public:

      /**
       * Compute the conditional distribution
       * @param mu The mean vector
       * @param sigma The covariance matrix
       * @param z The value of z
       */
      ConditionalDistributionAtZ(vec3<double> mu, mat3<double> sigma, double z) {

        using scitbx::matrix::multiply_transpose;

        // Partition the covariance matrix
        // | | S11 S12 | | S13 | |
        // | | S21 S22 | | S23 | |
        // |   -------     ---   |
        // | | S31 S32 | | S33 | |
        mat2<double> sigma11(
            sigma[0], sigma[1],
            sigma[3], sigma[4]);
        vec2<double> sigma12(sigma[2], sigma[5]);
        vec2<double> sigma21(sigma[6], sigma[7]);
        double sigma22 = sigma[8];

        // Partition the mean vector
        // | | mu1 | |
        // | | mu2 | |
        // |   --- | |
        // | | mu3 | |
        vec2<double> mu1(mu[0], mu[1]);
        double mu2 = mu[2];

        // Compute the mean of the conditional distribution
        mu_ = mu1 + sigma12 * (1/sigma22) * (z - mu2);

        // Compute the covariance matrix of the conditional distribution
        mat2<double> A;
        multiply_transpose(&sigma12[0], &sigma21[0], 2, 1, 2, &A[0]);
        sigma_ = sigma11 - A;

/*         if (sigma_.determinant() <= 0) { */
/*           std::cout << "SIGMA: "; */
/*           for (std::size_t i = 0; i < 9; ++i) { */
/*             std::cout << sigma[i] << ", "; */
/*           } */
/*           std::cout << std::endl; */
          
/*           std::cout << "MU: "; */
/*           for (std::size_t i = 0; i < 3; ++i) { */
/*             std::cout << mu[i] << ", "; */
/*           } */
/*           std::cout << std::endl; */
          
/*           std::cout << "SIGMA': "; */
/*           for (std::size_t i = 0; i < 4; ++i) { */
/*             std::cout << sigma_[i] << ", "; */
/*           } */
/*           std::cout << std::endl; */

/*           std::cout << "Z: " << z << std::endl; */
/*           std::cout << "DET: " << sigma_.determinant() << std::endl; */
/*           DIALS_ASSERT(false); */
/*         } */
      }

      /**
       * @returns The mean of the conditional distribution
       */
      vec2<double> mu() const {
        return mu_;
      }

      /**
       * @returns The covariance of the conditional distribution
       */
      mat2<double> sigma() const {
        return sigma_;
      }

    protected:

      vec2<double> mu_;
      mat2<double> sigma_;
    };

  }


  /**
   * A class to compute things associated with the multivariate normal on the
   * Ewald sphere surface.
   */
  class PotatoOnEwaldSphere {
  public:

    PotatoOnEwaldSphere(double wavelength, vec3<double> mu, mat3<double> sigma)
      : mu_(mu),
        sigma_(sigma),
        R_(compute_rotation_of_mu_onto_z(mu)),
        mup_(R_.transpose() * mu),
        sigmap_(R_.transpose() * sigma_ * R_) {
      DIALS_ASSERT(wavelength > 0);
      radius_ = 1.0 / wavelength;
    }

    /**
     * Return the likelihood of the point at this distance from the Ewald sphere
     * Note that here sigma comes from the covariance matrix so is the variance
     * not the standard deviation.
     */
    double log_likelihood() const {
      double mu = mup_[2];
      double sigma = sigmap_[8];
      return -0.5 * (std::log(sigma) + (1/sigma) * (radius_-mu)*(radius_-mu));
    }

    double conditional_likelihood(double x, double y) const {
      
      detail::ConditionalDistributionAtZ distribution(mup_, sigmap_, radius_);
      vec2<double> v(x, y);
      vec2<double> mu = distribution.mu();
      mat2<double> sigma = distribution.sigma();
      double D = sigma.determinant();
      if (D <= 0) {
        std::cout << "NOO" << std::endl;
        return -1e10;
      }
      return -0.5 * (std::log(sigma.determinant()) + (x-mu)*sigma.inverse()*(x-mu));
    }

    /**
     * Compute the scale factor for the reflection. This is calculated by taking
     * the marginal probability along mu, P(z) and then calculating the ratio
     * P(z=z0) / P(z=mu_z) where z0 is the radius of the ewald sphere.
     */
    double scale_factor() const {
      double mu = mup_[2];
      double sigma = sigmap_[8];
      return std::exp(-0.5*(radius_-mu) * (1/sigma) * (radius_-mu));
    }

    /**
     * Compute the centre of mass of the reflection on the Ewald sphere. This is
     * calculated by taking the conditional distribution on the tangent plane at
     * the ewald sphere surface and finding the mean of the marginal
     * distribution
     */
    vec3<double> centre_of_mass_on_ewald_sphere() const {

      // The conditional distribution at the Ewald sphere
      detail::ConditionalDistributionAtZ distribution(mup_, sigmap_, radius_);

      // Compute the vector to the centre of mass
      vec3<double> v(
          distribution.mu()[0],
          distribution.mu()[1],
          radius_);

      // Rotate the normalized vector back to the original coodinate system
      return R_ * (v.normalize() * radius_);
    }

    vec2<double> conditional_mean() const {
      return detail::ConditionalDistributionAtZ(mup_, sigmap_, radius_).mu();
    }
   
    mat2<double> conditional_sigma() const {
      return detail::ConditionalDistributionAtZ(mup_, sigmap_, radius_).sigma();
    }

  protected:

    /**
     * Compute the rotation matrix between the vector mu and the z axis
     * @param mu The mean offset vector from the centre of the Ewald sphere
     * @returns The rotation matrix
     */
    mat3<double> compute_rotation_of_mu_onto_z(vec3<double> mu) const {

      using scitbx::math::r3_rotation::axis_and_angle_as_matrix;

      const double TINY = 1e-7;

      // Compute the axes to aline
      vec3<double> m_axis = mu.normalize();
      vec3<double> z_axis(0, 0, 1);

      // Compute the rotation angle between z and mu. If the angle is very small
      // then just return the identity matrix.
      double angle = z_axis.angle_rad(m_axis).get();
      if (angle < TINY) {
        return mat3<double>(
            1, 0, 0,
            0, 1, 0,
            0, 0, 1);
      }

      // Compute the rotation axis
      vec3<double> r_axis = z_axis.cross(m_axis).normalize();

      // Return the rotation matrix
      return axis_and_angle_as_matrix(r_axis, angle, false);
    }

    vec3<double> mu_;
    mat3<double> sigma_;
    mat3<double> R_;
    vec3<double> mup_;
    mat3<double> sigmap_;
    double radius_;
  };



  BOOST_PYTHON_MODULE(dials_scratch_jmp_stills_potato_ext)
  {
    class_<PotatoOnEwaldSphere>("PotatoOnEwaldSphere", no_init)
      .def(init<
          double,
          vec3<double>,
          mat3<double> >())
      .def("scale_factor",
          &PotatoOnEwaldSphere::scale_factor)
      .def("centre_of_mass_on_ewald_sphere",
          &PotatoOnEwaldSphere::centre_of_mass_on_ewald_sphere)
      .def("log_likelihood",
          &PotatoOnEwaldSphere::log_likelihood)
      .def("conditional_mean",
          &PotatoOnEwaldSphere::conditional_mean)
      .def("conditional_sigma",
          &PotatoOnEwaldSphere::conditional_sigma)
      .def("conditional_likelihood",
          &PotatoOnEwaldSphere::conditional_likelihood)
      ;
  }

}}}
