
#ifndef DIALS_STILLS_SIMULATOR_H
#define DIALS_STILLS_SIMULATOR_H

#include <cmath>
#include <scitbx/constants.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/array_family/reflection_table.h>
#include <dxtbx/model/detector.h>
#include <dxtbx/model/scan.h>
#include <dxtbx/model/crystal.h>
#include <dials/algorithms/spot_prediction/pixel_to_miller_index.h>
#include <dials/algorithms/background/glm/robust_poisson_mean.h>
#include <dials/algorithms/integration/sum/summation.h>
#include <dials/model/data/mask_code.h>
#include <dials/model/data/shoebox.h>
#include <dials_scratch/jmp/stills/simulator.h>
#include <dials/error.h>


namespace dials {

  using scitbx::constants::pi;
  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::mat3;
  using dxtbx::model::Beam;
  using dxtbx::model::Detector;
  using dxtbx::model::Crystal;


  class StillsSimulator {
  public:

    StillsSimulator(
          const Beam &beam,
          const Detector &detector,
          const Crystal &crystal,
          mat3<double> sigma_wavelength_spread,
          mat3<double> sigma_angular_spread,
          mat3<double> sigma_rlp_mosaicity,
          std::size_t num_points)
      : beam_(beam),
        detector_(detector),
        crystal_(crystal),
        sigma_wavelength_spread_(sigma_wavelength_spread),
        sigma_angular_spread_(sigma_angular_spread),
        sigma_rlp_mosaicity_(sigma_rlp_mosaicity),
        num_points_(num_points) {
    }

    double normal_pdf(vec3<double> x, vec3<double> mu, mat3<double> sigma) const {
      mat3<double> sigma_inv = sigma.inverse();
      double A = 1.0 / (std::sqrt(std::pow(2*pi,3) *sigma.determinant()));
      double B = (x-mu)*sigma_inv*(x-mu);
      return A * std::exp(-0.5 * B);
    }

    double polygon_area(vec3<double> A, vec3<double> B, vec3<double> C, vec3<double> D) const {
      return ((A-B).cross(A-D)).length() / 2.0 + ((C-B).cross(C-D)).length() / 2.0;
    }

    double pixel_area(std::size_t panel, double s0_length, int x, int y) const {

      // Get beam vectors at each corner
      vec3<double> s00 = detector_[panel].get_pixel_lab_coord(
          vec2<double>(x,y)).normalize() * s0_length;
      vec3<double> s01 = detector_[panel].get_pixel_lab_coord(
          vec2<double>(x+1,y)).normalize() * s0_length;
      vec3<double> s10 = detector_[panel].get_pixel_lab_coord(
          vec2<double>(x,y+1)).normalize() * s0_length;
      vec3<double> s11 = detector_[panel].get_pixel_lab_coord(
          vec2<double>(x+1,y+1)).normalize() * s0_length;

      // Compute the pixel area
      return polygon_area(s00, s01, s11, s10);
    }

    double integrate_pixel(cctbx::miller::index<> h, std::size_t panel, int x, int y) const {
    
      // Get some stuff from the models
      mat3<double> A = crystal_.get_A();
      mat3<double> U = crystal_.get_U();
      vec3<double> s0 = beam_.get_s0();

      // The reciprocal lattice vector
      vec3<double> rlp = A * h;

      // Construct the rotated covariance matrix
      mat3<double> sigma_M = U * sigma_rlp_mosaicity_ * U.transpose();
      
      // Normalize the ewald sphere and boost the rlp and sigma
      double wavelength = 1.0 / s0.length();
      rlp = wavelength * rlp;
      s0 = s0.normalize();
      sigma_M = wavelength * wavelength * sigma_M;

      // Construct the full covariance matrix
      mat3<double> sigma = sigma_M;
      if (rlp.length() > 0) {

        // Construct the convolution. Both components scale with the length of the
        // reciprocal lattice vector so multiply this here.
        mat3<double> sigma_lw = (sigma_wavelength_spread_ + sigma_angular_spread_) * rlp.length();
    
        // The coordinate system at the end of the rlp
        vec3<double> s2 = s0 + rlp;
        vec3<double> e1 = s2.cross(s0).normalize();
        vec3<double> e2 = -e1.cross(rlp).normalize();
        vec3<double> e3 = rlp.normalize();
    
        // The change of basis matrix
        mat3<double> E(
            e1[0], e2[0], e3[0],
            e1[1], e2[1], e3[1],
            e1[2], e2[2], e3[2]);
        /* DIALS_ASSERT(E.is_r3_rotation_matrix()); */
      
        // Construct the rotated covariance matrices
        mat3<double> sigma_E = E * sigma_lw * E.transpose();
        sigma += sigma_E;

      }

      // Compute the pixel area
      double s0_length = s0.length();
      double area = pixel_area(panel, s0_length, x, y);

      // Do the integration over the pixel
      double N = num_points_;
      double I = 0;
      for (std::size_t j = 0; j < N; ++j) {
        for (std::size_t i = 0; i < N; ++i) {
          double xx = x + (i + 0.5) / N;
          double yy = y + (j + 0.5) / N;
          vec3<double> ss = detector_[panel].get_pixel_lab_coord(
              vec2<double>(xx, yy)).normalize() * s0_length;
          I += normal_pdf(ss, s0+rlp, sigma);
        }
      }

      return area * I / (N * N);
    }

  protected:

    Beam beam_;
    Detector detector_;
    Crystal crystal_;
    mat3<double> sigma_wavelength_spread_;
    mat3<double> sigma_angular_spread_;
    mat3<double> sigma_rlp_mosaicity_;
    std::size_t num_points_;
  };

}


#endif
