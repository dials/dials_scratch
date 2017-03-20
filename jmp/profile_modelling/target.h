/*
 * target.h
 *
 *  Copyright (C) 2017 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_PROFILE_MODELLING_TARGET_H
#define DIALS_ALGORITHMS_PROFILE_MODELLING_TARGET_H

#include <random>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <scitbx/array_family/tiny_types.h>
#include <scitbx/array_family/simple_tiny_io.h>
#include <scitbx/mat3.h>
#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/matrix/cholesky.h>
#include <scitbx/math/r3_rotation.h>
#include <scitbx/constants.h>
#include <dxtbx/model/beam.h>
#include <dxtbx/model/detector.h>
#include <dxtbx/model/goniometer.h>
#include <dxtbx/model/scan.h>
#include <dxtbx/model/experiment.h>
#include <dials/array_family/reflection_table.h>
#include <dials/model/data/mask_code.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using scitbx::mat3;
  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::af::double6;
  using scitbx::af::int6;
  using dxtbx::model::Beam;
  using dxtbx::model::Panel;
  using dxtbx::model::Goniometer;
  using dxtbx::model::Scan;
  using dxtbx::model::Experiment;
  using dials::af::reflection_table;
  using dials::model::Shoebox;
  using dials::model::Valid;
  using dials::model::Foreground;

  namespace detail {
 
    /**
     * Helper function to select 2 orthogonal vectors
     * @param x The input vector
     * @returns A vector orthogonal to x
     */
    inline
    vec3<double> select_orthogonal(vec3<double> x) {
      std::size_t m = 0;
      for (std::size_t i = 1; i < 3; ++i) {
        if (x[i] > x[m]) {
          m = i;
        }
      }
      std::size_t n = (m+1) % 3;
      vec3<double> y(0,0,0);
      y[n] = x[m];
      y[m] = -x[n];
      DIALS_ASSERT(x * y < 1e-7);
      DIALS_ASSERT(y.length() > 0);
      return y.normalize();
    }

    /**
     * Check if an angle is "small"
     * @param angle The angle
     * @returns Is the angle "small"
     */
    inline
    bool small_angle(double angle) {
      return std::abs(angle) < scitbx::deg_as_rad(10.0);
    }

  }

  /**
   * A class to model the reciprocal lattice point spread
   */
  class ReciprocalLatticePointSpread {

  public:
    
    ReciprocalLatticePointSpread(mat3<double> A)
      : A_(A) {
      set_parameters(double6(0,0,0,0,0,0));
    }

    ReciprocalLatticePointSpread(mat3<double> A, double6 parameters)
      : A_(A) {
      set_parameters(parameters);
    }

    mat3<double> get_A() const {
      return A_;
    }

    double6 get_parameters() const {
      return parameters_;
    }

    void set_parameters(double6 parameters) {
      // Set parameters
      parameters_ = parameters;

      // The parameters are given as 6 parameters in the lower triangular matrix
      // form and the covariance matrix is given by the Cholesky decomposition
      // COV = LL*
      mat3<double> L(
          parameters[0], 0.0, 0.0,
          parameters[1], parameters[2], 0.0,
          parameters[3], parameters[4], parameters[5]);
      covariance_ = A_.transpose() * L * L.transpose() * A_;
    }

    mat3<double> get_covariance() const {
      return covariance_;
    }
    
    void set_covariance(mat3<double> covariance) {
      typedef scitbx::matrix::cholesky::l_l_transpose_decomposition_in_place<double> decomposition;
      typedef decomposition::accessor_type accessor_type;

      // Check covariance is symmetric
      DIALS_ASSERT(covariance.is_symmetric(1e-7));

      // Get LL*
      mat3<double> Ai = A_.inverse();
      mat3<double> LL = Ai.transpose() * covariance * Ai;

      // Get the lower half of the covariance
      double6 L(LL[0], 
                LL[3], LL[4],
                LL[6], LL[7], LL[8]);

      // Perform an in place cholesky decompositon to get L for COV = LL*
      decomposition cholesky(scitbx::af::ref< double, accessor_type >(&L[0], accessor_type(3)));
      DIALS_ASSERT(!cholesky.failure.failed);

      // Set the parameters
      parameters_ = L;
    }

  private:

    double6 parameters_;
    mat3<double> A_;
    mat3<double> covariance_;
  };
  
  class MosaicBlockAngularSpread {
  public:

    MosaicBlockAngularSpread(vec3<double> r)
      : r_(r),
        parameter_(0.0) {
      DIALS_ASSERT(r_.length() > 0);
      set_parameter(0.0);
    }

    MosaicBlockAngularSpread(vec3<double> r, double parameter)
      : r_(r),
        parameter_(0.0) {
      DIALS_ASSERT(r_.length() > 0);
      set_parameter(parameter);
    }
    
    vec3<double> get_r() const {
      return r_;
    }

    double get_parameter() const {
      return parameter_;
    }

    void set_parameter(double parameter) {
     
      // Ensure the angle is small
      DIALS_ASSERT(detail::small_angle(parameter));

      // Save the parameters
      parameter_ = parameter;
      
      // The eigen values
      double variance = 2.0 * r_.length() * std::tan(parameter_/2.0);

      // Matrix where diagonals are eigen balues
      mat3<double> D = mat3<double>(
          0.0,      0.0,      0.0,
          0.0, variance,      0.0,
          0.0,      0.0, variance);

      // Get the matrix of eigen vectors
      mat3<double> L = eigen_vector_matrix();

      // Construct the covariance matrix
      covariance_ = L * D * L.inverse();
    }

    mat3<double> get_covariance() const {
      return covariance_;
    }

    void set_covariance(mat3<double> covariance) {

      // Get the matrix of eigen vectors
      mat3<double> L = eigen_vector_matrix();

      // Compute the Eigen value matrix and ensure that all the elements
      // except the first diagonal element are zero
      mat3<double> D = L.inverse() * covariance * L;
      DIALS_ASSERT(std::abs(D[0]) < 1e-7);
      DIALS_ASSERT(std::abs(D[1]) < 1e-7);
      DIALS_ASSERT(std::abs(D[2]) < 1e-7);
      DIALS_ASSERT(std::abs(D[3]) < 1e-7);
      DIALS_ASSERT(std::abs(D[5]) < 1e-7);
      DIALS_ASSERT(std::abs(D[6]) < 1e-7);
      DIALS_ASSERT(std::abs(D[7]) < 1e-7);
      DIALS_ASSERT(std::abs(D[8]-D[4]) < 1e-7);

      // Compute the parameter value
      double variance = D[4];
      double parameter = 2.0 * std::atan2(variance / 2.0, r_.length());

      // Set the parameter
      set_parameter(parameter);
    }

    mat3<double> eigen_vector_matrix() const {
      
      // The eigen vectors
      vec3<double> x = r_.normalize();

      // Find two vectors orthogonal
      vec3<double> y = detail::select_orthogonal(x);
      vec3<double> z = x.cross(y);

      // Matrix where each column is an eigenvector
      return mat3<double>(
          x[0], x[1], x[2],
          y[0], y[1], y[2],
          z[0], z[1], z[2]).transpose();
    }

  protected:

    vec3<double> r_;
    double parameter_;
    mat3<double> covariance_;

  };


  class WavelengthSpread {
  public:
    
    WavelengthSpread(vec3<double> s0, vec3<double> r)
      : s0_(s0),
        r_(r),
        parameter_(0.0) {
      DIALS_ASSERT(s0_.length() > 0);
      DIALS_ASSERT(r_.length() > 0);
      set_parameter(0.0);
    }

    WavelengthSpread(vec3<double> s0, vec3<double> r, double parameter)
      : s0_(s0),
        r_(r),
        parameter_(0.0) {
      DIALS_ASSERT(s0_.length() > 0);
      DIALS_ASSERT(r_.length() > 0);
      set_parameter(parameter);
    }
    
    vec3<double> get_s0() const {
      return s0_;
    }

    vec3<double> get_r() const {
      return r_;
    }

    double get_parameter() const {
      return parameter_;
    }

    void set_parameter(double parameter) {
      
      // Save the parameters
      parameter_ = parameter;
      
      // The eigen values
      double ratio = r_.length() / s0_.length();
      double variance = parameter_ * (ratio * ratio);

      // Matrix where diagonals are eigen balues
      mat3<double> D = mat3<double>(
          variance, 0.0, 0.0,
               0.0, 0.0, 0.0,
               0.0, 0.0, 0.0);
      
      // Get the matrix of eigen vectors
      mat3<double> L = eigen_vector_matrix();

      // Construct the covariance matrix
      covariance_ = L * D * L.inverse();
    }

    mat3<double> get_covariance() const {
      return covariance_;
    }

    void set_covariance(mat3<double> covariance) {

      // Get the matrix of eigen vectors
      mat3<double> L = eigen_vector_matrix();

      // Compute the Eigen value matrix and ensure that all the elements
      // except the first diagonal element are zero
      mat3<double> D = L.inverse() * covariance * L;
      for (std::size_t i = 1; i < 9; ++i) {
        DIALS_ASSERT(std::abs(D[i]) < 1e-7);
      }

      // Compute the parameter value
      double ratio = r_.length() / s0_.length();
      double variance = D[0];
      double parameter = variance / (ratio * ratio);

      // Set the parameter
      set_parameter(parameter);
    }
    
    mat3<double> eigen_vector_matrix() const {
      
      // The eigen vectors
      vec3<double> x = r_.normalize();

      // Find two vectors orthogonal
      vec3<double> y = detail::select_orthogonal(x);
      vec3<double> z = x.cross(y);

      // Matrix where each column is an eigenvector
      return mat3<double>(
          x[0], x[1], x[2],
          y[0], y[1], y[2],
          z[0], z[1], z[2]).transpose();
    }

  protected:

    vec3<double> s0_;
    vec3<double> r_;
    double parameter_;
    mat3<double> covariance_;

  };

  
  class CovarianceMatrix {
  public:

    CovarianceMatrix(
          mat3<double> A, 
          vec3<double> s0, 
          vec3<double> r,
          bool use_mosaic_block_angular_spread,
          bool use_wavelength_spread)
      : A_(A),
        s0_(s0),
        r_(r),
        use_mosaic_block_angular_spread_(use_mosaic_block_angular_spread),
        use_wavelength_spread_(use_wavelength_spread),
        reciprocal_lattice_point_spread_(A),
        mosaic_block_angular_spread_(r),
        wavelength_spread_(s0, r) {}

    mat3<double> get_A() const {
      return A_;
    }
    
    vec3<double> get_s0() const {
      return s0_;
    }

    vec3<double> get_r() const {
      return r_;
    }

    bool use_mosaic_block_angular_spread() const {
      return use_mosaic_block_angular_spread_;
    }
    
    bool use_wavelength_spread() const {
      return use_wavelength_spread_;
    }

    std::size_t num_parameters() const {
      return 6 // reciprocal_lattice_point_spread
          +  (use_mosaic_block_angular_spread() ? 1 : 0)
          +  (use_wavelength_spread() ? 1 : 0);
    }

    void set_parameters(const scitbx::af::const_ref<double> parameters) {
      
      DIALS_ASSERT(parameters.size() == num_parameters());
     
      // Set reciprocal lattice point parameters
      reciprocal_lattice_point_spread_.set_parameters(
        double6(
          parameters[0],
          parameters[1],
          parameters[2],
          parameters[3],
          parameters[4],
          parameters[5]));
      std::size_t counter = 6;

      // Setup the covariance matrix
      covariance_ = reciprocal_lattice_point_spread_.get_covariance();

      // Maybe use mosaic block angular spread
      if (use_mosaic_block_angular_spread()) {
        mosaic_block_angular_spread_.set_parameter(parameters[counter]);
        counter += 1;
        covariance_ += mosaic_block_angular_spread_.get_covariance();
      }

      // Maybe use wavelength spread
      if (use_wavelength_spread()) {
        wavelength_spread_.set_parameter(parameters[counter]);
        counter += 1;
        covariance_ += wavelength_spread_.get_covariance();
      }
    }

    scitbx::af::shared<double> get_parameters() const {
      scitbx::af::shared<double> result;

      // Ass reciprocal lattice point spread parameters
      result.insert(result.end(), 
          reciprocal_lattice_point_spread_.get_parameters().begin(),
          reciprocal_lattice_point_spread_.get_parameters().end());

      // Add mosaic block spread parameters
      if (use_mosaic_block_angular_spread()) {
        result.push_back(
            mosaic_block_angular_spread_.get_parameter());
      }

      // Add wavelength spread parameters
      if (use_wavelength_spread()) {
        result.push_back(
            wavelength_spread_.get_parameter());
      }

      // Return parametes
      return result;
    }

    mat3<double> get_covariance() {
      return covariance_;
    }

    const ReciprocalLatticePointSpread& get_reciprocal_lattice_point_spread() const {
      return reciprocal_lattice_point_spread_;
    }
    
    const MosaicBlockAngularSpread& get_mosaic_block_angular_spread() const {
      return mosaic_block_angular_spread_;
    }
    
    const WavelengthSpread& get_wavelength_spread() const {
      return wavelength_spread_;
    }
    
  protected:

    mat3<double> A_;
    vec3<double> s0_;
    vec3<double> r_;
    bool use_mosaic_block_angular_spread_;
    bool use_wavelength_spread_;
    mat3<double> covariance_;
    ReciprocalLatticePointSpread reciprocal_lattice_point_spread_;
    MosaicBlockAngularSpread mosaic_block_angular_spread_;
    WavelengthSpread wavelength_spread_;
  };

  class Model3D {
  public:

    Model3D(
          const Beam &beam,
          const Panel &panel,
          const Goniometer &gonio,
          const Scan &scan,
          mat3<double> sigma,
          vec3<double> r0)
      : s0_(beam.get_s0()),
        m2_(gonio.get_rotation_axis()),
        panel_(panel),
        scan_(scan),
        sigma_(sigma),
        sigma_inv_(sigma.inverse()),
        r0_(r0) {
      double two_pi_cubed = 
        (2*scitbx::constants::pi) *
        (2*scitbx::constants::pi) *
        (2*scitbx::constants::pi);
      double det = std::abs(sigma.determinant());
      DIALS_ASSERT(det > 0);
      normalizing_constant_ = (1.0 / std::sqrt(two_pi_cubed * det));
    }

    double f(double x, double y, double z) const {
      vec3<double> r = coord(x, y, z);
      double d = (r - r0_) * sigma_inv_ * (r - r0_);
      return normalizing_constant_ * std::exp(-0.5 * d);
    }

    double J(double x, double y, double z) const {
      vec3<double> dr_dx = (coord(x+0.5,y,z) - coord(x-0.5,y,z)) / 1.0;
      vec3<double> dr_dy = (coord(x,y+0.5,z) - coord(x,y-0.5,z)) / 1.0;
      vec3<double> dr_dz = (coord(x,y,z+0.5) - coord(x,y,z-0.5)) / 1.0;
      return std::abs(mat3<double>(
        dr_dx[0], dr_dy[0], dr_dz[0],
        dr_dx[1], dr_dy[1], dr_dz[1],
        dr_dx[2], dr_dy[2], dr_dz[2]).determinant());
    }

    vec3<double> coord(double x, double y, double z) const {
      vec3<double> s1 = panel_.get_pixel_lab_coord(vec2<double>(x,y)).normalize();
      s1 = s1 * s0_.length();
      double angle = scan_.get_angle_from_array_index(z);
      vec3<double> r = s1 - s0_;
      mat3<double> R = scitbx::math::r3_rotation::axis_and_angle_as_matrix(m2_, angle);
      return R.transpose() * r;
    }

  protected:

    vec3<double> s0_;
    vec3<double> m2_;
    Panel panel_;
    Scan scan_;
    mat3<double> sigma_;
    mat3<double> sigma_inv_;
    vec3<double> r0_;
    double normalizing_constant_;
  };

  class Model2D {
  public:

    Model2D(
          const Beam &beam,
          const Panel &panel,
          mat3<double> sigma,
          vec3<double> r0)
      : s0_(beam.get_s0()),
        panel_(panel),
        sigma_(sigma),
        sigma_inv_(sigma.inverse()),
        r0_(r0) {
      double two_pi_cubed = 
        (2*scitbx::constants::pi) *
        (2*scitbx::constants::pi) *
        (2*scitbx::constants::pi);
      double det = std::abs(sigma.determinant());
      DIALS_ASSERT(det > 0);
      normalizing_constant_ = (1.0 / std::sqrt(two_pi_cubed * det));
    }

    double f(double x, double y) const {
      vec3<double> r = coord(x, y);
      double d = (r - r0_) * sigma_inv_ * (r - r0_);
      return normalizing_constant_ * std::exp(-0.5 * d);
    }

    double J(double x, double y) const {
      vec3<double> dr_dx = (coord(x+0.5,y) - coord(x-0.5,y)) / 1.0;
      vec3<double> dr_dy = (coord(x,y+0.5) - coord(x,y-0.5)) / 1.0;
      return std::abs(mat3<double>(
        dr_dx[0], dr_dy[0], 0.0,
        dr_dx[1], dr_dy[1], 0.0,
        dr_dx[2], dr_dy[2], 0.0).determinant());
    }

    vec3<double> coord(double x, double y) const {
      vec3<double> s1 = panel_.get_pixel_lab_coord(vec2<double>(x,y)).normalize();
      s1 = s1 * s0_.length();
      vec3<double> r = s1 - s0_;
      return r;
    }

  protected:

    vec3<double> s0_;
    Panel panel_;
    mat3<double> sigma_;
    mat3<double> sigma_inv_;
    vec3<double> r0_;
    double normalizing_constant_;
  };

  class MLTarget2DSingle {
  public:

    MLTarget2DSingle(
          const Beam &beam,
          const Panel &panel,
          mat3<double> sigma,
          vec3<double> r0,
          int6 bbox,
        std::size_t num_integral)
      : model_(
          beam,
          panel,
          sigma,
          r0),
        bbox_(bbox),
        num_integral_(num_integral),
        uniform_(0.0, 1.0) {
      DIALS_ASSERT(num_integral > 0);
      DIALS_ASSERT(bbox[1] > bbox[0]);
      DIALS_ASSERT(bbox[3] > bbox[2]);
      DIALS_ASSERT(bbox[5] == bbox[4] + 1);
    }
    
    double log_likelihood(
        const scitbx::af::const_ref< float, scitbx::af::c_grid<3> > &data,
        const scitbx::af::const_ref< int, scitbx::af::c_grid<3> > &mask) {
      DIALS_ASSERT(mask.accessor().all_eq(data.accessor()));
      std::size_t xsize = (bbox_[1] - bbox_[0]);
      std::size_t ysize = (bbox_[3] - bbox_[2]);
      std::size_t zsize = (bbox_[5] - bbox_[4]);
      DIALS_ASSERT(data.accessor()[0] == zsize);
      DIALS_ASSERT(data.accessor()[1] == ysize);
      DIALS_ASSERT(data.accessor()[2] == xsize);
      int mask_code = Valid | Foreground;
      double logL = 0.0;
      double ntot = 0.0;
      double Ptot = 0.0;
      for (std::size_t j = 0; j < data.accessor()[1]; ++j) {
        for (std::size_t i = 0; i < data.accessor()[2]; ++i) {
          if ((mask(0,j,i) & mask_code) == mask_code) {
            double nj = data(0,j,i);
            double Pj = compute_P(j, i) + 1e-10;
            DIALS_ASSERT(nj >= 0);
            logL += nj * std::log(Pj);
            ntot += nj;
            Ptot += Pj;
          }
        }
      }
      DIALS_ASSERT(ntot > 0);
      DIALS_ASSERT(Ptot > 0);
      logL -= ntot * std::log(Ptot);
      return logL;
    }

    double compute_P(std::size_t j, std::size_t i) {
      int x0 = bbox_[0];
      int y0 = bbox_[2];
      int jj = y0 + j;
      int ii = x0 + i;
      double J = model_.J(ii+0.5, jj+0.5);
      double volume = 1.0;
      double sumf = 0;
      double scale = 1.0 / num_integral_;
      double num_integral_total = num_integral_ * num_integral_;
      for (std::size_t jint = 0; jint < num_integral_; ++jint) {
        for (std::size_t iint = 0; iint < num_integral_; ++iint) {
          double jj0 = jj + jint * scale;
          double ii0 = ii + iint * scale;
          double y = jj0 + uniform_(generator_) * scale;
          double x = ii0 + uniform_(generator_) * scale;
          double f = model_.f(x, y);
          sumf += f;
        }
      }
      return J * volume * sumf / num_integral_total;
    }


  protected:

    Model2D model_;
    int6 bbox_;
    std::size_t num_integral_;
    std::default_random_engine generator_;
    std::uniform_real_distribution<double> uniform_;
  };

  class MLTarget3DSingle {
  public:

    MLTarget3DSingle(
          const Beam &beam,
          const Panel &panel,
          const Goniometer &gonio,
          const Scan &scan,
          mat3<double> sigma,
          vec3<double> r0,
          int6 bbox,
          std::size_t num_integral)
      : model_(
          beam,
          panel,
          gonio,
          scan,
          sigma,
          r0),
        bbox_(bbox),
        num_integral_(num_integral),
        generator_(std::random_device()()),
        uniform_(0.0, 1.0) {
      DIALS_ASSERT(num_integral > 0);
      DIALS_ASSERT(bbox[1] > bbox[0]);
      DIALS_ASSERT(bbox[3] > bbox[2]);
      DIALS_ASSERT(bbox[5] > bbox[4]);
    }
    
    double log_likelihood(
        const scitbx::af::const_ref< float, scitbx::af::c_grid<3> > &data,
        const scitbx::af::const_ref< int, scitbx::af::c_grid<3> > &mask) {
      DIALS_ASSERT(mask.accessor().all_eq(data.accessor()));
      std::size_t xsize = (bbox_[1] - bbox_[0]);
      std::size_t ysize = (bbox_[3] - bbox_[2]);
      std::size_t zsize = (bbox_[5] - bbox_[4]);
      DIALS_ASSERT(data.accessor()[0] == zsize);
      DIALS_ASSERT(data.accessor()[1] == ysize);
      DIALS_ASSERT(data.accessor()[2] == xsize);
      int mask_code = Valid | Foreground;
      double logL = 0.0;
      double ntot = 0.0;
      double Ptot = 0.0;
      for (std::size_t k = 0; k < data.accessor()[0]; ++k) {
        for (std::size_t j = 0; j < data.accessor()[1]; ++j) {
          for (std::size_t i = 0; i < data.accessor()[2]; ++i) {
            if ((mask(k,j,i) & mask_code) == mask_code) {
              double Pj = compute_P(k, j, i) + 1e-10;
              double nj = data(k,j,i);
              DIALS_ASSERT(nj >= 0);
              logL += nj * std::log(Pj);
              ntot += nj;
              Ptot += Pj;
            }
          }
        }
      }
      DIALS_ASSERT(ntot > 0);
      DIALS_ASSERT(Ptot > 0);
      logL -= ntot * std::log(Ptot);
      return logL;
    }

    double compute_P(std::size_t k, std::size_t j, std::size_t i) {
      int x0 = bbox_[0];
      int y0 = bbox_[2];
      int z0 = bbox_[4];
      int kk = z0 + k;
      int jj = y0 + j;
      int ii = x0 + i;
      double J = model_.J(ii+0.5, jj+0.5, kk+0.5);
      double volume = 1.0;
      double sumf = 0;
      double scale = 1.0 / num_integral_;
      double num_integral_total = num_integral_ * num_integral_ * num_integral_;
      for (std::size_t kint = 0; kint < num_integral_; ++kint) {
        for (std::size_t jint = 0; jint < num_integral_; ++jint) {
          for (std::size_t iint = 0; iint < num_integral_; ++iint) {
            double kk0 = kk + kint * scale;
            double jj0 = jj + jint * scale;
            double ii0 = ii + iint * scale;
            double z = kk0 + uniform_(generator_) * scale;
            double y = jj0 + uniform_(generator_) * scale;
            double x = ii0 + uniform_(generator_) * scale;
            double f = model_.f(x, y, z);
            sumf += f;
          }
        }
      }
      return J * volume * sumf / num_integral_total;
    }


  protected:

    Model3D model_;
    int6 bbox_;
    std::size_t num_integral_;
    std::default_random_engine generator_;
    std::uniform_real_distribution<double> uniform_;
  };

  class MLTarget3D {
  public:

    MLTarget3D(
        const Experiment &experiment,
        const reflection_table &reflections,
        std::size_t num_integral,
        bool use_mosaic_block_angular_spread,
        bool use_wavelength_spread)
      : experiment_(experiment),
        reflections_(reflections),
        num_integral_(num_integral),
        use_mosaic_block_angular_spread_(use_mosaic_block_angular_spread),
        use_wavelength_spread_(use_wavelength_spread) {
    }

    double log_likelihood(const scitbx::af::const_ref<double> &parameters) const {

      // Check the input
      DIALS_ASSERT(experiment_.get_beam());
      DIALS_ASSERT(experiment_.get_detector());
      DIALS_ASSERT(experiment_.get_goniometer());
      DIALS_ASSERT(experiment_.get_scan());
      DIALS_ASSERT(reflections_.contains("shoebox"));
      DIALS_ASSERT(reflections_.contains("s1"));

      // Get arrays from reflections
      af::const_ref< Shoebox<> > sbox = 
        reflections_.get< Shoebox<> >("shoebox").const_ref();
      af::const_ref< cctbx::miller::index<> > h = 
        reflections_.get< cctbx::miller::index<> >("miller_index").const_ref();

      // Get the beam vector
      mat3<double> A = experiment_.get_crystal()->get_A();
      vec3<double> s0 = experiment_.get_beam()->get_s0();

      // Compute the sum of the logL from each reflection
      double logL = 0.0;
      for (std::size_t i = 0; i < reflections_.size(); ++i) {

        // Grad some data for the reflection
        scitbx::af::const_ref<float, scitbx::af::c_grid<3> > data = sbox[i].data.const_ref();
        scitbx::af::const_ref<int,   scitbx::af::c_grid<3> > mask = sbox[i].mask.const_ref();
        int6 bbox = sbox[i].bbox;
        std::size_t panel = sbox[i].panel;
        vec3<double> r0 = A * h[i];

        // Initialize the covariance matrix class
        CovarianceMatrix cov(
            A,
            s0,
            r0, 
            use_mosaic_block_angular_spread_,
            use_wavelength_spread_);

        // Set the parameters to generate the covariance matrix
        cov.set_parameters(parameters);

        // Create the target function class
        MLTarget3DSingle target(
            (*experiment_.get_beam()),
            (*experiment_.get_detector())[panel],
            (*experiment_.get_goniometer()),
            (*experiment_.get_scan()),
            cov.get_covariance(),
            r0,
            bbox,
            num_integral_);        

        // Compute the loglikelihood for the reflection
        logL += target.log_likelihood(data, mask);
      }
      return logL;
    }

  protected:

    Experiment experiment_;
    reflection_table reflections_;
    std::size_t num_integral_;
    bool use_mosaic_block_angular_spread_;
    bool use_wavelength_spread_;
  };

}} // namespace dials::algorithms

#endif // DIALS_ALGORITHMS_PROFILE_MODELLING_TARGET_H

