/*
 * target_forward.h
 *
 *  Copyright (C) 2017 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_PROFILE_MODELLING_TARGET_FORWARD_H
#define DIALS_ALGORITHMS_PROFILE_MODELLING_TARGET_FORWARD_H

#include <dials/algorithms/spot_prediction/ray_predictor.h>
#include <dials/array_family/reflection_table.h>
#include <dials/error.h>
#include <dials/model/data/mask_code.h>
#include <dxtbx/model/beam.h>
#include <dxtbx/model/detector.h>
#include <dxtbx/model/experiment.h>
#include <dxtbx/model/goniometer.h>
#include <dxtbx/model/scan.h>
#include <random>
#include <scitbx/array_family/accessors/c_grid.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/simple_tiny_io.h>
#include <scitbx/array_family/tiny_types.h>
#include <scitbx/constants.h>
#include <scitbx/mat3.h>
#include <scitbx/math/r3_rotation.h>
#include <scitbx/matrix/cholesky.h>
#include <scitbx/vec2.h>
#include <scitbx/vec3.h>

namespace dials {
namespace algorithms {

using dials::af::reflection_table;
using dials::model::Foreground;
using dials::model::Shoebox;
using dials::model::Valid;
using dxtbx::model::Beam;
using dxtbx::model::Experiment;
using dxtbx::model::Goniometer;
using dxtbx::model::Panel;
using dxtbx::model::Scan;
using scitbx::mat3;
using scitbx::vec2;
using scitbx::vec3;
using scitbx::af::double6;
using scitbx::af::int6;

/**
 * Class for reciprocal lattice point distriubtion
 */
class ReciprocalLatticePointDistribution {
public:
  ReciprocalLatticePointDistribution(miller_index h0, double6 parameters)
      : h0_(h0), normal_(0.0, 1.0) {
    set_parameters(parameters);
  }

  /**
   * @returns The central miller index
   */
  miller_index h0() const { return h0_; }

  /**
   * @returns The parameters
   */
  double6 get_parameters() const { return parameters_; }

  /**
   * Set the parameters
   */
  void set_parameters(double6 parameters) {
    // Set parameters
    parameters_ = parameters;

    // The parameters are given as 6 parameters in the lower triangular matrix
    // form and the covariance matrix is given by the Cholesky decomposition
    // COV = LL*
    mat3<double> L(parameters[0], 0.0, 0.0, parameters[1], parameters[2], 0.0,
                   parameters[3], parameters[4], parameters[5]);
    covariance_ = L * L.transpose();
  }

  /**
   * Get the covariance matrix
   */
  mat3<double> get_covariance() const { return covariance_; }

  /**
   * Set the covariance matrix
   */
  void set_covariance(mat3<double> covariance) {
    typedef scitbx::matrix::cholesky::l_l_transpose_decomposition_in_place<
        double>
        decomposition;
    typedef decomposition::accessor_type accessor_type;

    // Check covariance is symmetric
    DIALS_ASSERT(covariance.is_symmetric(1e-7));

    // Get LL*
    mat3<double> LL = covariance;

    // Get the lower half of the covariance
    double6 L(LL[0], LL[3], LL[4], LL[6], LL[7], LL[8]);

    // Perform an in place cholesky decompositon to get L for COV = LL*
    decomposition cholesky(
        scitbx::af::ref<double, accessor_type>(&L[0], accessor_type(3)));
    DIALS_ASSERT(!cholesky.failure.failed);

    // Set the parameters
    parameters_ = L;
  }

  vec3<double> sample() {
    // The parameters are given as 6 parameters in the lower triangular matrix
    // form and the covariance matrix is given by the Cholesky decomposition
    // COV = LL*
    mat3<double> L(parameters_[0], 0.0, 0.0, parameters_[1], parameters_[2],
                   0.0, parameters_[3], parameters_[4], parameters_[5]);

    // Generate a normal random variate
    return L * vec3<double>(normal_(generator_), normal_(generator_),
                            normal_(generator_)) +
           vec3<double>(h0_[0], h0_[1], h0_[2]);
  }

private:
  miller_index h0_;
  double6 parameters_;
  mat3<double> covariance_;
  std::default_random_engine generator_;
  std::normal_distribution<double> normal_;
};

/**
 * A class to represent the wavelength distribution
 */
class WavelengthDistribution {
public:
  WavelengthDistribution(double wavelength, double parameters)
      : wavelength_(wavelength) {
    set_parameters(parameters);
  }

  /**
   * @returns The wavelength
   */
  double wavelength() const { return wavelength_; }

  /**
   * @returns The parameters
   */
  double get_parameters() const { return parameters_; }

  /**
   * Set the parameters
   */
  void set_parameters(double parameters) {
    DIALS_ASSERT(parameters > 0);
    parameters_ = parameters;
  }

  double sample() {
    std::normal_distribution<double> normal(wavelength_, parameters_);
    return normal(generator_);
  }

protected:
  double wavelength_;
  double parameters_;
  std::default_random_engine generator_;
};

class BeamVectorDistribution {
public:
  BeamVectorDistribution(vec3<double> s0, double parameters)
      : s0_(s0), normalized_s0_(s0.normalize()),
        wavelength_distribution_(1.0 / s0.length(), parameters) {}

  /**
   * @returns The wavelength
   */
  vec3<double> s0() const { return s0_; }

  /**
   * @returns The parameters
   */
  double get_parameters() const {
    return wavelength_distribution_.get_parameters();
  }

  /**
   * Set the parameters
   */
  void set_parameters(double parameters) {
    wavelength_distribution_.set_parameters(parameters);
  }

  vec3<double> sample() {
    double wavelength = wavelength_distribution_.sample();
    return normalized_s0_ / wavelength;
  }

private:
  vec3<double> s0_;
  vec3<double> normalized_s0_;
  WavelengthDistribution wavelength_distribution_;
};

/**
 * Maximum likelihood target function
 */
class MLTargetForward3DSingle {
public:
  /**
   * Init
   */
  MLTargetForward3DSingle(
      const Panel &panel, const Goniometer &gonio, const Scan &scan,
      bool entering, int6 bbox,
      ReciprocalLatticePointDistribution reciprocal_lattice_point_distribution,
      BeamVectorDistribution beam_vector_distribution, std::size_t num_sample)
      : m2_(gonio.get_rotation_axis()),
        fixed_rotation_(gonio.get_fixed_rotation()),
        setting_rotation_(gonio.get_setting_rotation()), panel_(panel),
        scan_(scan), entering_(entering), bbox_(bbox),
        reciprocal_lattice_point_distribution_(
            reciprocal_lattice_point_distribution),
        beam_vector_distribution_(beam_vector_distribution),
        model_(simulate()) {
    DIALS_ASSERT(num_sample > 0);
    DIALS_ASSERT(bbox[1] > bbox[0]);
    DIALS_ASSERT(bbox[3] > bbox[2]);
    DIALS_ASSERT(bbox[5] > bbox[4]);
  }

  /**
   * Compute the log likelihood
   */
  double log_likelihood(
      const scitbx::af::const_ref<float, scitbx::af::c_grid<3>> &data,
      const scitbx::af::const_ref<int, scitbx::af::c_grid<3>> &mask) {
    DIALS_ASSERT(mask.accessor().all_eq(data.accessor()));
    std::size_t xsize = (bbox_[1] - bbox_[0]);
    std::size_t ysize = (bbox_[3] - bbox_[2]);
    std::size_t zsize = (bbox_[5] - bbox_[4]);
    DIALS_ASSERT(data.accessor()[0] == zsize);
    DIALS_ASSERT(data.accessor()[1] == ysize);
    DIALS_ASSERT(data.accessor()[2] == xsize);
    scitbx::af::const_ref<double, scitbx::af::c_grid<3>> model =
        model_.const_ref();
    int mask_code = Valid | Foreground;
    double logL = 0.0;
    double ntot = 0.0;
    double Ptot = 0.0;
    for (std::size_t k = 0; k < data.accessor()[0]; ++k) {
      for (std::size_t j = 0; j < data.accessor()[1]; ++j) {
        for (std::size_t i = 0; i < data.accessor()[2]; ++i) {
          if ((mask(k, j, i) & mask_code) == mask_code) {
            double Pj = model(k, j, i) + 1e-10;
            double nj = data(k, j, i);
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

  /**
   * Simulate
   */
  scitbx::af::versa<double, scitbx::af::c_grid<3>> simulate() {
    std::size_t xsize = (bbox_[1] - bbox_[0]);
    std::size_t ysize = (bbox_[3] - bbox_[2]);
    std::size_t zsize = (bbox_[5] - bbox_[4]);
    int x0 = bbox_[0];
    int y0 = bbox_[2];
    int z0 = bbox_[4];
    if (model_.size() == 0) {

      // Allocate array
      scitbx::af::c_grid<3> grid(zsize, ysize, xsize);
      model_ = scitbx::af::versa<double, scitbx::af::c_grid<3>>(grid);
      scitbx::af::ref<double, scitbx::af::c_grid<3>> model = model_.ref();

      // Create the ray predictor
      dials::algorithms::RayPredictor ray_predictor(m2_, fixed_rotation_,
                                                    setting_rotation_);

      // Loop through number of rays
      for (std::size_t i = 0; i < num_sample_; ++i) {

        // Sample from distributions
        vec3<double> h = reciprocal_lattice_point_distribution_.sample();
        vec3<double> s0 = beam_vector_distribution_.sample();

        try {
          Ray ray = ray_predictor(A_, s0, h, entering_);
          double z = scan_.get_array_index_from_angle(ray.angle);
          vec2<double> pxy = panel_.get_ray_intersection_px(ray.s1);
          int i = (int)std::floor((pxy[0] - x0));
          int j = (int)std::floor((pxy[1] - y0));
          int k = (int)std::floor((z - z0));
          if (k >= 0 && k < zsize && j >= 0 && j < ysize && i >= 0 &&
              i < xsize) {
            model(k, j, i) += 1;
          }
        } catch (dials::error) {
          continue;
        }
      }

      // Normalize the model
      for (std::size_t i = 0; i < model.size(); ++i) {
        model[i] /= num_sample_;
      }

    } else {
      DIALS_ASSERT(model_.accessor()[0] == zsize);
      DIALS_ASSERT(model_.accessor()[1] == ysize);
      DIALS_ASSERT(model_.accessor()[2] == xsize);
    }
    return model_;
  }

protected:
  vec3<double> m2_;
  mat3<double> fixed_rotation_;
  mat3<double> setting_rotation_;
  mat3<double> A_;
  Panel panel_;
  Scan scan_;
  miller_index h0_;
  bool entering_;
  int6 bbox_;
  std::size_t num_sample_;
  ReciprocalLatticePointDistribution reciprocal_lattice_point_distribution_;
  BeamVectorDistribution beam_vector_distribution_;
  scitbx::af::versa<double, scitbx::af::c_grid<3>> model_;
};

} // namespace algorithms
} // namespace dials

#endif // DIALS_ALGORITHMS_PROFILE_MODELLING_TARGET_FORWARD_H
