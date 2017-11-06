/*
 * model.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */


#ifndef DLSTBX_ALGORITHMS_PROFILE_MODEL_NAVE_MODEL_H
#define DLSTBX_ALGORITHMS_PROFILE_MODEL_NAVE_MODEL_H

#include <scitbx/vec3.h>
#include <scitbx/mat3.h>
#include <scitbx/math/r3_rotation.h>
#include <scitbx/array_family/small.h>
#include <cctbx/miller.h>
#include <dials/algorithms/spot_prediction/rotation_angles.h>
#include <dxtbx/model/ray_intersection.h>

namespace af = scitbx::af;

namespace dlstbx { namespace algorithms {

  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::mat3;
  using scitbx::math::r3_rotation::axis_and_angle_as_matrix;
  using dxtbx::model::plane_ray_intersection;
  using dials::algorithms::rotation_angles;

  template <typename T>
  T sqr(T a) {
    return a * a;
  }

  /**
   * A class to represent the nave profile model
   */
  class Model {
  public:

    /**
     * Initialise the model for a reflection
     * @param D the detector d matrix
     * @param A the UB matrix
     * @param s0 The incident beam vector
     * @param m2 The rotation axis
     * @param s1 The diffracted beam vector
     * @param phi The rotation angle
     * @param sig_s The mosiac block size parameters
     * @param sig_a The unit_cell size parameters
     * @param sig_w The angular spread parameters
     */
    Model(mat3<double> D,
          mat3<double> A,
          vec3<double> s0,
          vec3<double> m2,
          vec3<double> s1,
          double phi,
          vec3<double> sig_s,
          vec3<double> sig_a,
          vec3<double> sig_w)
      : D_(D),
        A_(A),
        s0_(s0),
        m2_(m2.normalize()),
        s1_(s1.normalize() * s0.length()),
        phi0_(phi),
        rlp_(R(-phi) * (s1 - s0)) {

      // Check the input
      DIALS_ASSERT(s0_.length() > 0);
      DIALS_ASSERT(sig_s.const_ref().all_ge(0));
      DIALS_ASSERT(sig_a.const_ref().all_ge(0));
      DIALS_ASSERT(sig_w.const_ref().all_ge(0));

      // The covariance matrix for the mosaic block size component in the
      // reciprocal lattice coordinate system
      mat3<double> sigma_s(
        sqr(sig_s[0]), 0, 0,
        0, sqr(sig_s[1]), 0,
        0, 0, sqr(sig_s[2])
      );

      // The covariance matrix for the spread in unit cell size component in the
      // reciprocal lattice coordinate system
      vec3<double> h = A.inverse() * rlp_;
      mat3<double> sigma_a(
        sqr(h[0]*sig_a[0]), 0, 0,
        0, sqr(h[1]*sig_a[1]), 0,
        0, 0, sqr(h[1]*sig_a[2])
      );

      // Compute the covariance in orthognal coordinate system
      sigma_s = A * sigma_s * A.transpose();
      sigma_a = A * sigma_a * A.transpose();

      // Select two vectors orthogonal to the rlp
      vec3<double> rn = rlp_.normalize();
      vec3<double> v1 = std::abs(rn[0]) > std::abs(rn[2])
        ? vec3<double>(-rn[1], rn[0], 0.0).normalize()
        : vec3<double>(0.0, -rn[2], rn[1]).normalize();
      vec3<double> v2 = rn.cross(v1).normalize();
      vec3<double> v3 = rn.cross(v2).normalize();

      // Construct an eigenvector matrix
      mat3<double> U(
        v2[0], v3[0], rn[0],
        v2[1], v3[1], rn[1],
        v2[2], v3[2], rn[2]
      );

      // Compute the angular spread
      double w = ((A * mat3<double>(
        sig_w[0], 0, 0,
        0, sig_w[1], 0,
        0, 0, sig_w[2])) * rlp_).length();

      // Construct an eigenvalue matrix
      mat3<double> V(
        w*w, 0, 0,
        0, w*w, 0,
        0,   0, 0
      );

      // Compute the covariance matrix for the angular spread of mosaic blocks
      // in the orthogonal lab coordinate system using the eigenvectors and
      // eigenvalues to produce a 2d gaussian in the plane normal to the rlp.
      mat3<double> sigma_w = U*V*U.transpose();

      // The full covariance matrix and its inverse
      sigma_ = sigma_s + sigma_a + sigma_w;
      sigma_inv_ = sigma_.inverse();
    }

    /** @returns The detector d matrix */
    mat3<double> D() const {
      return D_;
    }

    /** @returns The A matrix */
    mat3<double> A() const {
      return A_;
    }

    /** @returns The incident beam vector */
    vec3<double> s0() const {
      return s0_;
    }

    /** @returns The rotation axis */
    vec3<double> m2() const {
      return m2_;
    }

    /** @returns The diffracted beam vector */
    vec3<double> s1() const {
      return s1_;
    }

    /** @returns The rotation axis */
    double phi0() const {
      return phi0_;
    }

    /** @returns The reciprocal lattice point */
    vec3<double> rlp() const {
      return rlp_;
    }

    /** @returns The covariance matrix */
    mat3<double> sigma() const {
      return sigma_;
    }

    /** @returns The inverse covariance matrix */
    mat3<double> sigma_inv() const {
      return sigma_inv_;
    }

    /** @returns The rotation matrix about m2 */
    mat3<double> R(double phi) const {
      return axis_and_angle_as_matrix(m2_, phi);
    }

    /** @returns The rlp at a point */
    vec3<double> r(double x, double y, double phi) const {
      vec3<double> v = D_ * vec3<double>(x, y, 1.0);
      double slen = s0_.length();
      double vlen = v.length();
      DIALS_ASSERT(vlen > 0);
      return R(phi).transpose() * (v * slen / vlen - s0_);
    }

    /** @returns The mahanabonis distance */
    double Dm(double x, double y, double phi) const {
      vec3<double> dh = r(x, y, phi) - rlp_;
      return dh * sigma_inv_ * dh;
    }

    /** @returns The probability density */
    double P(double x, double y, double phi) const {
      return std::exp(-0.5 * Dm(x, y, phi));
    }

  private:

    mat3<double> D_;
    mat3<double> A_;
    vec3<double> s0_;
    vec3<double> m2_;
    vec3<double> s1_;
    double phi0_;
    vec3<double> rlp_;
    mat3<double> sigma_;
    mat3<double> sigma_inv_;
  };

}}

#endif // DLSTBX_ALGORITHMS_PROFILE_MODEL_NAVE_MODEL_H
