/*
 * ext.cc
 *
 *  Copyright (C) 2018 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dials/algorithms/profile_model/gaussian_rs/coordinate_system.h>
#include <dials/algorithms/shoebox/mask_code.h>
#include <dials/array_family/reflection_table.h>
#include <dials/error.h>
#include <dials/model/data/shoebox.h>
#include <dxtbx/model/experiment.h>
#include <scitbx/constants.h>
#include <scitbx/mat2.h>
#include <scitbx/mat3.h>
#include <scitbx/math/r3_rotation.h>
#include <scitbx/matrix/eigensystem.h>
#include <scitbx/matrix/multiply.h>
#include <scitbx/sym_mat2.h>
#include <scitbx/vec2.h>
#include <scitbx/vec3.h>

using namespace boost::python;

namespace dials {
namespace algorithms {
namespace boost_python {

using boost::math::erf;
using scitbx::vec2;
using scitbx::constants::pi;

vec2<double> truncated_normal_first_and_second_moments(double mu,
                                                       double sigma) {
  /* double E_p = 0; */
  /* double E_p2 = 0; */
  /* if (mu > 1) { */
  /*   E_p = 1.0; */
  /*   E_p2 = E_p*E_p; */
  /* } else if (mu < 0) { */
  /*   E_p = 0.0; */
  /*   E_p2 = E_p*E_p; */
  /* } else { */
  /*   E_p = mu; */
  /*   E_p2 = E_p*E_p; */
  /* } */
  double alpha = (0 - mu) / sigma;
  double beta = (1 - mu) / sigma;
  double e1 = (1 / std::sqrt(2 * pi)) * std::exp(-0.5 * alpha * alpha);
  double e2 = (1 / std::sqrt(2 * pi)) * std::exp(-0.5 * beta * beta);
  double e3 = 0.5 * (1 + erf(alpha / std::sqrt(2)));
  double e4 = 0.5 * (1 + erf(beta / std::sqrt(2)));
  double Z = e4 - e3;
  double E_p = 0;
  double E_p2 = 0;
  if (Z > 1e-5) {
    E_p = mu + (e1 - e2) * sigma / Z;
    E_p2 =
        mu * mu + sigma * sigma - ((mu + 1) * e2 - (mu + 0) * e1) * sigma / Z;
  } else if (mu > 1) {
    E_p = 1.0;
    E_p2 = 1.0;
  } else {
    E_p = 0.0;
    E_p2 = 0.0;
  }
  DIALS_ASSERT(E_p >= 0 && E_p <= 1);
  DIALS_ASSERT(E_p2 >= E_p * E_p);
  return vec2<double>(E_p, E_p2);
}

vec2<double> compute_E(double rho0, double var_rho0, double mu, double I,
                       double V) {
  double var_rho1 = 1.0 / (1 / var_rho0 + mu * mu / V);
  double rho1 = (rho0 / var_rho0 + I * mu / V) * var_rho1;
  return truncated_normal_first_and_second_moments(rho1, std::sqrt(var_rho1));
}

void estimate_partiality(af::const_ref<double> I, af::const_ref<double> V,
                         af::const_ref<double> S, af::const_ref<double> M,
                         double rho0, double var_rho0, af::ref<double> E_p,
                         af::ref<double> E_p2) {
  DIALS_ASSERT(I.size() == V.size());
  DIALS_ASSERT(I.size() == S.size());
  DIALS_ASSERT(I.size() == M.size());
  DIALS_ASSERT(I.size() == E_p.size());
  DIALS_ASSERT(I.size() == E_p2.size());
  for (std::size_t i = 0; i < I.size(); ++i) {
    vec2<double> E = compute_E(rho0, var_rho0, S[i] * M[i], I[i], V[i]);
    E_p[i] = E[0];
    E_p2[i] = E[1];
  }
}

af::shared<double>
estimate_intensity(af::const_ref<double> I, af::const_ref<double> V,
                   af::const_ref<double> S, af::const_ref<double> E_p,
                   af::const_ref<double> E_p2, af::ref<double> M,
                   boost::python::list groups) {
  af::shared<double> Iest;
  for (std::size_t j = 0; j < boost::python::len(groups); ++j) {
    af::const_ref<std::size_t> indices =
        boost::python::extract<af::const_ref<std::size_t>>(groups[j]);
    double sum1 = 0;
    double sum2 = 0;
    for (std::size_t i = 0; i < indices.size(); ++i) {
      std::size_t idx = indices[i];
      DIALS_ASSERT(V[idx] > 0);
      sum1 += I[idx] * S[idx] * E_p[idx] / V[idx];
      sum2 += S[idx] * S[idx] * E_p2[idx] / V[idx];
    }
    DIALS_ASSERT(sum2 > 0);
    double II = sum1 / sum2;
    for (std::size_t i = 0; i < indices.size(); ++i) {
      std::size_t idx = indices[i];
      M[idx] = II;
    }
    Iest.push_back(II);
  }
  return Iest;
}

af::shared<double>
estimate_scale_factor(af::const_ref<double> I, af::const_ref<double> V,
                      af::const_ref<double> M, af::const_ref<double> E_p,
                      af::const_ref<double> E_p2, af::ref<double> S,
                      boost::python::list groups) {
  af::shared<double> Sest;
  for (std::size_t j = 0; j < boost::python::len(groups); ++j) {
    af::const_ref<std::size_t> indices =
        boost::python::extract<af::const_ref<std::size_t>>(groups[j]);
    double sum1 = 0;
    double sum2 = 0;
    for (std::size_t i = 0; i < indices.size(); ++i) {
      std::size_t idx = indices[i];
      DIALS_ASSERT(V[idx] > 0);
      sum1 += I[idx] * M[idx] * E_p[idx] / V[idx];
      sum2 += M[idx] * M[idx] * E_p2[idx] / V[idx];
    }
    DIALS_ASSERT(sum2 > 0);
    double SS = sum1 / sum2;
    for (std::size_t i = 0; i < indices.size(); ++i) {
      std::size_t idx = indices[i];
      S[idx] = SS;
    }
    Sest.push_back(SS);
  }
  return Sest;
}

BOOST_PYTHON_MODULE(dials_scratch_jmp_merge_ext) {
  def("estimate_partiality", &estimate_partiality);
  def("estimate_intensity", &estimate_intensity);
  def("estimate_scale_factor", &estimate_scale_factor);
}

} // namespace boost_python
} // namespace algorithms
} // namespace dials
