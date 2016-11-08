/*
 * array_gotchas.h
 *
 *  Copyright (C) (2016) STFC Rutherford Appleton Laboratory, UK.
 *
 *  Author: David Waterman
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */

#ifndef DIALS_SCRATCH_EXAMPLES_ARRAY_GOTCHAS_H
#define DIALS_SCRATCH_EXAMPLES_ARRAY_GOTCHAS_H

#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/ref.h>
#include <dials/error.h>
#include <cctbx/miller.h>

// For debugging
// #include <scitbx/array_family/simple_io.h>
// then can do things like this:
// std::cout << something.const_ref() << "\n";

namespace dials_scratch { namespace examples {

  typedef cctbx::miller::index<> miller_index;

  class BadBucket {
  public:
    BadBucket(scitbx::af::const_ref< miller_index > hkl)
      : hkl_(hkl) {}
    // const_ref getter. No conversion to Python type!
    scitbx::af::const_ref<miller_index> get_const_ref_hkl() const {
      return hkl_;
    };
    // shared getter, makes a copy of the data, but still subject to
    // externally-determined lifetime of hkl_.
    scitbx::af::shared<miller_index> get_shared_hkl() const {
      return scitbx::af::shared<miller_index>(hkl_.begin(), hkl_.end());
    };

  private:
    scitbx::af::const_ref< miller_index > hkl_;
  };

}} // namespace dials_scratch::examples

#endif // DIALS_SCRATCH_EXAMPLES_ARRAY_GOTCHAS_H
