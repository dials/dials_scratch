/*
 * use_of_sparse.h
 *
 *  Copyright (C) (2016) STFC Rutherford Appleton Laboratory, UK.
 *
 *  Author: David Waterman
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */

#ifndef DIALS_SCRATCH_EXAMPLES_USE_OF_SPARSE_H
#define DIALS_SCRATCH_EXAMPLES_USE_OF_SPARSE_H
#include <scitbx/sparse/matrix.h>
#include <dials/error.h>

namespace dials_scratch { namespace examples {

  using scitbx::sparse::matrix;

  // check return of a sparse matrix to Python
  matrix<double> create_sparse_matrix() {
    matrix<double> result(5, 3);
    // set a value
    result(0, 0) =  1.;
    return result;
  }

}} // namespace dials_scratch::examples

#endif // DIALS_SCRATCH_EXAMPLES_USE_OF_SPARSE_H
