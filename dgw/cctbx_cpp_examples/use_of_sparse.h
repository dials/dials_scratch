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
#include <scitbx/sparse/vector.h>
#include <dials/error.h>

namespace dials_scratch { namespace examples {

  using scitbx::sparse::matrix;
  using scitbx::sparse::vector;

  // check return of a sparse matrix to Python
  matrix<double> create_sparse_matrix() {
    matrix<double> result(5, 3);
    // set a value
    result(0, 0) =  1.;
    return result;
  }

  // test iteration over the non-zero elements of a sparse matrix
  double mat_sum(matrix<double> m) {

    // for algorithms like sum, must call compact first to ensure that each
    // elt of the matrix is only defined once
    m.compact();

    // outer loop iterate over the columns
    double s = 0.0;
    for (std::size_t j=0; j < m.n_cols(); j++) {

      // inner loop iterate over the non-zero elements of the column
      for (matrix<double>::row_iterator p=m.col(j).begin(); p != m.col(j).end(); ++p)
      {
        //std::size_t i = p.index();
        s += *p;
      }
    }
    return s;
  }

}} // namespace dials_scratch::examples

#endif // DIALS_SCRATCH_EXAMPLES_USE_OF_SPARSE_H
