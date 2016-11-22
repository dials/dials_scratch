#!/usr/bin/env cctbx.python

"""Examples of use of the scitbx sparse matrix and vector type"""

from scitbx import sparse
from dials_scratch_cctbx_cpp_examples_ext import create_sparse_matrix

if __name__ == '__main__':

  mat = create_sparse_matrix()
  assert mat.non_zeroes == 1
  assert mat.n_cols, mat.n_rows == (3, 5)
  for i in range(mat.n_rows):
    for j in range(mat.n_cols):
      if not mat.is_structural_zero(i, j):
        assert (i, j) == (0, 0)
