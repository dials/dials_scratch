from __future__ import division
from scitbx import matrix
from dials.array_family import flex

def rotate_mat3_double(R, A):
  '''
  Helper function to rotate a flex.mat3_double array of matrices

  '''
  accessor = A.accessor()
  RAR = flex.mat3_double([R*matrix.sqr(a)*R.transpose() for a in A])
  RAR.reshape(accessor)
  return RAR


def compute_change_of_basis_operation(s0, s2):
  '''
  Compute the change of basis operation that puts the s2 vector along the z axis

  '''
  e1 = s2.cross(s0).normalize()
  e2 = s2.cross(e1).normalize()
  e3 = s2.normalize()
  R = matrix.sqr(
    e1.elems +
    e2.elems +
    e3.elems)
  return R
