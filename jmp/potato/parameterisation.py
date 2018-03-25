from __future__ import division
from scitbx import matrix
from dials.array_family import flex
from dials.algorithms.refinement.parameterisation.crystal_parameters import CrystalUnitCellParameterisation
from dials.algorithms.refinement.parameterisation.crystal_parameters import CrystalOrientationParameterisation


class SimpleMosaicityParameterisation(object):
  '''
  A simple mosaicity parameterisation that uses 6 parameters to describe a
  multivariate normal reciprocal lattice profile. Sigma is enforced as positive
  definite by parameterising using the cholesky decomposition.

  M = | b1 0  0  |
      | b2 b3 0  |
      | b4 b5 b6 |

  S = M*M^T

  '''

  def __init__(self, params=None):
    '''
    Initialise with the parameters

    '''
    if params is not None:
      assert len(params) == self.num_parameters()
      self.params = params
    else:
      self.params = flex.double(self.num_parameters(), 0)

  def num_parameters(self):
    '''
    Get the number of parameters

    '''
    return 6

  def set_parameters(self, params):
    '''
    Set the parameters

    '''
    assert len(params) == self.num_parameters()
    self.params = params

  def parameters(self):
    '''
    Return the parameters

    '''
    return self.params

  def sigma(self):
    '''
    Compute the covariance matrix of the MVN from the parameters

    '''
    M = matrix.sqr((
      self.params[0], 0, 0,
      self.params[1], self.params[2], 0,
      self.params[3], self.params[4], self.params[5]))
    return M*M.transpose()

  def first_derivatives(self):
    '''
    Compute the first derivatives of Sigma w.r.t the parameters

    '''
    b1, b2, b3, b4, b5, b6 = self.params

    dSdb1 = (
      2*b1,b2,b4,
      b2,0,0,
      b4,0,0)

    dSdb2 = (
      0,b1,0,
      b1,2*b2,b4,
      0,b4,0)

    dSdb3 = (
      0,0,0,
      0,2*b3,b5,
      0,b5,0)

    dSdb4 = (
      0,0,b1,
      0,0,b2,
      b1,b2,2*b4)

    dSdb5 = (
      0,0,0,
      0,0,b3,
      0,b3,2*b5)

    dSdb6 = (
      0,0,0,
      0,0,0,
      0,0,2*b6)

    return flex.mat3_double([dSdb1, dSdb2, dSdb3, dSdb4, dSdb5, dSdb6])

  def second_derivatives(self):
    '''
    Compute the second derivatives of Sigma w.r.t the parameters

    '''
    b1, b2, b3, b4, b5, b6 = self.params

    zero = (
      0, 0, 0,
      0, 0, 0,
      0, 0, 0)

    d11 = (
      2, 0, 0,
      0, 0, 0,
      0, 0, 0)
    d12 = (
      0, 1, 0,
      1, 0, 0,
      0, 0, 0)
    d13 = zero
    d14 = (
      0, 0, 1,
      0, 0, 0,
      1, 0, 0)
    d15 = zero
    d16 = zero

    d21 = d12
    d22 = (
      0, 0, 0,
      0, 2, 0,
      0, 0, 0)
    d23 = zero
    d24 = (
      0, 0, 0,
      0, 0, 1,
      0, 1, 0)
    d25 = zero
    d26 = zero

    d31 = zero
    d32 = zero
    d33 = (
      0, 0, 0,
      0, 2, 0,
      0, 0, 0)
    d34 = zero
    d35 = (
      0, 0, 0,
      0, 0, 1,
      0, 1, 0)
    d36 = zero

    d41 = d14
    d42 = d24
    d43 = zero
    d44 = (
      0, 0, 0,
      0, 0, 0,
      0, 0, 2)
    d45 = zero
    d46 = zero

    d51 = zero
    d52 = zero
    d53 = d35
    d54 = zero
    d55 = (
      0, 0, 0,
      0, 0, 0,
      0, 0, 2)
    d56 = zero

    d61 = zero
    d62 = zero
    d63 = zero
    d64 = zero
    d65 = zero
    d66 = (
      0, 0, 0,
      0, 0, 0,
      0, 0, 2)

    d2 = flex.mat3_double([
      d11, d12, d13, d14, d15, d16,
      d21, d22, d23, d24, d25, d26,
      d31, d32, d33, d34, d35, d36,
      d41, d42, d43, d44, d45, d46,
      d51, d52, d53, d54, d55, d56,
      d61, d62, d63, d64, d65, d66])
    d2.reshape(flex.grid(6,6))

    return d2


class WavelengthSpreadParameterisation(object):
  '''
  A simple mosaicity parameterisation that uses 1 parameters to describe a
  multivariate normal wavelength spread. Sigma is enforced as positive
  definite by parameterising using the cholesky decomposition.

  L = | 0 0 0  |
      | 0 0 0  |
      | 0 0 l1 |

  S = L*L^T

  '''

  def __init__(self, params=None):
    '''
    Initialise with the parameters

    '''
    if params is not None:
      assert len(params) == self.num_parameters()
      self.params = params
    else:
      self.params = flex.double(self.num_parameters(), 0)

  def num_parameters(self):
    '''
    Get the number of parameters

    '''
    return 1

  def set_parameters(self, params):
    '''
    Set the parameters

    '''
    assert len(params) == self.num_parameters()
    self.params = params

  def parameters(self):
    '''
    Return the parameters

    '''
    return self.params

  def sigma(self):
    '''
    Compute the covariance matrix of the MVN from the parameters

    '''
    M = matrix.sqr((
      0, 0, 0,
      0, 0, 0,
      0, 0, self.params[0]))
    return M*M.transpose()

  def first_derivatives(self):
    '''
    Compute the first derivatives of Sigma w.r.t the parameters

    '''
    l1 = self.params[0]

    dSdl1 = (
      0,0,0,
      0,0,0,
      0,0,2*l1)

    return flex.mat3_double([dSdl1])

  def second_derivatives(self):
    '''
    Compute the second derivatives of Sigma w.r.t the parameters

    '''
    l1 = self.params[0]

    d11 = (
      0,0,0,
      0,0,0,
      0,0,2)

    d2 = flex.mat3_double([d11])
    d2.reshape(flex.grid(1,1))

    return d2


class AngularMosaicityParameterisation(object):
  '''
  A simple mosaicity parameterisation that uses 3 parameters to describe a
  multivariate normal angular mosaic spread. Sigma is enforced as positive
  definite by parameterising using the cholesky decomposition.

  W = | w1 0  0  |
      | w2 w3 0  |
      | 0  0  0 |

  S = W*W^T

  '''

  def __init__(self, params=None):
    '''
    Initialise with the parameters

    '''
    if params is not None:
      assert len(params) == self.num_parameters()
      self.params = params
    else:
      self.params = flex.double(self.num_parameters(), 0)

  def num_parameters(self):
    '''
    Get the number of parameters

    '''
    return 3

  def set_parameters(self, params):
    '''
    Set the parameters

    '''
    assert len(params) == self.num_parameters()
    self.params = params

  def parameters(self):
    '''
    Return the parameters

    '''
    return self.params

  def sigma(self):
    '''
    Compute the covariance matrix of the MVN from the parameters

    '''
    M = matrix.sqr((
      self.params[0], 0, 0,
      self.params[1], self.params[2], 0,
      0, 0, 0))
    return M*M.transpose()

  def first_derivatives(self):
    '''
    Compute the first derivatives of Sigma w.r.t the parameters

    '''
    b1, b2, b3 = self.params

    d1 = (
      2*b1,b2,0,
      b2,0,0,
      0,0,0)

    d2 = (
      0,b1,0,
      b1,2*b2,0,
      0,0,0)

    d3 = (
      0,0,0,
      0,2*b3,0,
      0,0,0)

    return flex.mat3_double([d1, d2, d3])

  def second_derivatives(self):
    '''
    Compute the second derivatives of Sigma w.r.t the parameters

    '''
    b1, b2, b3 = self.params

    zero = (
      0, 0, 0,
      0, 0, 0,
      0, 0, 0)

    d11 = (
      2, 0, 0,
      0, 0, 0,
      0, 0, 0)
    d12 = (
      0, 1, 0,
      1, 0, 0,
      0, 0, 0)
    d13 = zero

    d21 = d12
    d22 = (
      0, 0, 0,
      0, 2, 0,
      0, 0, 0)
    d23 = zero

    d31 = zero
    d32 = zero
    d33 = (
      0, 0, 0,
      0, 2, 0,
      0, 0, 0)

    d2 = flex.mat3_double([
      d11, d12, d13,
      d21, d22, d23,
      d31, d32, d33])
    d2.reshape(flex.grid(3,3))

    return d2


class ModelState(object):
  '''
  A class to keep track of the model state

  '''

  def __init__(self,
               crystal,
               fix_rlp_mosaicity=False,
               fix_wavelength_spread=False,
               fix_angular_mosaicity=False,
               fix_unit_cell=False,
               fix_orientation=False):
    '''
    Initialise the model state

    '''

    # Save the crystal model
    self.crystal = crystal

    # The U and P parameterisation
    self.U_parameterisation = CrystalOrientationParameterisation(self.crystal)
    self.B_parameterisation = CrystalUnitCellParameterisation(self.crystal)

    # The M, L and W parameterisations
    self.M_parameterisation = SimpleMosaicityParameterisation()
    self.L_parameterisation = WavelengthSpreadParameterisation()
    self.W_parameterisation = AngularMosaicityParameterisation()

    # Set the flags to fix parameters
    self._is_rlp_mosaicity_fixed = fix_rlp_mosaicity
    self._is_wavelength_spread_fixed = fix_wavelength_spread
    self._is_angular_mosaicity_fixed = fix_angular_mosaicity
    self._is_unit_cell_fixed = fix_unit_cell
    self._is_orientation_fixed = fix_orientation

  def is_orientation_fixed(self):
    '''
    Return whether the orientation is fixed

    '''
    return self._is_orientation_fixed

  def is_unit_cell_fixed(self):
    '''
    Return whether the unit celll is fixed

    '''
    return self._is_unit_cell_fixed

  def is_rlp_mosaicity_fixed(self):
    '''
    Return whether the rlp mosaicity is fixed

    '''
    return self._is_rlp_mosaicity_fixed

  def is_wavelength_spread_fixed(self):
    '''
    Return whether the wavelength spread is fixed

    '''
    return self._is_wavelength_spread_fixed

  def is_angular_mosaicity_fixed(self):
    '''
    Return whether the angular mosaicity is fixed

    '''
    return self._is_angular_mosaicity_fixed

  def get_U(self):
    '''
    Get the crystal U matrix

    '''
    return self.crystal.get_U()

  def get_B(self):
    '''
    Get the crystal B matrix

    '''
    return self.crystal.get_B()

  def get_A(self):
    '''
    Get the crystal A matrix

    '''
    return self.crystal.get_A()

  def get_M(self):
    '''
    Get the Sigma M matrix

    '''
    return self.M_parameterisation.sigma()

  def get_L(self):
    '''
    Get the Sigma L matrix

    '''
    return self.L_parameterisation.sigma()

  def get_W(self):
    '''
    Get the Sigma W matrix

    '''
    return self.W_parameterisation.sigma()

  def get_U_params(self):
    '''
    Get the U parameters

    '''
    return flex.double(self.U_parameterisation.get_param_vals())

  def get_B_params(self):
    '''
    Get the B parameters

    '''
    return flex.double(self.B_parameterisation.get_param_vals())

  def get_M_params(self):
    '''
    Get the M parameters

    '''
    return self.M_parameterisation.parameters()

  def get_L_params(self):
    '''
    Get the L parameters

    '''
    return self.L_parameterisation.parameters()

  def get_W_params(self):
    '''
    Get the W parameters

    '''
    return self.W_parameterisation.parameters()

  def set_U_params(self, params):
    '''
    Set the U parameters

    '''
    return self.U_parameterisation.set_param_vals(params)

  def set_B_params(self, params):
    '''
    Set the B parameters

    '''
    return self.B_parameterisation.set_param_vals(params)

  def set_M_params(self, params):
    '''
    Set the M parameters

    '''
    return self.M_parameterisation.set_parameters(params)

  def set_L_params(self, params):
    '''
    Set the L parameters

    '''
    return self.L_parameterisation.set_parameters(params)

  def set_W_params(self, params):
    '''
    Set the W parameters

    '''
    return self.W_parameterisation.set_parameters(params)

  def num_U_params(self):
    '''
    Get the number of U parameters

    '''
    return len(self.get_U_params())

  def num_B_params(self):
    '''
    Get the number of B parameters

    '''
    return len(self.get_B_params())

  def num_M_params(self):
    '''
    Get the number of M parameters

    '''
    return len(self.get_M_params())

  def num_L_params(self):
    '''
    Get the number of L parameters

    '''
    return len(self.get_L_params())

  def num_W_params(self):
    '''
    Get the number of W parameters

    '''
    return len(self.get_W_params())

  def get_dU_dp(self):
    '''
    Get the first derivatives of U w.r.t its parameters

    '''
    return flex.mat3_double(self.U_parameterisation.get_ds_dp())

  def get_dB_dp(self):
    '''
    Get the first derivatives of B w.r.t its parameters

    '''
    return flex.mat3_double(self.B_parameterisation.get_ds_dp())

  def get_dM_dp(self):
    '''
    Get the first derivatives of M w.r.t its parameters

    '''
    return self.M_parameterisation.first_derivatives()

  def get_dL_dp(self):
    '''
    Get the first derivatives of L w.r.t its parameters

    '''
    return self.L_parameterisation.first_derivatives()

  def get_dW_dp(self):
    '''
    Get the first derivatives of W w.r.t its parameters

    '''
    return self.W_parameterisation.first_derivatives()

  def active_parameters(self):
    '''
    Get the active parameters in order: U, B, M, L, W

    '''
    active_params = flex.double()
    if not self._is_orientation_fixed:
      active_params.extend(self.get_U_params())
    if not self._is_unit_cell_fixed:
      active_params.extend(self.get_B_params())
    if not self._is_rlp_mosaicity_fixed:
      active_params.extend(self.get_M_params())
    if not self._is_wavelength_spread_fixed:
      active_params.extend(self.get_L_params())
    if not self._is_angular_mosaicity_fixed:
      active_params.extend(self.get_W_params())
    return active_params

  def set_active_parameters(self, params):
    '''
    Set the active parameters in order: U, B, M, L, W

    '''
    if not self._is_orientation_fixed:
      temp   = params[:self.num_U_params()]
      params = params[self.num_U_params():]
      self.set_U_params(temp)
    if not self._is_unit_cell_fixed:
      temp   = params[:self.num_B_params()]
      params = params[self.num_B_params():]
      self.set_B_params(temp)
    if not self._is_rlp_mosaicity_fixed:
      temp   = params[:self.num_M_params()]
      params = params[self.num_M_params():]
      self.set_M_params(temp)
    if not self._is_wavelength_spread_fixed:
      temp   = params[:self.num_L_params()]
      params = params[self.num_L_params():]
      self.set_L_params(temp)
    if not self._is_angular_mosaicity_fixed:
      temp   = params[:self.num_W_params()]
      self.set_W_params(temp)



class ReflectionModel(object):

  def __init__(self, state, s0, h):

    # Get a load of matrices
    A = state.get_A()
    U = state.get_U()
    B = state.get_B()
    M = state.get_M()
    L = state.get_L()
    W = state.get_W()

    # Compute the reciprocal lattice vector
    self._h = h
    self._r = A*h
    t = self.r.length()

    # Define rotation for L and W sigma components
    q1 = self._r.cross(s0).normalize()
    q2 = self._r.cross(q1).normalize()
    q3 = self._r.normalize()
    Q = matrix.sqr(
      q1.elems +
      q2.elems +
      q3.elems)

    # Compute the covariance matrix
    self._sigma = M + self.r.length()*Q*(L + W)*Q.transpose()

    # Get the derivatives of the model state
    dU_dp = state.get_dU_dp()
    dB_dp = state.get_dB_dp()
    dM_dp = state.get_dM_dp()
    dL_dp = state.get_dL_dp()
    dW_dp = state.get_dW_dp()

    # Compute the derivatives of r w.r.t the parameters
    dr_dp_m = flex.vec3_double((0,0,0)) * len(state.num_M_params())
    dr_dp_l = flex.vec3_double((0,0,0)) * len(state.num_L_params())
    dr_dp_w = flex.vec3_double((0,0,0)) * len(state.num_W_params())
    dr_dp_u = dU_dp*B*h
    dr_dp_b = U*dB_dp*h

    # Derivative of |r| w.r.t U and B parameters
    dt_dp_u = r.dot(dr_dp_u) / t
    dt_dp_b = r.dot(dr_dp_b) / t

    # Derivatives of q1, q2 and q3 w.r.t U and B parameters
    r = self.r
    rs0 = r.cross(s0)
    rq1 = r.cross(q1)
    drs0_dp_u = dr_dp_u.cross(s0)
    drs0_dp_b = dr_dp_b.cross(s0)
    dq1_dp_u = drs0_dp_u/rs0.length() - rs0*rs0.dot(drs0_dp_u)/rs0.length()**3
    dq1_dp_b = drs0_dp_b/rs0.length() - rs0*rs0.dot(drs0_dp_b)/rs0.length()**3
    drq1_dp_u = dr_dp_u.cross(q1) + r.cross(dq1_dp_u)
    drq1_dp_b = dr_dp_b.cross(q1) + r.cross(dq1_dp_b)
    dq2_dp_u = drq1_dp_u/rq1.length() - rq1*rq1.dot(drq1_dp_u)/rq1.length()**3
    dq2_dp_b = drq1_dp_b/rq1.length() - rq1*rq1.dot(drq1_dp_b)/rq1.length()**3
    dq3_dp_u = dr_dp_u/r.length() - r*r.dot(dr_dp_u)/r.length()**3
    dq3_dp_b = dr_dp_b/r.length() - r*r.dot(dr_dp_b)/r.length()**3

    # Derivatives of Q w.r.t U and B parameters
    dQ_dp_u = matrix.sqr(
      dq1_dp_u.elems +
      dq2_dp_u.elems +
      dq3_dp_u.elems)
    dQ_dp_b = matrix.sqr(
      dq1_dp_b.elems +
      dq2_dp_b.elems +
      dq3_dp_b.elems)

    # Compute derivatives of Sigma w.r.t the parameters
    ds_dp_m = dM_dp
    ds_dp_l = self.r.length()*Q*dL_dp*Q.transpose()
    ds_dp_w = self.r.length()*Q*dW_dp*Q.transpose()
    ds_dp_u = dt_dp_u*Q*(L+W)*Q.transpose() \
            + t*dQ_dp_u*(L+W)*Q.transpose() \
            + t*Q*(L+W)*dQ_dp_u.transpose()
    ds_dp_b = dt_dp_b*Q*(L+W)*Q.transpose() \
            + t*dQ_dp_b*(L+W)*Q.transpose() \
            + t*Q*(L+W)*dQ_dp_b.transpose()

    # Create array with derivatives of r
    self._dr_dp = flex.vec3_double()
    self._dr_dp.extend(dr_dp_m)
    self._dr_dp.extend(dr_dp_l)
    self._dr_dp.extend(dr_dp_w)
    self._dr_dp.extend(dr_dp_u)
    self._dr_dp.extend(dr_dp_b)

    # Create array with derivatives of S
    self._ds_dp = flex.mat3_double()
    self._ds_dp.extend(ds_dp_m)
    self._ds_dp.extend(ds_dp_l)
    self._ds_dp.extend(ds_dp_w)
    self._ds_dp.extend(ds_dp_u)
    self._ds_dp.extend(ds_dp_b)

  def sigma(self):
    return self._sigma

  def r(self):
    return self._r

  def get_dS_dp(self):
    pass

  def get_dr_dp(self):
    return self._dr_dp
