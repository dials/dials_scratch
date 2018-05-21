""""
Classes that each define a smoothly varying component of a scaling model.

These classes use a gaussian smoother (1D, 2D or 3D) to calculate the
inverse scale factors and derivatives with respect to the component
parameters.
"""
import numpy as np
from scitbx import sparse
from dials.array_family import flex
from dials_scaling_helpers_ext import row_multiply
from dials_refinement_helpers_ext import GaussianSmoother as GS1D
from dials_refinement_helpers_ext import GaussianSmoother2D as GS2D
from dials_refinement_helpers_ext import GaussianSmoother3D as GS3D
from dials_scratch_scaling_ext import elementwise_square
from dials_scratch.jbe.scaling_code.model.components.scale_components import \
  ScaleComponentBase

# The following gaussian smoother classes make the implementation
# consistent with that used in dials.refinement.

class GaussianSmoother1D(GS1D):
  """A 1D Gaussian smoother."""

  def value_weight(self, x, value):
    """Return the value, weight and sumweight at a single point."""
    result = super(GaussianSmoother1D, self).value_weight(x,
      flex.double(value))
    return (result.get_value(), result.get_weight(), result.get_sumweight())

  def multi_value_weight(self, x, value):
    """Return the value, weight and sumweight at multiple points."""
    result = super(GaussianSmoother1D, self).multi_value_weight(
      flex.double(x),
      flex.double(value))
    return (result.get_value(), result.get_weight(), result.get_sumweight())

  def positions(self):
    """Return the smoother positions."""
    return list(super(GaussianSmoother1D, self).positions())

class GaussianSmoother2D(GS2D):
  """A 2D Gaussian smoother."""

  def value_weight(self, x, y, value):
    """Return the value, weight and sumweight at a single point."""
    result = super(GaussianSmoother2D, self).value_weight(x, y,
      flex.double(value))
    return (result.get_value(), result.get_weight(), result.get_sumweight())

  def multi_value_weight(self, x, y, value):
    """Return the value, weight and sumweight at multiple points."""
    result = super(GaussianSmoother2D, self).multi_value_weight(
      flex.double(x), flex.double(y),
      flex.double(value))
    return (result.get_value(), result.get_weight(), result.get_sumweight())

  def x_positions(self):
    """Return the smoother x-positions."""
    return list(super(GaussianSmoother2D, self).x_positions())

  def y_positions(self):
    """Return the smoother y-positions."""
    return list(super(GaussianSmoother2D, self).y_positions())

class GaussianSmoother3D(GS3D):
  """A 3D Gaussian smoother."""

  def value_weight(self, x, y, z, value):
    """Return the value, weight and sumweight at a single point."""
    result = super(GaussianSmoother3D, self).value_weight(x, y, z,
      flex.double(value))
    return (result.get_value(), result.get_weight(), result.get_sumweight())

  def multi_value_weight(self, x, y, z, value):
    """Return the value, weight and sumweight at multiple points."""
    result = super(GaussianSmoother3D, self).multi_value_weight(
      flex.double(x), flex.double(y), flex.double(z),
      flex.double(value))
    return (result.get_value(), result.get_weight(), result.get_sumweight())

  def x_positions(self):
    """Return the smoother x-positions."""
    return list(super(GaussianSmoother3D, self).x_positions())

  def y_positions(self):
    """Return the smoother y-positions."""
    return list(super(GaussianSmoother3D, self).y_positions())

  def z_positions(self):
    """Return the smoother z-positions."""
    return list(super(GaussianSmoother3D, self).z_positions())


class SmoothMixin(object):
  """Mixin class for smooth scale factor components.

  This uses a Gaussian smoother to calculate scales and derivatives
  based on the parameters and a have a set of normalised_values
  associated with the data."""

  def __init__(self):
    self._Vr = 1.0
    self._smoother = None

  @property
  def value(self):
    """Extra access to the parameters for the gaussian smoother."""
    return self._parameters

  @property
  def smoother(self):
    """The Gaussian smoother."""
    return self._smoother

  @staticmethod
  def nparam_to_val(n_params):
    """Convert the number of parameters to the required input value
    for the smoother."""
    assert n_params >= 2, '''cannot initialise a smooth scale factor
      for <2 parameters.'''
    if n_params == 2 or n_params == 3:
      return n_params - 1
    return n_params - 2

class SmoothScaleComponent1D(ScaleComponentBase, SmoothMixin):
  """A smoothly varying scale component in one dimension."""

  def __init__(self, initial_values, col_name, parameter_esds=None):
    super(SmoothScaleComponent1D, self).__init__(initial_values,
      parameter_esds)
    self._normalised_values = []
    self._col_name = col_name

  @property
  def normalised_values(self):
    """This is a list of the relevant data needed to calculate the
    inverse scale factors, normalised to give 'normalised coordinate
    values' that fit in the range of the smoother parameters, which
    are defined as a 1D array at normalised coordinates separated by
    a spacing of 1."""
    return self._normalised_values

  @property
  def col_name(self):
    """The column name to use to obtain normalised coordinates from a
    reflection table."""
    return self._col_name

  def update_reflection_data(self, reflection_table, selection=None,
    block_selections=None):
    """Set the normalised coordinate values and configure the smoother."""
    self._normalised_values = []
    self._inverse_scales = []
    self._n_refl = []
    normalised_values = reflection_table[self._col_name]
    if selection:
      normalised_values = normalised_values.select(selection)
    # Make sure zeroed correctly.
    normalised_values = normalised_values - min(normalised_values)
    phi_range_deg = [int(min(normalised_values)//1),
                     int(max(normalised_values)//1)+1]
    self._smoother = GaussianSmoother1D(phi_range_deg,
      self.nparam_to_val(self._n_params))
    if block_selections:
      block_selection_list = block_selections
      #norm_vals_permuted = normalised_values.select(permuted)
      for i, sel in enumerate(block_selection_list):
        #self._normalised_values.append(norm_vals_permuted.select(sel))
        self._normalised_values.append(normalised_values.select(sel))
        self._inverse_scales.append(flex.double(
          self._normalised_values[i].size(), 1.0))
        self._n_refl.append(self.inverse_scales[i].size())
    else:
      self._normalised_values.append(normalised_values)
      self._inverse_scales.append(flex.double(normalised_values.size(), 1.0))
      self._n_refl.append(self._inverse_scales[0].size())

  def calculate_scales_and_derivatives(self, curvatures=False):
    self._inverse_scales = []
    self._derivatives = []
    self._curvatures = []
    for block_id in range(len(self._n_refl)):#len of the list, not num of refl
      value, weight, sumweight = self._smoother.multi_value_weight(
        self._normalised_values[block_id], self.value)
      inv_sw = 1. / sumweight
      dv_dp = row_multiply(weight, inv_sw)
      self._inverse_scales.append(value)
      self._derivatives.append(dv_dp)
      if curvatures:
        self._curvatures.append(sparse.matrix(
          self._inverse_scales[block_id].size(), self.n_params))

  def calculate_scales(self):
    """"Only calculate the scales if needed, for performance."""
    self._inverse_scales = []
    for block_id in range(len(self._n_refl)):#len of the list, not numof refl
      value, _, _ = self._smoother.multi_value_weight(
        self._normalised_values[block_id], self.value)
      self._inverse_scales.append(value)


class SmoothBScaleComponent1D(SmoothScaleComponent1D):
  '''Subclass of SmoothScaleComponent1D to implement a smoothly
  varying B-factor correction.'''

  def __init__(self, initial_values, col_name, parameter_esds=None):
    super(SmoothBScaleComponent1D, self).__init__(initial_values, col_name,
      parameter_esds)
    self._d_values = []

  @property
  def d_values(self):
    """The current set of d-values associated with this component."""
    return self._d_values

  def update_reflection_data(self, reflection_table, selection=None,
    block_selections=None):
    super(SmoothBScaleComponent1D, self).update_reflection_data(
      reflection_table, selection, block_selections)
    #self._d_values = reflection_table['d']
    self._d_values = []
    if selection:
      d_values = reflection_table['d'].select(selection)
    else:
      d_values = reflection_table['d']
    if block_selections:
      block_selection_list = block_selections
      #perumted_d_values = d_values.select(perumted)
      for sel in block_selection_list:
        #self._d_values.append(perumted_d_values.select(sel))
        self._d_values.append(d_values.select(sel))
    else:
      self._d_values.append(d_values)

  def calculate_scales_and_derivatives(self, curvatures=False):
    super(SmoothBScaleComponent1D, self).calculate_scales_and_derivatives(
      curvatures=curvatures)
    for block_id in range(len(self._n_refl)):#len of the list, not numb of refl
      self._inverse_scales[block_id] = flex.double(np.exp(
        self._inverse_scales[block_id] /(2.0 * (self._d_values[block_id]**2))))
      self._derivatives[block_id] = row_multiply(self._derivatives[block_id],
        self._inverse_scales[block_id] / (2.0 * (self._d_values[block_id]**2)))
      if curvatures:
        self._curvatures[block_id] = row_multiply(elementwise_square(
          self._derivatives[block_id]), 1.0/self._inverse_scales[block_id])

  def calculate_scales(self):
    super(SmoothBScaleComponent1D, self).calculate_scales()
    for block_id in range(len(self._n_refl)):#len of the list, not numb of refl
      self._inverse_scales[block_id] = flex.double(np.exp(
        self._inverse_scales[block_id] /(2.0 * (self._d_values[block_id]**2))))


class SmoothScaleComponent2D(ScaleComponentBase, SmoothMixin):
  """Implementation of a 2D array-based smoothly varying scale factor.

  A 2d array of parameters is defined, and the scale factor at fractional
  coordinates is calculated as smoothly varying based on the distance to
  the nearby parameters as calculated in the GaussianSmoother2D. The
  initial values are passed as a 1D array, and shape is a 2-tuple
  indicating the number of parameters in each dimension."""

  def __init__(self, initial_values, shape, col_names, parameter_esds=None):
    assert len(initial_values) == (shape[0] * shape[1]), '''The shape
    information to initialise a 2D smoother is inconsistent with the length
    of the initial parameter list.'''
    super(SmoothScaleComponent2D, self).__init__(initial_values, parameter_esds)
    self._n_x_params = shape[0]
    self._n_y_params = shape[1]
    self._col_names = col_names
    self._normalised_x_values = None
    self._normalised_y_values = None

  @property
  def col_names(self):
    """The column names used to obtain normalised coordinates from a
    reflection table."""
    return self._col_names

  @property
  def n_x_params(self):
    """The number of parameters that parameterise the x-component."""
    return self._n_x_params

  @property
  def n_y_params(self):
    """The number of parameters that parameterise the y-component."""
    return self._n_y_params

  @property
  def normalised_x_values(self):
    """The normalised coordinate values in the first dimension."""
    return self._normalised_x_values

  @property
  def normalised_y_values(self):
    """The normalised coordinate values in the second dimension."""
    return self._normalised_y_values

  def update_reflection_data(self, reflection_table, selection=None,
    block_selections=None):
    '''control access to setting all of reflection data at once'''

    self._normalised_x_values = []
    self._normalised_y_values = []
    self._inverse_scales = []
    self._n_refl = []
    if selection:
      reflection_table = reflection_table.select(selection)
    normalised_x_values = reflection_table[self._col_names[0]]
    normalised_y_values = reflection_table[self._col_names[1]]
    normalised_x_values = normalised_x_values - min(normalised_x_values)
    normalised_y_values = normalised_y_values - min(normalised_y_values)
    x_range = [int(min(normalised_x_values)//1),
               int(max(normalised_x_values)//1)+1]
    y_range = [int(min(normalised_y_values)//1),
               int(max(normalised_y_values)//1)+1]
    self._smoother = GaussianSmoother2D(x_range, self.nparam_to_val(
      self._n_x_params), y_range, self.nparam_to_val(self._n_y_params))
    if block_selections:
      for i, sel in enumerate(block_selections):
        self._normalised_x_values.append(normalised_x_values.select(sel))
        self._normalised_y_values.append(normalised_y_values.select(sel))
        self._inverse_scales.append(flex.double(
          self._normalised_x_values[i].size(), 1.0))
        self._n_refl.append(self.inverse_scales[i].size())
    else:
      self._normalised_x_values.append(normalised_x_values)
      self._normalised_y_values.append(normalised_y_values)
      self._inverse_scales.append(flex.double(normalised_x_values.size(), 1.0))
      self._n_refl.append(self._inverse_scales[0].size())

  def calculate_scales_and_derivatives(self, curvatures=False):
    self._inverse_scales = []
    self._derivatives = []
    self._curvatures = []
    for block_id in range(len(self._n_refl)):
      value, weight, sumweight = self._smoother.multi_value_weight(
        self._normalised_x_values[block_id],
        self._normalised_y_values[block_id], self.value)
      inv_sw = 1. / sumweight
      dv_dp = row_multiply(weight, inv_sw)
      self._inverse_scales.append(value)
      self._derivatives.append(dv_dp)
      if curvatures:
        self._curvatures.append(sparse.matrix(
          self._inverse_scales[block_id].size(), self.n_params))

  def calculate_scales(self):
    """"Only calculate the scales if needed, for performance."""
    self._inverse_scales = []
    for block_id in range(len(self._n_refl)):
      value, _, _ = self._smoother.multi_value_weight(
        self._normalised_x_values[block_id],
        self._normalised_y_values[block_id], self.value)
      self._inverse_scales.append(value)


class SmoothScaleComponent3D(ScaleComponentBase, SmoothMixin):
  """Implementation of a 3D array-based smoothly varying scale factor.

  A 3d array of parameters is defined, and the scale factor at fractional
  coordinates is calculated as smoothly varying based on the distance to
  the nearby parameters as calculated in the GaussianSmoother3D. The
  initial values are passed as a 1D array, and shape is a 3-tuple
  indicating the number of parameters in each dimension."""

  def __init__(self, initial_values, shape, col_names, parameter_esds=None):
    assert len(initial_values) == (shape[0] * shape[1] * shape[2]), '''The
    shape information to initialise a 3D smoother is inconsistent with the
    length of the initial parameter list.'''
    super(SmoothScaleComponent3D, self).__init__(initial_values,
      parameter_esds)
    self._n_x_params = shape[0]
    self._n_y_params = shape[1]
    self._n_z_params = shape[2]
    self._normalised_x_values = None
    self._normalised_y_values = None
    self._normalised_z_values = None
    self._col_names = col_names

  @property
  def col_names(self):
    """The column names used to obtain normalised coordinates from a
    reflection table."""
    return self._col_names

  @property
  def n_x_params(self):
    """The number of parameters that parameterise the x-component."""
    return self._n_x_params

  @property
  def n_y_params(self):
    """The number of parameters that parameterise the y-component."""
    return self._n_y_params

  @property
  def n_z_params(self):
    """The number of parameters that parameterise the z-component."""
    return self._n_z_params

  @property
  def normalised_x_values(self):
    """The normalised coordinate values in the first dimension."""
    return self._normalised_x_values

  @property
  def normalised_y_values(self):
    """The normalised coordinate values in the second dimension."""
    return self._normalised_y_values

  @property
  def normalised_z_values(self):
    """The normalised coordinate values in the third dimension."""
    return self._normalised_z_values

  def update_reflection_data(self, reflection_table, selection=None,
    block_selections=None):
    '''control access to setting all of reflection data at once'''
    self._normalised_x_values = []
    self._normalised_y_values = []
    self._normalised_z_values = []
    self._inverse_scales = []
    self._n_refl = []
    if selection:
      reflection_table = reflection_table.select(selection)
    normalised_x_values = reflection_table[self._col_names[0]]
    normalised_y_values = reflection_table[self._col_names[1]]
    normalised_z_values = reflection_table[self._col_names[2]]
    """Set the normalised coordinate values and configure the smoother."""
    normalised_x_values = normalised_x_values - min(normalised_x_values)
    normalised_y_values = normalised_y_values - min(normalised_y_values)
    normalised_z_values = normalised_z_values - min(normalised_z_values)
    x_range = [int(min(normalised_x_values)//1),
               int(max(normalised_x_values)//1)+1]
    y_range = [int(min(normalised_y_values)//1),
               int(max(normalised_y_values)//1)+1]
    z_range = [int(min(normalised_z_values)//1),
               int(max(normalised_z_values)//1)+1]
    self._smoother = GaussianSmoother3D(x_range, self.nparam_to_val(
      self._n_x_params), y_range, self.nparam_to_val(self._n_y_params),
      z_range, self.nparam_to_val(self._n_z_params))
    if block_selections:
      for i, sel in enumerate(block_selections):
        self._normalised_x_values.append(normalised_x_values.select(sel))
        self._normalised_y_values.append(normalised_y_values.select(sel))
        self._normalised_z_values.append(normalised_z_values.select(sel))
        self._inverse_scales.append(flex.double(
          self._normalised_x_values[i].size(), 1.0))
        self._n_refl.append(self.inverse_scales[i].size())
    else:
      self._normalised_x_values.append(normalised_x_values)
      self._normalised_y_values.append(normalised_y_values)
      self._normalised_z_values.append(normalised_z_values)
      self._inverse_scales.append(flex.double(normalised_x_values.size(), 1.0))
      self._n_refl.append(self._inverse_scales[0].size())


  def calculate_scales_and_derivatives(self, curvatures=False):
    self._inverse_scales = []
    self._derivatives = []
    self._curvatures = []
    for block_id in range(len(self._n_refl)):
      value, weight, sumweight = self._smoother.multi_value_weight(
        self._normalised_x_values[block_id], self._normalised_y_values[block_id],
        self._normalised_z_values[block_id], self.value)
      inv_sw = 1. / sumweight
      dv_dp = row_multiply(weight, inv_sw)
      self._inverse_scales.append(value)
      self._derivatives.append(dv_dp)
      if curvatures:
        self._curvatures.append(sparse.matrix(
          self._inverse_scales[block_id].size(), self.n_params))

  def calculate_scales(self):
    """"Only calculate the scales if needed, for performance."""
    self._inverse_scales = []
    for block_id in range(len(self._n_refl)):
      value, _, _ = self._smoother.multi_value_weight(
        self._normalised_x_values[block_id], self._normalised_y_values[block_id],
        self._normalised_z_values[block_id], self.value)
      self._inverse_scales.append(value)
