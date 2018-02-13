from __future__ import absolute_import, division, print_function

import libtbx.pkg_utils

libtbx.pkg_utils.define_entry_points({
  'dxtbx.scaling_model_ext': [
  'aimless = dials_scratch.jbe.scaling_code.model.scaling_model_ext:AimlessScalingModelExt',
  'KB = dials_scratch.jbe.scaling_code.model.scaling_model_ext:KBScalingModelExt',
  'xscale = dials_scratch.jbe.scaling_code.model.scaling_model_ext:XscaleScalingModelExt'
  ]})
