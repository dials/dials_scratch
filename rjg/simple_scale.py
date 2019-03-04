from __future__ import division
from __future__ import print_function

import iotbx.phil
from dials.util.options import OptionParser
from cctbx.array_family import flex

help_message = '''
'''

phil_scope = iotbx.phil.parse("""
d_min = None
  .type = float(value_min=0)
""", process_includes=True)


def run(args):
  import libtbx.load_env
  usage = "%s [options]" %libtbx.env.dispatcher_name

  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    check_format=False,
    epilog=help_message)

  params, options, args = parser.parse_args(show_diff_phil=True,
                                            return_unhandled=True)

  assert len(args) == 1
  from iotbx.reflection_file_reader import any_reflection_file

  intensities = None
  batches = None

  f = args[0]

  reader = any_reflection_file(f)
  arrays = reader.as_miller_arrays(merge_equivalents=False)
  for ma in arrays:
    print(ma.info().labels)
    if ma.info().labels == ['I', 'SIGI']:
      intensities = ma
    elif ma.info().labels == ['IMEAN', 'SIGIMEAN']:
      intensities = ma
    elif ma.info().labels == ['I(+)', 'SIGI(+)', 'I(-)', 'SIGI(-)']:
      intensities = ma
    elif ma.info().labels == ['BATCH']:
      batches = ma

  assert intensities is not None
  assert batches is not None

  if params.d_min is not None:
    intensities = intensities.resolution_filter(d_min=params.d_min).set_info(
      intensities.info())

  refined = refinery(intensities, batches)
  refined.plot_scales()

  scaled_unmerged = refined.apply_scales()
  scales = refined.get_scales()

  assert reader.file_type() == 'ccp4_mtz'
  mtz_obj = reader.file_content()

  mtz_dataset = mtz_obj.crystals()[0].datasets()[0]

  labels = intensities.info().labels
  i_column = mtz_dataset.columns()[mtz_dataset.column_labels().index(labels[0])]
  sigi_column = mtz_dataset.columns()[mtz_dataset.column_labels().index(labels[1])]
  selection_valid = flex.bool(scaled_unmerged.size(), True)
  i_column.set_values(
    scaled_unmerged.data().as_float(), selection_valid=selection_valid)
  sigi_column.set_values(
    scaled_unmerged.sigmas().as_float(), selection_valid=selection_valid)

  mtz_dataset.add_column('SCALEUSED', 'R').set_values(scales.as_float())
  mtz_out = 'scaled_unmerged.mtz'
  print('Writing scaled unmerged intensities to %s' %mtz_out)
  mtz_obj.write(mtz_out)
  return


import scitbx.lbfgs

class refinery(object):

  def __init__(self, unmerged_intensities, batches):
    self.unmerged_intensities = unmerged_intensities
    self.batches = batches
    self.minb = flex.min(self.batches.data())
    self.maxb = flex.max(self.batches.data())
    n_scale_factors = self.maxb-self.minb + 1
    self.x = flex.double(n_scale_factors, 1)
    scitbx.lbfgs.run(target_evaluator=self)

  def compute_functional_and_gradients(self):
    unmerged_intensities = self.apply_scales()
    merging = unmerged_intensities.merge_equivalents()
    merged_intensities = merging.array()

    from cctbx import miller
    mi = miller.match_multi_indices(
      merged_intensities.indices(), unmerged_intensities.indices())

    f = 0
    g = flex.double(self.x.size())
    for i in range(unmerged_intensities.size()):
      w = 1/unmerged_intensities.sigmas()[i]
      j = self.batches.data()[i] - self.minb
      k = self.x[j]
      p = mi.pairs()[i]
      mean_I = merged_intensities.data()[p[0]]
      unmerged_I = self.unmerged_intensities.data()[i]
      delta = unmerged_I - k * mean_I
      f += (w * delta**2)
      g[j] += (- 2 * w * mean_I * delta)

    #print f

    return f, g

  def apply_scales(self):
    scales = self.get_scales()
    unmerged_intensities = self.unmerged_intensities.customized_copy(
      data=self.unmerged_intensities.data()*scales)
    return unmerged_intensities

  def get_scales(self):
    scales = flex.double(self.unmerged_intensities.size())
    for i in range(scales.size()):
      scales[i] = self.x[self.batches.data()[i] - self.minb]
    return 1/scales

  def plot_scales(self):
    from matplotlib import pyplot
    scales = 1/self.x
    pyplot.scatter(range(len(scales)), list(scales))
    pyplot.show()

  def callback_after_step(self, minimizer):
    if 0:
      self.plot_scales()
    #print "LBFGS step"


if __name__ == '__main__':
  import sys
  from libtbx.utils import show_times_at_exit
  show_times_at_exit()
  run(sys.argv[1:])
