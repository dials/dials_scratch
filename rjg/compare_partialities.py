from __future__ import division
from __future__ import print_function
import libtbx.phil

help_message = '''
'''

phil_scope= libtbx.phil.parse("""
""", process_includes=True)

import matplotlib
matplotlib.use('Agg')

def run(args):

  from cctbx.array_family import flex
  from dials.util.options import OptionParser
  from dials.util.options import flatten_reflections
  import libtbx.load_env

  usage = "%s [options] reflections_1.pickle reflections_2.pickle" %(
    libtbx.env.dispatcher_name)

  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    read_reflections=True,
    epilog=help_message)

  params, options, args = parser.parse_args(
    show_diff_phil=True, return_unhandled=True)
  reflections = flatten_reflections(params.input.reflections)

  if flex.max(reflections[0]['id']) > 0:
    reflections = list(reversed(reflections))
  assert flex.max(reflections[0]['id']) == 0

  assert len(reflections) == 2
  partialities = []
  intensities = []
  sigmas = []
  ids = []
  xyz = []

  # only want fully-recorded reflections in full dataset
  #reflections[0] = reflections[0].select(reflections[0]['partiality'] > 0.99)
  print(reflections[0].size())
  # only want partial reflections in sliced dataset
  #reflections[1] = reflections[1].select(reflections[1]['partiality'] < 0.99)
  print(reflections[1].size())

  for refl in reflections:
    #sel = refl.get_flags(refl.flags.integrated_sum)
    sel = refl.get_flags(refl.flags.integrated)
    sel &= refl['intensity.sum.value'] > 0
    sel &= refl['intensity.sum.variance'] > 0
    refl = refl.select(sel)
    hkl = refl['miller_index']
    partiality = refl['partiality']
    intensity = refl['intensity.sum.value']
    vari = refl['intensity.sum.variance']
    assert vari.all_gt(0)
    sigi = flex.sqrt(vari)
    intensities.append(intensity)
    partialities.append(partiality)
    sigmas.append(sigi)
    ids.append(refl['id'])
    xyz.append(refl['xyzcal.px'])

  from annlib_ext import AnnAdaptor as ann_adaptor
  ann = ann_adaptor(xyz[0].as_double().as_1d(), 3)
  ann.query(xyz[1].as_double().as_1d())
  distances = flex.sqrt(ann.distances)
  matches = distances < 2 # pixels
  isel0 = flex.size_t(list(ann.nn.select(matches)))
  isel1 = flex.size_t(list(matches.iselection()))

  p0 = partialities[0].select(isel0)
  p1 = partialities[1].select(isel1)
  i0 = intensities[0].select(isel0)
  i1 = intensities[1].select(isel1)

  print((p0 > p1).count(True), (p0 < p1).count(True))

  h0 = flex.histogram(p0, data_min=0, data_max=1, n_slots=20)
  h1 = flex.histogram(p1, data_min=0, data_max=1, n_slots=20)
  h0.show()
  h1.show()

  from matplotlib import pyplot

  perm0 = flex.sort_permutation(p0)
  perm1 = flex.sort_permutation(p1)
  fig, axes = pyplot.subplots(nrows=2, sharex=True)
  axes[0].plot(p0.select(perm0), flex.int_range(p0.size()))
  axes[1].plot(p1.select(perm1), flex.int_range(p1.size()))
  axes[1].set_xlabel('Partiality')
  for ax in axes: ax.set_ylabel('Cumulative frequency')
  for ax in axes: ax.set_yscale('log')
  pyplot.savefig('sorted_partialities.png')
  pyplot.clf()

  blue = '#3498db'
  fig, axes = pyplot.subplots(nrows=2, sharex=True)
  axes[0].bar(h0.slot_centers(), h0.slots(), width=h0.slot_width(),
              align='center', color=blue, edgecolor=blue)
  axes[1].bar(h1.slot_centers(), h1.slots(), width=h1.slot_width(),
              align='center', color=blue, edgecolor=blue)
  axes[1].set_xlabel('Partiality')
  for ax in axes: ax.set_ylabel('Frequency')
  for ax in axes: ax.set_yscale('log')
  pyplot.savefig('partiality_histogram.png')
  #pyplot.show()
  pyplot.clf()

  pyplot.scatter(p0, p1, s=5, alpha=0.3, marker='+')
  pyplot.xlabel('Partiality (full)')
  pyplot.ylabel('Partiality (sliced)')
  pyplot.savefig('partiality_full_vs_sliced.png')
  pyplot.clf()

  pyplot.scatter(i0, i1, s=5, alpha=0.3, marker='+')
  pyplot.xlim(flex.min(i0), flex.max(i0))
  pyplot.ylim(flex.min(i1), flex.max(i1))
  pyplot.xlabel('Intensity (full)')
  pyplot.ylabel('Intensity (sliced)')
  pyplot.xscale('log')
  pyplot.yscale('log')
  pyplot.savefig('intensity_full_vs_sliced.png')
  pyplot.clf()

  i_ratio = i1/i0
  p_ratio = p1/p0
  pyplot.scatter(p_ratio, i_ratio, s=5, alpha=0.3, marker='+')
  pyplot.ylim(flex.min(i_ratio), flex.max(i_ratio))
  pyplot.yscale('log')
  pyplot.xlabel('P(full)/P(sliced)')
  pyplot.ylabel('I(full)/I(sliced)')
  pyplot.savefig('partiality_ratio_vs_intensity_ratio.png')
  pyplot.clf()

if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
