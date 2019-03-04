"""
Usage:
  compare_aimless_scales integrated_scaled.pickle scaled_unmerged.mtz
"""
from __future__ import print_function

import sys
import math
from cctbx.array_family import flex
from annlib_ext import AnnAdaptor as ann_adaptor
from dials.util.options import OptionParser, flatten_reflections
from libtbx import phil
import cPickle as pickle

def load_data(filename):
  data_file = open(filename)
  data = pickle.load(data_file)
  data_file.close()
  return data

def meansd(values):
    assert(len(values) > 3)

    mean = sum(values) / len(values)
    var = sum([(v - mean) * (v - mean) for v in values]) / (len(values) - 1)

    return mean, math.sqrt(var)

def cc(a, b):

    assert(len(a) == len(b))

    ma, sa = meansd(a)
    mb, sb = meansd(b)

    r = (1 / (len(a) - 1)) * sum([((a[j] - ma) / sa) * ((b[j] - mb) / sb)
                                  for j in range(len(a))])

    return r

def cc_slope(a, b):

    assert(len(a) == len(b))
    fa = flex.double(a)
    fb = flex.double(b)
    aa = fa * fa
    bb = fb * fb
    ab = fa * fb
    N = len(a)
    sum_ab = flex.sum(ab)
    sum_a = flex.sum(fa)
    sum_b = flex.sum(fb)
    sum_aa = flex.sum(aa)
    sum_bb = flex.sum(bb)
    slope = (N * sum_ab - sum_a * sum_b) / (N * sum_aa - sum_a**2)
    corr  = (N * sum_ab - sum_a * sum_b) / (math.sqrt(N * sum_aa - sum_a**2) *
             math.sqrt(N * sum_bb - sum_b**2))

    return corr,slope

def R(calc, obs):

    assert(len(calc) == len(obs))

    scale = sum(obs) / sum(calc)
    print("scale difference",scale)

    return sum([math.fabs(math.fabs(o) - math.fabs(scale * c)) \
                for c, o in zip(calc, obs)]) / \
                sum([math.fabs(o) for o in obs])

def read_aimless(unmerged_mtz):
  #mtzfile = aimless_scaled_file
  from iotbx.reflection_file_reader import any_reflection_file

  reader = any_reflection_file(unmerged_mtz)
  assert reader.file_type() == 'ccp4_mtz'
  arrays = reader.as_miller_arrays(merge_equivalents=False)

  intensities = None
  scales = None

  for ma in arrays:
    print(ma.info().labels)
    if ma.info().labels == ['I', 'SIGI']:
      intensities = ma
    elif ma.info().labels == ['SCALEUSED']:
      scales = ma
    elif ma.info().labels == ['XDET']:
      xdet = ma
    elif ma.info().labels == ['YDET']:
      ydet = ma
    elif ma.info().labels == ['ROT']:
      rot = ma
    elif ma.info().labels == ['LP']:
      lp = ma

  assert intensities is not None
  assert scales is not None

  indices = reader.file_content().extract_original_index_miller_indices()
  #indices = intensities.indices()
  observations = []

  for hkl, x,y,z, Isigma, LP, scale in zip(indices, xdet.data(), ydet.data(), rot.data(),
  intensities, lp.data(), scales.data()):
    observations.append((hkl, (x,y,z), [Isigma[1:2][0], Isigma[2:][0]], 1/scale))
  return observations


def read_dials_scaled(reflections):

  hkls = reflections['miller_index']
  xs = reflections['xyzobs.px.value'].parts()[0]
  ys = reflections['xyzobs.px.value'].parts()[1]
  zs = reflections['phi']
  Is = (reflections['intensity.prf.value'] * reflections['lp']
  / (reflections['dqe'] * reflections['inverse_scale_factor']))
  sigmas = ((reflections['intensity.prf.variance']**0.5) * reflections['lp']
  / (reflections['dqe'] * reflections['inverse_scale_factor']))
  scales = reflections['inverse_scale_factor']

  observations = []
  for hkl, x, y, z, I, sigma, scale in zip(hkls, xs, ys, zs, Is, sigmas, scales):
    observations.append((hkl, (x,y,z), [I, sigma], scale))
  return observations


def main(args):
  reflections = load_data(args[0])
  aimless_scaled_file = args[1]

  dials_hkl_xyz_isigi = read_dials_scaled(reflections)

  print('Read %d observations from %s' % (len(dials_hkl_xyz_isigi), args[0]))

  aimless_hkl_xyz_isigi = read_aimless(aimless_scaled_file)

  print('Read %d observations from %s' % \
        (len(aimless_hkl_xyz_isigi), aimless_scaled_file))

  # treat aimless as reference, dials as query
  reference = flex.double()
  query = flex.double()

  for hkl, xyz, isigi, scale in dials_hkl_xyz_isigi:
    query.append(xyz[0])
    query.append(xyz[1])
    query.append(xyz[2])
  #print hkl, xyz, isigi

  #extract xyz positions
  for hkl, xyz, isigi, scale in aimless_hkl_xyz_isigi:
    reference.append(xyz[0])
    reference.append(xyz[1])
    reference.append(xyz[2])
  #print hkl, xyz, isigi

  ann = ann_adaptor(data = reference, dim = 3, k = 1)
  ann.query(query)

  i_s_dials = []
  i_s_aimless = []

  #go through matching the reflections, then appending the (scaled) I's to a list
  for j in range(len(dials_hkl_xyz_isigi)):
    c = ann.nn[j]
    if aimless_hkl_xyz_isigi[c][0] == dials_hkl_xyz_isigi[j][0]:
      #print aimless_hkl_xyz_isigi[c][0], dials_hkl_xyz_isigi[j][0]
      i_s_dials.append(dials_hkl_xyz_isigi[j][3])
      i_s_aimless.append(aimless_hkl_xyz_isigi[c][3])
    #else:
      #print aimless_hkl_xyz_isigi[c][0], dials_hkl_xyz_isigi[j][0]

  #calculate correlation between I's
  print('Matched %d observations' % len(i_s_dials))
  correlation_coefficient = cc(i_s_dials, i_s_aimless)
  R_factor = R(i_s_dials, i_s_aimless)
  print('CC: %.6f' % correlation_coefficient)
  print('R:  %.3f' % R_factor)
  import matplotlib.pyplot as plt
  y_ideal = [x for x in i_s_dials]


  plt.figure(figsize=(10,7))
  plt.scatter(i_s_dials, i_s_aimless, s=0.1)
  plt.plot(i_s_dials, y_ideal, color='r')
  plt.xlabel('Inverse scale factor in DIALS')
  plt.ylabel('Inverse scale factor in aimless')
  plt.axes().set_aspect('equal')
  plt.title('''Comparison inverse scale factors from aimless and
    dials.aimless_scaling, CC = %.5f, R = %.3f''' % (correlation_coefficient, R_factor))
  plt.savefig('Aimless_DIALS_comparison.png')
  plt.show()
  return correlation_coefficient, R_factor

if __name__=="__main__":
  #this expects two arguments - a integrated_scaled pickle file from DIALS
  #and a scaled.mtz file from aimless.
  import sys
  args = sys.argv[1:]
  main(args)
