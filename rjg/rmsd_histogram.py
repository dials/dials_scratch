from __future__ import division

import iotbx.phil
from dials.util.options import OptionParser
from dials.util.options import flatten_reflections
from dials.array_family import flex
import math

help_message = '''
'''

phil_scope = iotbx.phil.parse('''
min_reflections = 1
  .type = int(value_min=1)
''', process_includes=True)


def run(args):
  import libtbx.load_env
  usage = "%s reflections.pickle [options]" %libtbx.env.dispatcher_name

  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    read_reflections=True,
    epilog=help_message)

  params, options = parser.parse_args(show_diff_phil=True)
  reflections = flatten_reflections(params.input.reflections)
  if len(reflections) == 0:
    parser.print_help()
    return
  reflections = reflections[0]
  reflections = reflections.select(
    reflections.get_flags(reflections.flags.used_in_refinement))

  from dials.algorithms.shoebox import MaskCode
  bg_code = MaskCode.Valid | MaskCode.BackgroundUsed
  fg_code = MaskCode.Valid | MaskCode.Foreground
  strong_code = MaskCode.Valid | MaskCode.Strong

  reflection_size = flex.size_t()
  for shoebox in reflections['shoebox']:
    reflection_size.append(shoebox.mask.count(fg_code))

  ref_by_npix = {}
  for npix in range(1, flex.max(reflection_size)):
    ref_by_npix[npix] = reflections.select(reflection_size == npix)

  n_pixels = flex.size_t()
  rmsds = flex.vec3_double()
  mean_var = flex.vec3_double()
  n_ref = flex.size_t()

  for npix in range(1, flex.max(reflection_size)):
    refl = reflections.select(reflection_size == npix)
    if refl.size() < params.min_reflections: continue
    sel = refl['xyzcal.mm'].norms() > 0
    if sel.count(True) == 0: continue
    refl = refl.select(sel)
    xc, yc, zc = refl['xyzcal.mm'].parts()
    xo, yo, zo = refl['xyzobs.mm.value'].parts()

    rmsdx = math.sqrt(flex.mean(flex.pow2(xc-xo)))
    rmsdy = math.sqrt(flex.mean(flex.pow2(yc-yo)))
    rmsdz = math.sqrt(flex.mean(flex.pow2(zc-zo)))

    #print '%i (%i refs): (%.2f, %.2f, %.2f)' %(
      #npix, refl.size(), rmsdx, rmsdy, rmsdz)

    n_pixels.append(npix)
    rmsds.append((rmsdx, rmsdy, rmsdz))
    n_ref.append(refl.size())
    mean_var.append(refl['xyzobs.mm.variance'].mean())

  from matplotlib import pyplot
  pyplot.style.use('ggplot')
  rmsdx, rmsdy, rmsdz = rmsds.parts()
  varx, vary, varz = mean_var.parts()
  fig = pyplot.figure()
  ax = fig.add_subplot(111)
  linex, = ax.plot(n_pixels, rmsdx, label='rmsd_x')
  liney, = ax.plot(n_pixels, rmsdy, label='rmsd_y')
  ax.plot(n_pixels, flex.sqrt(varx), label='sd_x', color=linex.get_color(), linestyle='-.')
  ax.plot(n_pixels, flex.sqrt(varx), label='sd_y', color=liney.get_color(), linestyle=':')
  ax2 = ax.twinx()
  ax2.plot(n_pixels, n_ref, label='# reflections', color='grey', linestyle='--')
  ax.set_xlabel('# pixels')
  ax.set_ylabel('rmsd (mm)')
  ax2.set_ylabel('# reflections')
  ax.legend(loc='best')

  pyplot.show()


if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
