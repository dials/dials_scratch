from __future__ import division

import iotbx.phil
from dials.util.options import OptionParser
from dials.util.options import flatten_reflections
from dials.array_family import flex
import math

help_message = '''
'''

phil_scope = iotbx.phil.parse('''
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

  print reflections.size()


  x = flex.size_t()
  rmsds = flex.vec3_double()

  for npix in range(1, flex.max(reflection_size)):
    refl = reflections.select(reflection_size == npix)
    if refl.size() == 0: continue
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

    x.append(npix)
    rmsds.append((rmsdx, rmsdy, rmsdz))

  from matplotlib import pyplot
  y1, y2, y3 = rmsds.parts()
  pyplot.plot(x, y1)
  pyplot.plot(x, y2)
  #pyplot.plot(x, y3)
  pyplot.show()


if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
