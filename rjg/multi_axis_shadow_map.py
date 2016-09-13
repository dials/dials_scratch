# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1
# LIBTBX_SET_DISPATCHER_NAME dev.dials.check_strategy

from __future__ import division
from scitbx.array_family import flex
import libtbx.phil

help_message = '''

'''

phil_scope= libtbx.phil.parse('''
step_size = 10
  .type = float(value_min=0)
''')


def run(args):

  from dials.util.options import OptionParser
  from dials.util.options import flatten_datablocks
  import libtbx.load_env

  usage = "%s [options] datablock.json" %(
    libtbx.env.dispatcher_name)

  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    read_datablocks=True,
    check_format=True,
    epilog=help_message)

  params, options = parser.parse_args(show_diff_phil=True)
  datablocks = flatten_datablocks(params.input.datablock)

  if len(datablocks) == 0:
    parser.print_help()
    exit(0)

  assert len(datablocks) == 1
  imagesets = datablocks[0].extract_imagesets()

  imageset = imagesets[0]
  goniometer = imageset.get_goniometer()
  detector = imageset.get_detector()
  masker = imageset.reader().get_format().get_goniometer_shadow_masker()
  angles = goniometer.get_angles()
  names = goniometer.get_names()
  phi = angles[0]

  step = params.step_size
  kappa_values = flex.double(libtbx.utils.frange(0, 360, step=step))
  omega_values = flex.double(libtbx.utils.frange(0, 360, step=step))

  grid = flex.grid(kappa_values.size(), omega_values.size())
  n_px_shadowed = flex.double(grid, 0)
  n_px_tot = flex.double(grid, 0)

  assert len(angles) == 3
  for i, kappa in enumerate(kappa_values):
    for j, omega in enumerate(omega_values):
      masker.goniometer.set_angles((phi, kappa, omega))
      shadow = masker.project_extrema(detector, scan_angle=omega)
      for p_id in range(len(detector)):
        px_x, px_y = detector[p_id].get_image_size()
        n_px_tot[i, j] += px_x * px_y
        if shadow[p_id].size() < 4:
          continue
        n_px_shadowed[i, j] += polygon_area(shadow[p_id])

  fraction_shadowed = n_px_shadowed/n_px_tot

  from matplotlib import pyplot as plt
  plt.imshow(fraction_shadowed.as_numpy_array() * 100, interpolation='bicubic')
  plt.xlabel('%s angle (degrees)' %names[2])
  plt.ylabel('%s angle (degrees)' %names[1])
  plt.axes().set_xticklabels(["%.0f" %(step * t) for t in plt.xticks()[0]])
  plt.axes().set_yticklabels(["%.0f" %(step * t) for t in plt.yticks()[0]])
  cbar = plt.colorbar()
  cbar.set_label('Shadowed area (%)')
  plt.savefig('shadow_map.png')


def polygon_area(points):
  #http://mathworld.wolfram.com/PolygonArea.html
  x0, y0 = points.parts()
  x1 = x0[1:]
  x1.append(x0[0])
  y1 = y0[1:]
  y1.append(y0[0])

  return 0.5 * abs(flex.sum(x0 * y1 - x1 * y0))



if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
