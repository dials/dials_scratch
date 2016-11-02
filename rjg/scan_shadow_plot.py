# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

from __future__ import division
from scitbx.array_family import flex
import libtbx.phil
from dials_scratch.rjg.multi_axis_shadow_map import polygon_area

help_message = '''

'''

phil_scope= libtbx.phil.parse('''
oscillation_range = None
  .type = floats(size=2)
step_size = 1
  .type = float(value_min=0)
y_max = None
  .type = float(value_min=0, value_max=100)
output {
  plot = scan_shadow_plot.png
    .type = path
  json = None
    .type = path
}
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
  scan = imageset.get_scan()
  masker = imageset.reader().get_format().get_goniometer_shadow_masker()
  angles = goniometer.get_angles()
  names = goniometer.get_names()
  scan_axis = goniometer.get_scan_axis()
  phi = angles[0]

  if params.step_size is not None:
    step = params.step_size
  else:
    step = scan.get_oscillation()[1]
  if params.oscillation_range is not None:
    start, end = params.oscillation_range
  else:
    start, end = scan.get_oscillation_range()
  scan_points = flex.double(libtbx.utils.frange(start, end, step=step))

  n_px_shadowed = flex.double(scan_points.size(), 0)
  n_px_tot = flex.double(scan_points.size(), 0)

  assert len(angles) == 3
  for i, scan_angle in enumerate(scan_points):
    shadow = masker.project_extrema(detector, scan_angle=scan_angle)
    for p_id in range(len(detector)):
      px_x, px_y = detector[p_id].get_image_size()
      n_px_tot[i] += px_x * px_y
      if shadow[p_id].size() < 4:
        continue
      n_px_shadowed[i] += polygon_area(shadow[p_id])

  fraction_shadowed = n_px_shadowed/n_px_tot

  if params.output.json is not None:
    print 'Writing json output to %s' %params.output.json
    d = {
      'scan_points': list(scan_points),
      'fraction_shadowed': list(fraction_shadowed),
    }
    import json
    with open(params.output.json, 'wb') as f:
      json.dump(d, f)

  if params.output.plot is not None:
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    plt.style.use('ggplot')
    plt.plot(scan_points.as_numpy_array(), fraction_shadowed.as_numpy_array() * 100)
    plt.xlabel('%s angle (degrees)' %names[scan_axis])
    plt.ylabel('Shadowed area (%)')
    if params.y_max is not None:
      plt.ylim(0, params.y_max)
    else:
      plt.ylim(0, plt.ylim()[1])
    plt.axes().set_xticklabels(["%.0f" %(step * t) for t in plt.xticks()[0]])
    plt.tight_layout()
    plt.savefig(params.output.plot)
    fig = plt.gcf()
    size = fig.get_size_inches()
    fig.set_size_inches((size[0]/2, size[1]/2))
    plt.tight_layout()
    plt.savefig('scan_shadow_plot_small.png')




if __name__ == '__main__':
  import sys
  run(sys.argv[1:])

