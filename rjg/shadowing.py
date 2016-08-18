from __future__ import division

from mpl_toolkits.mplot3d import Axes3D
from cctbx.array_family import flex

import iotbx.phil

help_message = '''

Examples::

'''

phil_scope = iotbx.phil.parse(
'''
height = 50
  .type = float
  .help = "Height of cone around goniometer phi axis in mm"
radius = 18
  .type = float
  .help = "Radius cone around goniometer phi axis in mm"
angle = None
  .type = float
''')

def run(args):

  from dials.util.options import OptionParser
  from dials.util.options import flatten_datablocks
  from libtbx.utils import Sorry
  import libtbx.load_env

  usage = "%s [options] experiments.json" %libtbx.env.dispatcher_name

  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    read_datablocks=True,
    check_format=False,
    epilog=help_message)

  params, options = parser.parse_args(show_diff_phil=True)
  datablocks = flatten_datablocks(params.input.datablock)
  if len(datablocks) == 0:
    parser.print_help()
    return

  imagesets = []
  for datablock in datablocks:
    imagesets.extend(datablock.extract_imagesets())

  for imageset in imagesets:
    import math

    height = params.height # mm
    radius = params.radius # mm

    theta = flex.double([range(360)]) * math.pi/180
    y = radius * flex.cos(theta) # x
    z = radius * flex.sin(theta) # y
    x = flex.double(theta.size(), height) # z

    coords = flex.vec3_double(zip(x, y, z))
    coords.insert(0, (0,0,0))

    from scitbx import matrix
    gonio = imageset.get_goniometer()
    print gonio
    scan = imageset.get_scan()
    beam = imageset.get_beam()
    detector = imageset.get_detector()

    axes = gonio.get_axes()
    angles = gonio.get_angles()
    names = gonio.get_names()

    for i in range(len(axes)):
      axis = axes[i]
      if i == gonio.get_scan_axis() and params.angle is not None:
        angle = params.angle
      else:
        angle = angles[i]
      print names[i], axis, angle
      rotation = matrix.col(axis).axis_and_angle_as_r3_rotation_matrix(angle, deg=True)
      coords = rotation.elems * coords

    mask = [flex.bool(flex.grid(p.get_image_size()), False)
            for p in detector.iter_panels()]

    shadow_boundary = [flex.vec2_double() for i in range(len(detector))]

    norm_coords = coords.deep_copy()
    lengths = coords.norms()
    sel = lengths > 0
    norm_coords.set_selected(sel, coords.select(sel)/lengths.select(sel))
    norm_coords *= 1/beam.get_wavelength()

    for c in norm_coords:
      p_id = detector.get_panel_intersection(c)
      if p_id == -1: continue
      px = detector[p_id].get_ray_intersection_px(c)
      shadow_boundary[p_id].append(px)

    import matplotlib.pyplot as plt

    fig = plt.figure()

    x, y, z = coords.parts()
    plt.scatter(x.as_numpy_array(), y.as_numpy_array())
    plt.axes().set_aspect('equal')
    plt.xlabel('x (gonio axis)')
    plt.ylabel('y (perpendicular to beam)')
    plt.show()

    plt.scatter(y.as_numpy_array(), z.as_numpy_array())
    plt.axes().set_aspect('equal')
    plt.xlabel('y (perpendicular to beam)')
    plt.ylabel('z (towards beam))')
    plt.show()

    plt.scatter(z.as_numpy_array(), x.as_numpy_array())
    plt.axes().set_aspect('equal')
    plt.xlabel('z (towards beam)')
    plt.ylabel('x (gonio axis)')
    plt.show()

    for p_id in range(len(detector)):
      x, y = shadow_boundary[p_id].parts()
      fig = plt.figure()
      plt.scatter(x.as_numpy_array(), y.as_numpy_array(), c='r')
      plt.axes().set_aspect('equal')
      plt.xlim(0, detector[p_id].get_image_size()[0])
      plt.ylim(0, detector[p_id].get_image_size()[0])
      plt.gca().invert_yaxis()
      plt.show()

if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
