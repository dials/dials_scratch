#!/usr/bin/env python

# LIBTBX_SET_DISPATCHER_NAME dev.dials.frame_orientations

"""
Print a table of the orientation for every image of a dataset. The
orientation is expressed as a zone axis (a direction referenced to the direct
lattice) [uvw] giving the beam direction with respect to the crystal lattice.
Take into account any scan-varying models.

Usage: dev.dials.frame_orientations.py refined_experiments.json refined.pickle
"""

from __future__ import division, print_function, absolute_import
import sys
from libtbx.utils import Sorry
from dials.util.options import (flatten_reflections, flatten_experiments,
    OptionParser)
from libtbx.table_utils import simple_table
from scitbx import matrix

class Script(object):
  '''A class for running the script.'''

  def __init__(self):
    '''Initialise the script.'''
    from libtbx.phil import parse
    import libtbx.load_env

    # The phil scope
    phil_scope = parse('''
      scale = unit *ewald_sphere_radius
        .type = choice
        .help = "Choose the scale for the direction vector in orthogonal"
                "coordinates prior to transformation into fractional"
                "coordinates [uvw]"
    ''', process_includes=True)

    # The script usage
    import __main__
    usage = ("usage: dials.python {0} refined_experiments.json "
             "refined.pickle").format(__main__.__file__)

    # Create the parser
    self.parser = OptionParser(
      usage=usage,
      phil=phil_scope,
      read_experiments=True,
      epilog=__doc__)

    return

  def run(self):
    '''Execute the script.'''

    # Parse the command line
    self.params, options = self.parser.parse_args(show_diff_phil=True)

    if not self.params.input.experiments:
      self.parser.print_help()
      sys.exit()

    # Try to load the models
    experiments = flatten_experiments(self.params.input.experiments)
    nexp = len(experiments)
    if nexp == 0:
      print("No Experiments found in the input")
      self.parser.print_help()
      return

    header = ["Exp\nid", "Image", "Beam direction (xyz)", "Zone axis [uvw]"]
    rows = []

    for iexp, exp in enumerate(experiments):
      dat = extract_experiment_data(exp, self.params.scale)
      images = dat['images']
      directions = dat['directions']
      zone_axes = dat['zone_axes']
      for i, d, z in zip(images, directions, zone_axes):
        row = [str(iexp),
               str(i),
               "{:.3f} {:.3f} {:.3f}".format(*d.elems),
               "{:.3f} {:.3f} {:.3f}".format(*z.elems)]
        rows.append(row)

    st = simple_table(rows, header)
    print(st.format())

    return

def extract_experiment_data(exp, scale=''):
  """Extract lists of the image number, beam direction and zone axis from an
  experiment"""
  crystal = exp.crystal
  beam = exp.beam
  scan = exp.scan
  gonio = exp.goniometer

  if scale == 'ewald_sphere_radius':
    scale = 1.0 / beam.get_wavelength()
  else:
    scale = 1.0

  image_range = scan.get_image_range()
  images = range(image_range[0], image_range[1] + 1)

  if beam.num_scan_points > 0:
    # There is one more scan point than the number of images. For simplicity,
    # omit the final scan point, to leave a list of the beam directions at
    # the _beginning_ of each image. Likewise for gonio and crystal, below.
    directions = []
    for i in range(beam.num_scan_points - 1):
      s0 = matrix.col(beam.get_s0_at_scan_point(i))
      directions.append(s0.normalize() * scale)
  else:
    directions = [matrix.col(beam.get_unit_s0()) * scale for i in images]

  if gonio.num_scan_points > 0:
    S_mats = []
    for i in range(gonio.num_scan_points - 1):
      S = matrix.sqr(gonio.get_setting_rotation_at_scan_point(i))
      S_mats.append(S)
  else:
    S_mats = [matrix.sqr(gonio.get_setting_rotation()) for i in images]

  F_mats = [matrix.sqr(gonio.get_fixed_rotation()) for i in images]
  array_range = scan.get_array_range()
  R_mats = []
  axis = matrix.col(gonio.get_rotation_axis_datum())
  for i in range(*array_range):
    phi = scan.get_angle_from_array_index(i, deg=True)
    R = matrix.sqr(axis.axis_and_angle_as_r3_rotation_matrix(phi, deg=True))
    R_mats.append(R)

  if crystal.num_scan_points > 0:
    UB_mats = []
    for i in range(crystal.num_scan_points - 1):
      UB = matrix.sqr(crystal.get_A_at_scan_point(i))
      UB_mats.append(UB)
  else:
    UB_mats = [matrix.sqr(crystal.get_A()) for i in images]

  assert len(directions) == len(S_mats) == len(F_mats) == len(R_mats) == len(UB_mats)

  setting_rotation = matrix.sqr(gonio.get_setting_rotation())
  rotation_axis = matrix.col(gonio.get_rotation_axis_datum())

  # Construct full setting matrix for each image
  SRFUB = []
  for S, R, F, UB in zip(S_mats, R_mats, F_mats, UB_mats):
    SRFUB.append(S * R * F * UB)

  # SFRUB is the orthogonalisation matrix for the reciprocal space laboratory
  # frame. We want the real space fractionalisation matrix, which is its
  # transpose (https://dials.github.io/documentation/conventions.html)
  frac_mats = [m.transpose() for m in SRFUB]

  zone_axes = [F * d for F, d in zip(frac_mats, directions)]

  return {'images':images,
          'directions':directions,
          'zone_axes':zone_axes}

if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
