from __future__ import division

import iotbx.phil

help_message = '''
'''

phil_scope = iotbx.phil.parse("""\
grid = None
  .type = ints(size=2, value_min=1)
stereographic_projections = False
  .type = bool
""", process_includes=True)


def run(args):
  from dials.util.options import OptionParser
  from scitbx.array_family import flex
  import libtbx.load_env

  usage = "%s [options] find_spots.json" %(
    libtbx.env.dispatcher_name)

  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    epilog=help_message)

  params, options, args = parser.parse_args(
    show_diff_phil=True, return_unhandled=True)

  assert len(args) == 1
  json_file = args[0]
  import json

  with open(json_file, 'rb') as f:
    results = json.load(f)

  n_indexed = flex.double()
  fraction_indexed = flex.double()
  n_spots = flex.double()
  n_lattices = flex.double()
  crystals = []
  image_names = flex.std_string()

  for r in results:
    n_spots.append(r['n_spots_total'])
    image_names.append(str(r['image']))
    if 'n_indexed' in r:
      n_indexed.append(r['n_indexed'])
      fraction_indexed.append(r['fraction_indexed'])
      n_lattices.append(len(r['lattices']))
      for d in r['lattices']:
        from dxtbx.serialize.crystal import from_dict
        crystals.append(from_dict(d['crystal']))
    else:
      n_indexed.append(0)
      fraction_indexed.append(0)
      n_lattices.append(0)

  import matplotlib
  matplotlib.use('Agg')
  from matplotlib import pyplot

  blue = '#3498db'
  red = '#e74c3c'

  marker = 'o'
  alpha = 0.5
  lw = 0

  plot = True
  table = True
  grid = params.grid

  from libtbx import group_args
  from dials.algorithms.peak_finding.per_image_analysis \
       import plot_stats, print_table

  estimated_d_min = flex.double()
  d_min_distl_method_1 = flex.double()
  d_min_distl_method_2 = flex.double()
  n_spots_total = flex.int()
  n_spots_no_ice = flex.int()
  total_intensity = flex.double()

  for d in results:
    estimated_d_min.append(d['estimated_d_min'])
    d_min_distl_method_1.append(d['d_min_distl_method_1'])
    d_min_distl_method_2.append(d['d_min_distl_method_2'])
    n_spots_total.append(d['n_spots_total'])
    n_spots_no_ice.append(d['n_spots_no_ice'])
    total_intensity.append(d['total_intensity'])

  stats = group_args(n_spots_total=n_spots_total,
                     n_spots_no_ice=n_spots_no_ice,
                     n_spots_4A=None,
                     total_intensity=total_intensity,
                     estimated_d_min=estimated_d_min,
                     d_min_distl_method_1=d_min_distl_method_1,
                     d_min_distl_method_2=d_min_distl_method_2,
                     noisiness_method_1=None,
                     noisiness_method_2=None)

  if plot:
    plot_stats(stats)
  if table:
    print_table(stats)

  if params.stereographic_projections and len(crystals):
    from dxtbx.datablock import DataBlockFactory
    datablocks = DataBlockFactory.from_filenames(
      [image_names[0]], verbose=False)
    assert len(datablocks) == 1
    imageset = datablocks[0].extract_imagesets()[0]
    s0 = imageset.get_beam().get_s0()
    # XXX what if no goniometer?
    rotation_axis = imageset.get_goniometer().get_rotation_axis()

    indices = ((1,0,0), (0,1,0), (0,0,1))
    for i, index in enumerate(indices):

      from cctbx import crystal, miller
      from cctbx.array_family import flex
      from scitbx import matrix
      miller_indices = flex.miller_index([index])
      symmetry = crystal.symmetry(
        unit_cell=crystals[0].get_unit_cell(),
        space_group=crystals[0].get_space_group())
      miller_set = miller.set(symmetry, miller_indices)
      d_spacings = miller_set.d_spacings()
      d_spacings = d_spacings.as_non_anomalous_array().expand_to_p1()
      d_spacings = d_spacings.generate_bijvoet_mates()
      miller_indices = d_spacings.indices()

      # plane normal
      d0 = matrix.col(s0).normalize()
      d1 = d0.cross(matrix.col(rotation_axis)).normalize()
      d2 = d1.cross(d0).normalize()
      reference_poles = (d0, d1, d2)

      from dials.command_line.stereographic_projection import stereographic_projection
      projections = []

      for cryst in crystals:
        reciprocal_space_points = list(cryst.get_U() * cryst.get_B()) * miller_indices.as_vec3_double()
        projections.append(stereographic_projection(
          reciprocal_space_points, reference_poles))

        #from dials.algorithms.indexing.compare_orientation_matrices import \
        #  difference_rotation_matrix_and_euler_angles
        #R_ij, euler_angles, cb_op = difference_rotation_matrix_and_euler_angles(
        #  crystals[0], cryst)
        #print max(euler_angles)

      from dials.command_line.stereographic_projection import plot_projections
      plot_projections(projections, filename='projections_%s.png' %('hkl'[i]))

  def plot_grid(values, grid, file_name, cmap=pyplot.cm.Reds,
                vmin=None, vmax=None):
    values = values.as_double()
    # At DLS, fast direction appears to be largest direction
    if grid[0] > grid[1]:
      values.reshape(flex.grid(reversed(grid)))
      values = values.matrix_transpose()
    else:
      values.reshape(flex.grid(grid))

    fig = pyplot.figure()
    pyplot.pcolormesh(
      values.as_numpy_array(), cmap=cmap, vmin=vmin, vmax=vmax)
    pyplot.xlim(0, grid[1])
    pyplot.ylim(0, grid[0])
    pyplot.gca().invert_yaxis()
    pyplot.axes().set_aspect('equal')
    pyplot.colorbar()
    pyplot.savefig(file_name)
    pyplot.clf()

  if grid is not None:
    grid = tuple(reversed(grid))
    plot_grid(n_spots_total, grid, 'grid_spot_count_total.png')
    plot_grid(n_spots_no_ice, grid, 'grid_spot_count_no_ice.png')
    plot_grid(total_intensity, grid, 'grid_total_intensity.png')
    if flex.max(n_indexed) > 0:
      plot_grid(n_indexed, grid, 'grid_n_indexed.png')
      plot_grid(fraction_indexed, grid, 'grid_fraction_indexed.png')

    for i, d_min in enumerate((estimated_d_min, d_min_distl_method_1, d_min_distl_method_2)):
      from cctbx import uctbx
      d_star_sq = uctbx.d_as_d_star_sq(d_min)
      d_star_sq.set_selected(d_star_sq == 1, 0)
      vmin = flex.min(d_star_sq.select(d_star_sq > 0))
      vmax = flex.max(d_star_sq)

      vmin = flex.min(d_min.select(d_min > 0))
      vmax = flex.max(d_min)
      cmap = pyplot.cm.Reds_r
      d_min.set_selected(d_min <= 0, vmax)

      if i == 0:
        plot_grid(d_min, grid, 'grid_d_min.png', cmap=cmap, vmin=vmin, vmax=vmax)
      else:
        plot_grid(
          d_min, grid, 'grid_d_min_method_%i.png' %i, cmap=cmap, vmin=vmin, vmax=vmax)

  if flex.max(n_indexed) > 0:
    pyplot.scatter(n_spots, n_indexed, marker=marker, alpha=alpha, c=blue, lw=lw)
    xlim = pyplot.xlim()
    ylim = pyplot.ylim()
    pyplot.plot([0, max(n_spots)], [0, max(n_spots)], c=red)
    pyplot.xlim(0, xlim[1])
    pyplot.ylim(0, ylim[1])
    pyplot.xlabel('# spots')
    pyplot.ylabel('# indexed')
    pyplot.savefig('n_spots_vs_n_indexed.png')
    pyplot.clf()

    pyplot.scatter(
    n_spots, fraction_indexed, marker=marker, alpha=alpha, c=blue, lw=lw)
    pyplot.xlim(0, pyplot.xlim()[1])
    pyplot.ylim(0, pyplot.ylim()[1])
    pyplot.xlabel('# spots')
    pyplot.ylabel('Fraction indexed')
    pyplot.savefig('n_spots_vs_fraction_indexed.png')
    pyplot.clf()

    pyplot.scatter(
      n_indexed, fraction_indexed, marker=marker, alpha=alpha, c=blue, lw=lw)
    pyplot.xlim(0, pyplot.xlim()[1])
    pyplot.ylim(0, pyplot.ylim()[1])
    pyplot.xlabel('# indexed')
    pyplot.ylabel('Fraction indexed')
    pyplot.savefig('n_indexed_vs_fraction_indexed.png')
    pyplot.clf()

    pyplot.scatter(
      n_spots, n_lattices, marker=marker, alpha=alpha, c=blue, lw=lw)
    pyplot.xlim(0, pyplot.xlim()[1])
    pyplot.ylim(0, pyplot.ylim()[1])
    pyplot.xlabel('# spots')
    pyplot.ylabel('# lattices')
    pyplot.savefig('n_spots_vs_n_lattices.png')
    pyplot.clf()

  pyplot.scatter(
    estimated_d_min, d_min_distl_method_1, marker=marker, alpha=alpha, c=blue, lw=lw)
  xlim = pyplot.xlim()
  ylim = pyplot.ylim()
  m = max(max(estimated_d_min), max(d_min_distl_method_1))
  pyplot.plot([0, m], [0, m], c=red)
  pyplot.xlim(0, xlim[1])
  pyplot.ylim(0, ylim[1])
  pyplot.xlabel('estimated_d_min')
  pyplot.ylabel('d_min_distl_method_1')
  pyplot.savefig('d_min_vs_distl_method_1.png')
  pyplot.clf()

  pyplot.scatter(
    estimated_d_min, d_min_distl_method_2, marker=marker, alpha=alpha, c=blue, lw=lw)
  xlim = pyplot.xlim()
  ylim = pyplot.ylim()
  m = max(max(estimated_d_min), max(d_min_distl_method_2))
  pyplot.plot([0, m], [0, m], c=red)
  pyplot.xlim(0, xlim[1])
  pyplot.ylim(0, ylim[1])
  pyplot.xlabel('estimated_d_min')
  pyplot.ylabel('d_min_distl_method_2')
  pyplot.savefig('d_min_vs_distl_method_2.png')
  pyplot.clf()

  pyplot.scatter(
    d_min_distl_method_1, d_min_distl_method_2, marker=marker, alpha=alpha, c=blue, lw=lw)
  xlim = pyplot.xlim()
  ylim = pyplot.ylim()
  m = max(max(d_min_distl_method_1), max(d_min_distl_method_2))
  pyplot.plot([0, m], [0, m], c=red)
  pyplot.xlim(0, xlim[1])
  pyplot.ylim(0, ylim[1])
  pyplot.xlabel('d_min_distl_method_1')
  pyplot.ylabel('d_min_distl_method_2')
  pyplot.savefig('distl_method_1_vs_distl_method_2.png')
  pyplot.clf()

  pyplot.scatter(
    n_spots, estimated_d_min, marker=marker, alpha=alpha, c=blue, lw=lw)
  pyplot.xlim(0, pyplot.xlim()[1])
  pyplot.ylim(0, pyplot.ylim()[1])
  pyplot.xlabel('# spots')
  pyplot.ylabel('estimated_d_min')
  pyplot.savefig('n_spots_vs_d_min.png')
  pyplot.clf()

  pyplot.scatter(
    n_spots, d_min_distl_method_1, marker=marker, alpha=alpha, c=blue, lw=lw)
  pyplot.xlim(0, pyplot.xlim()[1])
  pyplot.ylim(0, pyplot.ylim()[1])
  pyplot.xlabel('# spots')
  pyplot.ylabel('d_min_distl_method_1')
  pyplot.savefig('n_spots_vs_distl_method_1.png')
  pyplot.clf()

  pyplot.scatter(
    n_spots, d_min_distl_method_2, marker=marker, alpha=alpha, c=blue, lw=lw)
  pyplot.xlim(0, pyplot.xlim()[1])
  pyplot.ylim(0, pyplot.ylim()[1])
  pyplot.xlabel('# spots')
  pyplot.ylabel('d_min_distl_method_2')
  pyplot.savefig('n_spots_vs_distl_method_2.png')
  pyplot.clf()


if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
