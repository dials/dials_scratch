from __future__ import division

import iotbx.phil

help_message = '''
'''

phil_scope = iotbx.phil.parse("""\
grid = None
  .type = ints(size=2, value_min=1)
""", process_includes=True)


def run(args):
  from dials.util.options import OptionParser
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

  n_indexed = []
  fraction_indexed = []
  n_spots = []
  n_lattices = []

  for r in results:
    n_spots.append(r['n_spots_total'])
    if 'n_indexed' in r:
      n_indexed.append(r['n_indexed'])
      fraction_indexed.append(r['fraction_indexed'])
      n_lattices.append(len(r['lattices']))
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

  from scitbx.array_family import flex
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

  def plot_grid(values, grid, file_name):
    values = values.deep_copy()
    values.reshape(flex.grid(grid))

    fig = pyplot.figure()
    pyplot.pcolormesh(values.as_numpy_array(), cmap=pyplot.cm.Reds)
    pyplot.axes().set_aspect('equal')
    pyplot.savefig(file_name)
    pyplot.clf()


  if grid is not None:
    grid = tuple(reversed(grid))
    plot_grid(n_spots_no_ice, grid, 'grid_spot_count_no_ice.png')
    plot_grid(total_intensity, grid, 'grid_total_intensity.png')

    for i, d_min in enumerate((estimated_d_min, d_min_distl_method_1, d_min_distl_method_2)):
      from cctbx import uctbx
      d_star_sq = uctbx.d_as_d_star_sq(d_min)
      d_star_sq.set_selected(d_star_sq == 1, 0)
      if i == 0:
        plot_grid(d_star_sq, grid, 'grid_d_min.png')
      else:
        plot_grid(d_star_sq, grid, 'grid_d_min_method_%i.png' %i)

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
  pyplot.plot([0, max(estimated_d_min)], [0, max(d_min_distl_method_1)], c=red)
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
  pyplot.plot([0, max(estimated_d_min)], [0, max(d_min_distl_method_2)], c=red)
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
  pyplot.plot([0, max(d_min_distl_method_1)], [0, max(d_min_distl_method_2)], c=red)
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
