#!/usr/bin/env python
#
# LIBTBX_SET_DISPATCHER_NAME dev.dials.resolution_quality_estimate

from __future__ import division

help_message = '''
Script to estimate diffraction limit given a directory filled with integration results.
'''

from libtbx.phil import parse
phil_str = '''
  d_min = 2.0
    .type = float
  n_bins = 10
    .type = int
  sig_filter_sigma = 0.1
    .type = float
  best_count = 20
    .type = int
'''

phil_scope = parse(phil_str, process_includes=True)

class Script(object):
  '''A class for running the script.'''

  def __init__(self):
    '''Initialise the script.'''
    from dials.util.options import OptionParser
    import libtbx.load_env

    # The script usage
    usage = "usage: %s [options] [param.phil] directory name" % libtbx.env.dispatcher_name

    self.tag = None

    # Create the parser
    self.parser = OptionParser(
      usage=usage,
      phil=phil_scope,
      read_experiments=False,
      read_reflections=False,
      epilog=help_message
      )

  def run(self):
    '''Execute the script.'''
    import os, math
    from cctbx.crystal import symmetry
    from scitbx.array_family import flex
    from libtbx import table_utils, easy_pickle
    from xfel.command_line.cspad_cbf_metrology import find_files
    from dxtbx.model.experiment_list import ExperimentListFactory
    table_header = ["","","","I","IsigI","N >","RMSD","Cutoff"]
    table_header2 = ["Bin","Resolution Range","Completeness","","","cutoff","(um)",""]

    # Parse the command line
    params, options, all_paths = self.parser.parse_args(show_diff_phil=False, return_unhandled=True)
    exp_paths = []
    refl_paths = []
    for path in all_paths:
      exps, refs = find_files(path, "integrated")
      exp_paths.extend(exps)
      refl_paths.extend(refs)
    assert len(exp_paths) == len(refl_paths)

    best_data = {}
    best_limits = flex.double()
    for exp_path, refl_path in zip(exp_paths, refl_paths):
      experiments = ExperimentListFactory.from_json_file(exp_path, check_format=False)
      reflections = easy_pickle.load(refl_path)
      exp_name = os.path.basename(exp_path)
      if exp_name.startswith("idx-") and exp_name.endswith("_refined_experiments.json"):
        tag = exp_name.lstrip("idx-").rstrip("_refined_experiments.json")
      else:
        tag = "%s, %s"%(exp_path, refl_path)

      for exp_id, experiment in enumerate(experiments):
        print "*"*80
        print "Data table for", tag
        table_data = []
        table_data.append(table_header)
        table_data.append(table_header2)

        crystal = experiment.crystal
        refls = reflections.select(reflections['id'] == exp_id)
        sym = symmetry(unit_cell = crystal.get_unit_cell(), space_group = crystal.get_space_group())
        d = crystal.get_unit_cell().d(refls['miller_index'])
        mset = sym.miller_set(indices = refls['miller_index'].select(d>=params.d_min), anomalous_flag=False)
        binner = mset.setup_binner(n_bins=params.n_bins)
        acceptable_resolution_bins = []
        for i in binner.range_used():
          d_max, d_min = binner.bin_d_range(i)
          sel = (d <= d_max) & (d > d_min)
          sel &= refls['intensity.sum.value'] > 0
          bin_refls = refls.select(sel)
          n_refls = len(bin_refls)
          avg_i = flex.mean(bin_refls['intensity.sum.value']) if n_refls > 0 else 0
          avg_i_sigi = flex.mean(bin_refls['intensity.sum.value'] /
                                 flex.sqrt(bin_refls['intensity.sum.variance'])) if n_refls > 0 else 0
          acceptable_resolution_bins.append(avg_i_sigi >= params.sig_filter_sigma)

          bright_refls = bin_refls.select((bin_refls['intensity.sum.value']/flex.sqrt(bin_refls['intensity.sum.variance'])) >= params.sig_filter_sigma)
          n_bright = len(bright_refls)

          rmsd_obs = 1000*math.sqrt((bright_refls['xyzcal.mm']-bright_refls['xyzobs.mm.value']).sum_sq()/n_bright) if n_bright > 0 else 0

          table_row = []
          table_row.append("%3d"%i)
          table_row.append("%-13s"%binner.bin_legend(i_bin=i,show_bin_number=False,show_bin_range=False,
                                                     show_d_range=True, show_counts=False))
          table_row.append("%13s"%binner.bin_legend(i_bin=i,show_bin_number=False,show_bin_range=False,
                                                    show_d_range=False, show_counts=True))

          table_row.append("%.1f"%(avg_i))
          table_row.append("%.1f"%(avg_i_sigi))
          table_row.append("%3d"%n_bright)
          table_row.append("%.1f"%(rmsd_obs))
          table_data.append(table_row)

        acceptable_resolution_bins = [acceptable_resolution_bins[i] for i in xrange(len(acceptable_resolution_bins))
                                      if False not in acceptable_resolution_bins[:i+1]]

        for b, row in zip(acceptable_resolution_bins, table_data[2:]):
          if b:
            row.append("X")
        print table_utils.format(table_data,has_header=2,justify='center',delim=" ")
        print tag, "unit cell:", ", ".join(["%.2f"%p for p in crystal.get_unit_cell().parameters()]), crystal.get_space_group().info()

        if any(acceptable_resolution_bins):
          best_index = acceptable_resolution_bins.count(True)-1
          best_row = table_data[best_index+2]
          d_min = binner.bin_d_range(binner.range_used()[best_index])[1]
          if len(best_limits) < params.best_count:
            best_limits.append(d_min)
            best_data[tag] = d_min, best_row
          elif (d_min < best_limits).count(True) > 0:
            worst_d_min = flex.max(best_limits)
            for tag, data in best_data.iteritems():
              if worst_d_min == data[0]:
                best_data[tag] = d_min, best_row
                best_limits[flex.first_index(best_limits, worst_d_min)] = d_min
                break
          print tag, "best row:", " ".join(best_row)
        else:
          print "Data didn't pass cutoff"
    if len(best_limits) > 0:
      print "*"*80
      print "Top", len(best_limits)
      for tag, data in best_data.iteritems():
        print tag, " ".join(data[1])

if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
