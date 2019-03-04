# -*- coding: utf-8 -*-
#!/usr/bin/env python
#
# delta_psi_weights.py
#
#  Copyright (C) 2017 Lawrence Berkeley National Laboratory (LBNL)
#
#  Author: Aaron Brewster
#
#  This code is distributed under the X license, a copy of which is
#  included in the root directory of this package.
#
from __future__ import division
from __future__ import print_function
from dials.array_family import flex
from matplotlib import pyplot as plt
from libtbx.phil import parse
from cctbx.crystal import symmetry as crystal_symmetry
from cctbx.miller import set as miller_set
from cctbx import uctbx
from libtbx import easy_pickle

help_message = """

Description
  Weight delta psi by various methods for use in refinement. Weights are added to the reflection
  table using the key 'delpsical.weights'.

  summed_intensity: first, a wilson analysis is performed, where the binning is computed from the
  average unit cell of the input experiments. The delta psi weight for a reflection is the ratio of
  intensity.sum.value of that reflection to the mean of intensity.sum.value in that reflection's
  resolution bin.

Example:

  libtbx.python delta_psi_weights.py experiment.json reflections.pickle
"""

phil_scope = parse(
    """
  method = *summed_intensity
    .type = choice
    .help = Weight delta psi by the ratio of the observed intensity to the average intensity in the \
            reflection's resolution bin
  summed_intensity {
    scale_factor = 1.0
      .type = float
      .help = Multiplier applied to all weights after they are computed
  }
  show_weight_plots = False
    .type = bool
    .help = Show histograms of the weights per resolution bin
  output {
    reflections = weighted_reflections.pickle
      .type = str
      .help = output file name for the reflection table
  }
"""
)


class Script(object):
    """ Class to parse the command line options. """

    def __init__(self):
        """ Set the expected options. """
        from dials.util.options import OptionParser
        import libtbx.load_env

        # Create the option parser
        usage = (
            "usage: %s [options] /path/to/refined/json/file"
            % libtbx.env.dispatcher_name
        )
        self.parser = OptionParser(
            usage=usage,
            phil=phil_scope,
            read_experiments=True,
            read_reflections=True,
            epilog=help_message,
        )

    def run(self):
        """ Parse the options. """
        # Parse the command line arguments
        params, options = self.parser.parse_args(show_diff_phil=True)
        self.params = params
        assert (
            len(params.input.experiments) == len(params.input.reflections) == 1
        ), "Provide one experiment list and one reflection table"
        assert params.method == "summed_intensity"

        experiments = params.input.experiments[0].data
        reflections = params.input.reflections[0].data

        # Find the aveage unit cell for the crystals in the experiments provided
        weighted_relfs = flex.reflection_table()
        all_uc = [flex.double() for i in xrange(6)]
        space_group = None
        for expt_id, experiment in enumerate(experiments):
            refls = reflections.select(reflections["id"] == expt_id)
            unit_cell = experiment.crystal.get_unit_cell()
            for i in xrange(6):
                all_uc[i].append(unit_cell.parameters()[i])

            if space_group is None:
                space_group = experiment.crystal.get_space_group()
            else:
                assert (
                    space_group.type().lookup_symbol()
                    == experiment.crystal.get_space_group().type().lookup_symbol()
                )

        # Compute the average unit cell and build a miller array with it
        unit_cell = uctbx.unit_cell([flex.mean(data) for data in all_uc])
        cs = crystal_symmetry(unit_cell, space_group.type().lookup_symbol())
        ms = miller_set(cs, reflections["miller_index"], anomalous_flag=False)
        ma = ms.array(
            reflections["intensity.sum.value"]
            / flex.sqrt(reflections["intensity.sum.variance"])
        )

        ma.setup_binner(n_bins=10)
        binner = ma.binner()
        mean_i = flex.double()
        reflections["delpsical.weights"] = flex.double(len(reflections), 0)

        # Iterate through the bins and compute the Wilson plot, then use it compute the weights
        for i in binner.range_all():
            sel = binner.selection(i)
            if sel.count(True) == 0:
                mean_i.append(0)
                continue
            mean_i.append(flex.mean(reflections["intensity.sum.value"].select(sel)))
            reflections["delpsical.weights"].set_selected(
                sel,
                reflections["intensity.sum.value"].select(sel)
                * (params.summed_intensity.scale_factor / mean_i[i]),
            )

            if params.show_weight_plots:
                fig = plt.figure()
                plt.title(str(i))
                plt.hist(reflections["delpsical.weights"].select(sel))

        # Show unit cell distribution and mean I
        print("Average uc +/- std. deviation")
        labels = ["% 6s" % l for l in ["a", "b", "c", "alpha", "beta", "gamma"]]
        for label, data in zip(labels, all_uc):
            stats = flex.mean_and_variance(data)
            print(
                "%s % 6.1f +/- %6.1f"
                % (label, stats.mean(), stats.unweighted_sample_standard_deviation())
            )

        print("Mean I over all data")
        binner.show_data(mean_i, data_fmt="%.1f", show_unused=False)

        easy_pickle.dump(params.output.reflections, reflections)

        if params.show_weight_plots:
            plt.show()


if __name__ == "__main__":
    from dials.util import halraiser

    try:
        script = Script()
        script.run()
    except Exception as e:
        halraiser(e)
