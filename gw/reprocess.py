import sys

from libtbx.phil import parse

from dxtbx.model.experiment_list import ExperimentList
from dials.array_family import flex

from dials.util.options import OptionParser, reflections_and_experiments_from_files
from dials.algorithms.integration.integrator import create_integrator
from dials.command_line.integrate import working_phil as original_phil
from dials.algorithms.scaling.error_model.error_model import phil_scope as error_scope

# Local overrides for dials.integrate
phil_overrides = parse(
    """
output {
  reflections = reprocessed.refl
    .type = path
    .help = "Output filename for reprocessed data"
  experiments = reprocessed.expt
    .type = path
    .help = "Output filename for reprocessed models (for reference)"
}
"""
)

working_phil = original_phil.fetch(sources=[phil_overrides])


def reprocess(params, experiments, reflections):

    # only keep the integrated reflections
    reflections = reflections.select(
        reflections.get_flags(reflections.flags.integrated_sum)
    )

    # remove all the intensity columns
    for column in [k for k in reflections.keys() if k.startswith("intensity")]:
        del reflections[column]

    # actually perform the integration
    integrator = create_integrator(params, experiments, reflections)
    reflections = integrator.integrate()

    # apply the existing scale factors and error model to these reflections
    lp = reflections["lp"]
    qe = reflections["qe"]
    part = reflections["partiality"]
    scale = lp / (qe * part)

    reflections["intensity.scale.value"] = reflections["intensity.sum.value"] * scale
    reflections["intensity.scale.variance"] = (
        reflections["intensity.sum.variance"] * scale * scale
    )

    # get and apply the original error model: TODO make this more general
    scaling_model = experiments[0].scaling_model
    scaling_model.load_error_model(error_scope.extract())
    reflections[
        "intensity.scale.variance"
    ] = flex.double(scaling_model.error_model.update_variances(
        reflections["intensity.scale.variance"], reflections["intensity.scale.value"]
    ))

    reflections.as_file(params.output.reflections)
    experiments.as_file(params.output.experiments)


if __name__ == "__main__":

    parser = OptionParser(
        phil=working_phil,
        read_experiments=True,
        read_reflections=True,
    )

    imported = ExperimentList.from_file(sys.argv[1])

    params, options = parser.parse_args(args=sys.argv[2:], show_diff_phil=True)

    reflections, experiments = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )

    experiments[0].imageset = imported[0].imageset

    reflections = reflections[0]
    reflections = reflections.select(reflections.get_flags(reflections.flags.scaled))

    reprocess(params, experiments, reflections)
