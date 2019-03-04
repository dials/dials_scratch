def read_expt(filename):
    from dials.phil import ExperimentListConverters
    from dials.util.options import flatten_experiments

    converter = ExperimentListConverters(check_format=False)
    return flatten_experiments([converter.from_string(filename)])


def write_expt(experiments, filename):
    from dxtbx.model.experiment_list import ExperimentListDumper

    dump = ExperimentListDumper(experiments)
    dump.as_json(filename)


import sys

expts = read_expt(sys.argv[1])
expt = expts[0]

scan = expt.scan

epochs = scan.get_epochs()
exposure_times = scan.get_exposure_times()
image_range = scan.get_image_range()
oscillation = scan.get_oscillation()

current = 1 + image_range[1] - image_range[0]
turn = int(round(360.0 / oscillation[1]))
extra = turn - current

for j in range(extra):
    epochs.append(0.0)
    exposure_times.append(0.0)

image_range = image_range[0], image_range[1] + extra

scan.set_image_range(image_range)
scan.set_epochs(epochs)
scan.set_exposure_times(exposure_times)

write_expt(expts, sys.argv[2])
