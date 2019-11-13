from __future__ import division, print_function

import random

from dials.array_family import flex

from dials.util.options import OptionParser
from dials.util.options import flatten_experiments, flatten_reflections


def main():
    parser = OptionParser(read_experiments=True, read_reflections=True)

    params, options = parser.parse_args(show_diff_phil=True)

    experiments = flatten_experiments(params.input.experiments)
    reflections = flatten_reflections(params.input.reflections)

    data = reflections[0]
    expt = experiments[0]
    data = data.select(data.get_flags(data.flags.integrated))

    # data = data.select(data['d'] > 2.0)

    z = data["xyzcal.px"].parts()[2]

    shoebox = data["shoebox"]

    # make a new reflection table, this one empty - then add columns for
    # every event in this set (with some shuffle) - flags = 32 apparently,
    # id == id above, intensity.sum.value=1 variance=1 n_signal=1 panel=panel
    # xyzobs.px.value = pixel + random.random() - 0.5 variance = 1/12 in each
    # direction - then map to reciprocal space

    events = flex.vec3_double()

    for s in shoebox:
        d = s.data
        k0, j0, i0 = s.bbox[0], s.bbox[2], s.bbox[4]
        k1, j1, i1 = d.focus()
        for k in range(k1):
            for j in range(j1):
                for i in range(i1):
                    for n in range(int(d[k, j, i])):
                        if random.random() > 0.1:
                            continue
                        z = k + k0 + random.random()
                        y = j + j0 + random.random()
                        x = i + i0 + random.random()
                        events.append((z, y, x))

    rt = flex.reflection_table()
    rt["xyzobs.px.value"] = events
    variance = flex.double(events.size(), 1.0 / 12.0)
    rt["xyzobs.px.variance"] = flex.vec3_double(variance, variance, variance)
    rt["flags"] = flex.size_t(events.size(), 32)
    rt["id"] = flex.int(events.size(), 0)
    rt["panel"] = flex.size_t(events.size(), 0)
    rt["intensity.sum.value"] = flex.double(events.size(), 1)
    rt["intensity.sum.variance"] = flex.double(events.size(), 1)
    rt.as_file("events.refl")


main()
