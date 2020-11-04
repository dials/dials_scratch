from dials.array_family import flex
from dxtbx.model.experiment_list import ExperimentList


def combine_reflections(fins, fout):
    """Combine reflection files from fins into a single output file fout.
    Makes assumptions that the data are consistent from one to the next, no
    particular tests are performed at the moment."""

    d0 = flex.reflection_table.from_file(fins[0])
    if "intensity.prf.variance" in d0:
        d0 = d0.select(d0["intensity.prf.variance"] > 0)
    d0 = d0.select(d0["intensity.sum.variance"] > 0)

    for f in fins[1:]:
        d1 = flex.reflection_table.from_file(f)
        if "intensity.prf.variance" in d1:
            d1 = d1.select(d1["intensity.prf.variance"] > 0)
        d1 = d1.select(d1["intensity.sum.variance"] > 0)

        d1.experiment_identifiers()[0] = d0.experiment_identifiers()[0]

        matches = d0.match(d1)

        for i, j in zip(*matches):
            assert d0["miller_index"][i] == d1["miller_index"][j]

        # join up those reflections, extract them from the existing lists...
        # copy the unmatched reflections over to copy of the output data
        # then merge the partials and copy them in

        # detail: select(matches[]) reorders the reflections as well

        m0 = flex.bool(d0.size(), False)
        m0.set_selected(matches[0], True)
        d0p = d0.select(matches[0])
        d0 = d0.select(~m0)

        m1 = flex.bool(d1.size(), False)
        m1.set_selected(matches[1], True)
        d1p = d1.select(matches[1])
        d1 = d1.select(~m1)

        d0.extend(d1)

        prf = "intensity.prf.value" in d0p

        # combine d0p, d1p - copy the information from d1p[j] to d0p[j]
        for j in range(d0p.size()):
            d0p["partiality"][j] += d1p["partiality"][j]

            if prf:

                # weight profiles by (I/sig(I))^2 as done in dials.export -
                # should probably prune the reflection lists in here to only
                # include useful measurements before I get to this point...

                i0 = d0p["intensity.prf.value"][j]
                v0 = d0p["intensity.prf.variance"][j]
                w0 = i0 * i0 / v0

                i1 = d1p["intensity.prf.value"][j]
                v1 = d1p["intensity.prf.variance"][j]
                w1 = i1 * i1 / v1

                if w0 + w1 > 0:
                    i = (i0 * w0 + i1 * w1) / (w0 + w1)
                    v = (v0 * w0 + v1 * w1) / (w0 + w1)

                    d0p["intensity.prf.value"][j] = i
                    d0p["intensity.prf.variance"][j] = v

            d0p["intensity.sum.value"][j] += d1p["intensity.sum.value"][j]
            d0p["intensity.sum.variance"][j] += d1p["intensity.sum.variance"][j]

            d0p["background.sum.value"][j] += d1p["background.sum.value"][j]
            d0p["background.sum.variance"][j] += d1p["background.sum.variance"][j]

            # extend the bbox... sort out the flags...
            b = d0p["bbox"][j]
            _b = d1p["bbox"][j]
            d0p["bbox"][j] = (b[0], _b[1], b[2], b[3], b[4], b[5])

        d0.extend(d0p)

    d0.as_file(fout)


def combine_experiments(fins, fout):
    """Combine experiment files from fins to a single output file fout.
    Makes assumptions that the crystal model is consistent from one experiment
    to the next, that the detectors are all similar."""

    e0 = ExperimentList.from_file(fins[0])
    c0 = e0[0].crystal
    s0 = e0[0].scan
    mats = [c0.get_A_at_scan_point(i) for i in range(c0.num_scan_points)]
    for f in fins[1:]:
        e1 = ExperimentList.from_file(f)
        c1 = e1[0].crystal
        s1 = e1[0].scan
        s0 += s1
        mats.extend(
            [c1.get_A_at_scan_point(i + 1) for i in range(c1.num_scan_points - 1)]
        )
    c0.set_A_at_scan_points(mats)
    e0.as_file(fout)
