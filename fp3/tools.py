from dials.array_family import flex
from dxtbx.model.experiment_list import ExperimentList


def combine_partial_reflections(fins, fout):
    """Combine reflection files from fins into a single output file fout.
    Makes assumptions that the data are consistent from one to the next, no
    particular tests are performed at the moment."""

    d0 = flex.reflection_table.from_file(fins[0])
    for f in fins[1:]:
        d1 = flex.reflection_table.from_file(f)

        d1.experiment_identifiers()[0] = d0.experiment_identifiers()[0]

        matches = d0.match(d1)

        for i, j in zip(*matches):
            assert d0["miller_index"][i] == d1["miller_index"][j]

        # join up those reflections, extract them from the existing lists...
        m0 = flex.bool(d0.size(), False)
        m0.set_selected(matches[0], True)
        d0p = d0.select(m0)
        d0 = d0.select(~m0)

        m1 = flex.bool(d1.size(), False)
        m1.set_selected(matches[1], True)
        d1p = d1.select(m1)
        d1 = d1.select(~m1)

        d0.extend(d1)

        # combine d0p, d1p
        for j in range(d0p.size()):
            d0p["partial_id"][j] = 100000000 + j
            d1p["partial_id"][j] = 100000000 + j

        d01 = flex.reflection_table()
        d01.extend(d0p)
        d01.extend(d1p)

        d01 = sum_partial_reflections(d01)
        d0.extend(d01)

    d0.as_file(fout)


def combine_reflections(fins, fout):
    """Combine reflection files from fins into a single output file fout.
    Makes assumptions that the data are consistent from one to the next, no
    particular tests are performed at the moment."""

    d0 = flex.reflection_table.from_file(fins[0])
    for f in fins[1:]:
        d1 = flex.reflection_table.from_file(f)

        d1.experiment_identifiers()[0] = d0.experiment_identifiers()[0]
        d0.extend(d1)

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
