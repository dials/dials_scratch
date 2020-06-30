from dials.array_family import flex
from dxtbx.model.experiment_list import ExperimentList


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
    """Combine experiment files from find to a single output file fout. 
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
