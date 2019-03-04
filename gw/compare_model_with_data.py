from __future__ import division, print_function


def compare_data_with_model(cif_file, mtz_file):
    from iotbx import cif, mtz
    from scitbx.array_family import flex
    import math
    import random

    # read model, compute Fc, square to F^2
    model = cif.reader(file_path=cif_file).build_crystal_structures()["1"]
    ic = (
        model.structure_factors(anomalous_flag=True, d_min=0.55, algorithm="direct")
        .f_calc()
        .as_intensity_array()
    )

    # read experimental measurements
    m = mtz.object(mtz_file)
    mad = m.as_miller_arrays_dict(merge_equivalents=False)
    idata = mad[("HKL_base", "HKL_base", "I")].as_anomalous_array()
    match = idata.match_indices(ic)

    # pair up, extract to vanilla arrays for easier handling
    pairs = match.pairs()

    icalc = flex.double()
    iobs = flex.double()
    sobs = flex.double()

    for p in pairs:
        iobs.append(idata.data()[p[0]])
        sobs.append(idata.sigmas()[p[0]])
        icalc.append(ic.data()[p[1]])

    # estimate conversion scale - apply F^2
    icalc *= flex.sum(iobs) / flex.sum(icalc)

    d = (iobs - icalc) / sobs

    dh = flex.histogram(d, data_min=-6, data_max=6, n_slots=120)

    m = flex.sum(d) / d.size()
    s = math.sqrt(flex.sum(d * d) / d.size() - m * m)

    # mean and standard deviation -
    print(m, s)


if __name__ == "__main__":
    import sys

    compare_data_with_model(sys.argv[1], sys.argv[2])
