from __future__ import division, print_function

from libtbx import phil

phil_scope = phil.parse(
    """
  data = None
    .type = path
    .help = "Data to test (as pickle file)"
  reference = None
    .type = path
    .help = "Reference MTZ file (containing IMEAN, SIGIMEAN)"
  kb_scale = false
    .type = bool
    .help = "Apply kB scaling to data for comparision"
  filter_integrated = true
    .type = bool
    .help = "Filter integrated flag"
  filter_partiality = 0.0
    .type = float
    .help = "Filter partiality column"
  out = None
    .type = path
    .help = "Text file output (optional)"
"""
)


def compare_pickle_with_mtz(params):
    """Compare intensities in pickle file with the scaled and merged intensities
  provided through the mtz file."""

    from libtbx.utils import Sorry
    from iotbx import mtz
    from dials.array_family import flex
    import cPickle as pickle

    m = mtz.object(params.reference)
    d = pickle.load(open(params.data))

    # fnd the right reference data
    mad = m.as_miller_arrays_dict(merge_equivalents=False)
    i = None
    for k in mad.keys():
        if k[2] == "IMEAN":
            i = mad[k].as_non_anomalous_array().expand_to_p1()
    assert not i is None

    # print reference information
    print("Reference")
    i.show_summary()

    # clean up - remove non-integrated data - TODO look at just strong spots
    if params.filter_integrated:
        d = d.select(d.get_flags(d.flags.integrated))
    if params.filter_partiality:
        d = d.select(d["partiality"] > params.filter_partiality)

    # apply scale factors from input file, if present

    # if requested, scale data to reference using simple kB model
    if params.kb_scale:
        raise Sorry("No can do, implementing kB scaling is on the to do")

    # match data to reference
    from cctbx import miller

    mi = miller.set(
        crystal_symmetry=i.crystal_symmetry(),
        anomalous_flag=False,
        indices=d["miller_index"],
    ).expand_to_p1()
    match = mi.match_indices(i)
    pairs = match.pairs()

    i1 = flex.double()
    i2 = flex.double()
    scl = flex.double()

    # TODO remove outliers here from the paired up list => do not need to
    # worry about them biasing the correlation

    for p0, p1 in pairs:
        i1.append(d["intensity.sum.value"][p0])
        i2.append(i.data()[p1])
        if "partiality" in d:
            scl.append(d["partiality"][p0])
        else:
            scl.append(1.0)
    corr = flex.linear_correlation(i1, i2)
    corr2 = flex.linear_correlation(i1 / scl, i2)
    print(
        "Correlation:", corr.coefficient(), corr2.coefficient(), pairs.size(), d.size()
    )

    if params.out:
        with open(params.out, "w") as f:
            for iis in zip(i1, i2, scl):
                f.write("%f %f %f\n" % iis)


def run(args):
    from dials.util.options import OptionParser
    import libtbx.load_env

    usage = "%s [options] data=integrated.pickle reference=reference.mtz" % (
        libtbx.env.dispatcher_name
    )

    parser = OptionParser(usage=usage, phil=phil_scope)

    params, options, args = parser.parse_args(
        show_diff_phil=True, return_unhandled=True
    )

    compare_pickle_with_mtz(params)


if __name__ == "__main__":
    import sys

    run(sys.argv[1:])
