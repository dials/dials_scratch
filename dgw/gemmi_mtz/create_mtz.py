import sys
from dials.array_family import flex
from dials_gemmi_mtz_ext import create_mtz

if __name__ == "__main__":

    refl = flex.reflection_table.from_file(sys.argv[1])
    refl = refl.select(refl.get_flags(refl.flags.integrated))

    if not all(
        refl.has_key(k)
        for k in (
            "miller_index",
            "intensity.sum.value",
            "intensity.sum.variance",
            "intensity.prf.value",
            "intensity.prf.variance",
            "xyzobs.px.value",
        )
    ):
        sys.exit("missing required columns")

    # add a batch column
    zdet = refl["xyzobs.px.value"].parts()[2]
    refl["batch"] = flex.floor(zdet).iround() + 1

    create_mtz("foo", refl)
