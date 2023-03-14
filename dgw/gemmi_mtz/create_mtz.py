import sys
from dials.array_family import flex
from dials_gemmi_mtz_ext import create_mtz

if __name__ == "__main__":

    integrated = flex.reflection_table.from_file(sys.argv[1])

    # add a batch column
    zdet = integrated["xyzobs.px.value"].parts()[2]
    integrated["batch"] = flex.floor(zdet).iround() + 1

    create_mtz("foo", integrated)
