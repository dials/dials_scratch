from dials.array_family import flex
import cPickle as pickle


def fix_xy(reflections_in, reflections_out):
    reflections = pickle.load(open(reflections_in, "r"))

    # validate that input is from P12M@DLS => 24 panels

    assert flex.max(reflections["panel"]) == 23

    delta = 195 + 17
    y_offset = delta * reflections["panel"]

    # apply fixes
    x, y, z = reflections["xyzobs.px.value"].parts()
    y += y_offset.as_double()
    reflections["xyzobs.px.value"] = flex.vec3_double(x, y, z)

    x, y, z = reflections["xyzcal.px"].parts()
    y += y_offset.as_double()
    reflections["xyzcal.px"] = flex.vec3_double(x, y, z)

    # save - should probably do this "properly"
    pickle.dump(reflections, open(reflections_out, "wb"))

    return


if __name__ == "__main__":
    import sys

    fix_xy(sys.argv[1], sys.argv[2])
