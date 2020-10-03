import sys
from collections import namedtuple

import h5py
import numpy

from scitbx import matrix


def axis_rt(axis, setting=None):
    if axis.type == b"rotation":
        if setting:
            if hasattr(setting.positions, "__iter__"):
                a = setting.positions[0]
            else:
                a = setting.positions
        else:
            a = 0.0
        v = matrix.col(axis.vector)
        r = v.axis_and_angle_as_r3_rotation_matrix(a, deg=True)
        t = matrix.col(axis.offset)
    elif axis.type == b"translation":
        if setting:
            if hasattr(setting.positions, "__iter__"):
                t = setting.positions[0]
            else:
                t = setting.positions
        else:
            t = 0.0
        v = matrix.col(axis.vector)
        r = matrix.sqr((1, 0, 0, 0, 1, 0, 0, 0, 1))
        o = matrix.col(axis.offset)
        t = v * t + o

    rt = matrix.rt((r, t))
    return rt


def validate_nxmx(f):
    """Checks that this is an NXmx file, raises appropriate exception if not."""

    if not b"/entry/definition" in f:
        raise ValueError("/entry/definition missing")

    definition = f[b"/entry/definition"][()]
    if definition != numpy.string_("NXmx"):
        raise ValueError("definition not NXmx: %s" % str(definition))


def mm_distance_or_deg_angle(d):
    """Return the distance here in mm"""

    if not b"transformation_type" in d.attrs:
        raise ValueError("dataset has no transformation_type")

    if not b"units" in d.attrs:
        raise ValueError("units missing from dataset")

    if d.attrs[b"transformation_type"] == numpy.string_("translation"):
        mm_scale = {numpy.string_("m"): 1000.0, numpy.string_("mm"): 1.0}.get(
            d.attrs[b"units"]
        )

        return mm_scale * d[()]
    elif d.attrs[b"transformation_type"] == numpy.string_("rotation"):
        if d.attrs[b"units"] != numpy.string_("deg"):
            raise ValueError("only degrees supported as rotation unit")
        return d[()]
    else:
        raise ValueError("transformation_type not translation or rotation")


def detector_fast_slow_origin(f):
    """Construct a depends_on tree for the detector, assumed to be at
    /entry/instrument/detector with either attribute or dataset depends_on"""

    detector = f[b"/entry/instrument/detector"]

    # Commentary:
    #
    # In theory the detector dependency hierarchy can be derived from the
    # properties of what the detector depends on however this is not true
    # in real life, as the detector depends on only one axis (currently) and
    # has no offsets so is incorrectly defined. The fast and slow pixel
    # directions _are_ however robustly defined -> use these to define the
    # hierarchy then...

    depends_on = None

    if b"depends_on" in detector:
        depends_on = detector[b"depends_on"][()]
    elif b"depends_on" in detector.attrs:
        depends_on = detector.atrs[b"depends_on"]
    else:
        raise ValueError("no depends_on found in /entry/instrument/detector")

    # fast and slow directions, work on the basis that these have the dependency
    # attributes and take no prisoners if this is false

    fast = detector[b"module/fast_pixel_direction"]
    slow = detector[b"module/slow_pixel_direction"]

    # assert for the moment that these _do not_ have independent offsets
    nil = numpy.array((0.0, 0.0, 0.0))

    foff = fast.attrs.get(b"offset", nil)
    if not numpy.array_equal(foff, nil):
        raise ValueError("fast axis has offset")

    soff = slow.attrs.get(b"offset", nil)
    if not numpy.array_equal(soff, nil):
        raise ValueError("slow axis has offset")

    axis = namedtuple(
        "axis", ("id", "type", "equipment", "depends_on", "vector", "offset")
    )

    setting = namedtuple("setting", ("id", "positions"))

    axes = []
    settings = []

    # I don't know if this is really a constraint or not, quite possibly
    # one could depend on the other...
    if fast.attrs[b"depends_on"] != slow.attrs[b"depends_on"]:
        raise ValueError("fast and slow axes depend on different axes")

    depends_on = fast.attrs[b"depends_on"]

    while depends_on != numpy.string_("."):
        element = f[depends_on]

        # believe all offsets in m - want mm so
        offset = 1000 * element.attrs.get(b"offset", nil)
        vector = element.attrs[b"vector"]
        values = mm_distance_or_deg_angle(element)

        axes.append(
            axis(
                depends_on,
                element.attrs[b"transformation_type"],
                b"detector",
                element.attrs[b"depends_on"],
                vector,
                offset,
            )
        )
        settings.append(setting(depends_on, values))

        depends_on = element.attrs[b"depends_on"]

    # construct a tree from these
    tree = {}
    ht = {}
    for a in axes:
        ht[a.id] = a
        tree[a.id] = a.depends_on

    # convenient lookup for axis settings - means we have an indexed table
    sets = {}
    for s in settings:
        sets[s.id] = s

    # build up a handy dandy hash table of rt matrices
    rts = {}
    for a in axes:
        if not a.equipment in (b"goniometer", b"detector"):
            continue

        rts[a.id] = axis_rt(a, sets.get(a.id, None))

    origin = matrix.col(nil)
    fast_vector = list(fast.attrs[b"vector"])
    slow_vector = list(slow.attrs[b"vector"])
    depends = fast.attrs[b"depends_on"]

    while depends != b".":
        origin = rts[depends] * origin
        fast_vector = rts[depends].r * fast_vector
        slow_vector = rts[depends].r * slow_vector
        depends = ht[depends].depends_on

    return matrix.col(fast_vector), matrix.col(slow_vector), matrix.col(origin)


def mcstas_to_imgcif(v):
    """Convert vector v to imgCIF frame from mcstas"""

    rot = matrix.sqr((-1, 0, 0, 0, 1, 0, 0, 0, -1))

    return rot * v


def draw(filename):

    with h5py.File(filename, "r") as f:
        # belt + braces here
        validate_nxmx(f)

        # construct the tree for the detector
        f, s, o = detector_fast_slow_origin(f)

        f = mcstas_to_imgcif(f)
        s = mcstas_to_imgcif(s)
        o = mcstas_to_imgcif(o)

    print(f.elems)
    print(s.elems)
    print(o.elems)


if __name__ == "__main__":
    draw(sys.argv[1])
