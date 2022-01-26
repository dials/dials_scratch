import sys
from collections import namedtuple

import pycbf

from scitbx import matrix


def axis_rt(axis, setting=None):
    if axis.type == b"rotation":
        if setting:
            a = setting.angle_start_range_incr[0]
        else:
            a = 0.0
        v = matrix.col(axis.vector)
        r = v.axis_and_angle_as_r3_rotation_matrix(a, deg=True)
        t = matrix.col(axis.offset)
    elif axis.type == b"translation":
        if setting:
            t = setting.displacement_start_range_incr[0]
        else:
            t = 0.0
        v = matrix.col(axis.vector)
        r = matrix.sqr((1, 0, 0, 0, 1, 0, 0, 0, 1))
        o = matrix.col(axis.offset)
        t = v * t + o

    rt = matrix.rt((r, t))
    return rt


def get_fast_slow_axis_names(handle):
    axes = {}
    handle.find_category(b"array_structure_list")
    for i in range(handle.count_rows()):
        handle.find_column(b"precedence")
        # convert to 0-based counting
        p = int(handle.get_value()) - 1
        handle.find_column(b"axis_set_id")
        a = handle.get_value()
        axes[p] = a
        handle.next_row()

    for k in 0, 1:
        assert k in axes

    return axes


def get_canonical_fast_slow_origin(handle):
    detector = handle.construct_detector(0)

    _origin = matrix.col(detector.get_pixel_coordinates(0, 0))
    _fast = matrix.col(detector.get_detector_axis_fast())
    _slow = matrix.col(detector.get_detector_axis_slow())

    detector.__swig_destroy__(detector)
    del detector

    return _fast, _slow, _origin


def handle_to_settings(handle):

    # simple namedtuple to define the structure for settings
    setting = namedtuple(
        "setting", ("id", "angle_start_range_incr", "displacement_start_range_incr")
    )

    settings = []
    handle.find_category(b"diffrn_scan_axis")
    for i in range(handle.count_rows()):
        handle.find_column(b"axis_id")
        _id = handle.get_value()
        handle.find_column(b"angle_start")
        _angle_start = float(handle.get_value())
        handle.find_column(b"angle_range")
        _angle_range = float(handle.get_value())
        handle.find_column(b"angle_increment")
        _angle_incr = float(handle.get_value())

        handle.find_column(b"displacement_start")
        _displacement_start = float(handle.get_value())
        handle.find_column(b"displacement_range")
        _displacement_range = float(handle.get_value())
        handle.find_column(b"displacement_increment")
        _displacement_incr = float(handle.get_value())

        settings.append(
            setting(
                _id,
                (_angle_start, _angle_range, _angle_incr),
                (_displacement_start, _displacement_range, _displacement_incr),
            )
        )

        handle.next_row()
    return settings


def handle_to_axis_stack(handle):
    # simple namedtuple to define axes
    axis = namedtuple(
        "axis", ("id", "type", "equipment", "depends_on", "vector", "offset")
    )

    axes = []

    handle.find_category(b"axis")
    for i in range(handle.count_rows()):
        handle.find_column(b"equipment")
        _equipment = handle.get_value()
        handle.find_column(b"id")
        _id = handle.get_value()
        handle.find_column(b"type")
        _type = handle.get_value()
        _vector = []
        for i in range(3):
            handle.find_column(b"vector[%i]" % (i + 1))
            _vector.append(float(handle.get_value()))
        _offset = []
        for i in range(3):
            handle.find_column(b"offset[%i]" % (i + 1))
            try:
                _offset.append(float(handle.get_value()))
            except ValueError:
                _offset.append(0.0)
        handle.find_column(b"depends_on")
        _depends_on = handle.get_value()
        axes.append(
            axis(_id, _type, _equipment, _depends_on, tuple(_vector), tuple(_offset))
        )
        handle.next_row()

    return axes


def draw(cbf_filename):
    handle = pycbf.cbf_handle_struct()
    handle.read_widefile(cbf_filename.encode(), pycbf.MSG_DIGEST)

    # read out all the data structures
    stack = handle_to_axis_stack(handle)
    settings = handle_to_settings(handle)

    # construct a tree from these
    tree = {}
    ht = {}
    for s in stack:
        ht[s.id] = s
        tree[s.id] = s.depends_on

    # convenient lookup for axis settings - means we have an indexed table
    sets = {}
    for s in settings:
        sets[s.id] = s

    # the leaves are the things at the top of the tree i.e. those things
    # on which nothing depends
    leaves = set(s.id for s in stack) - set(tree.values())

    # now for each object in the stack show what it depends on - this should be
    # well defined until you hit the ground at "."
    for s in stack:
        t = s.id
        d = [t]
        while t != b".":
            t = tree[t]
            d.append(t)
        print("\n--> ".join(_d.decode() for _d in d))

    # build up a handy dandy hash table of rt matrices
    rts = {}
    for s in stack:
        if not s.equipment in (b"goniometer", b"detector"):
            continue

        rts[s.id] = axis_rt(s, sets.get(s.id, None))

    # derive the final fast, slow, origin for the detector

    axes = get_fast_slow_axis_names(handle)
    fast_name = axes[0]
    slow_name = axes[1]

    fast = matrix.col(ht[fast_name].vector)
    depends = ht[fast_name].depends_on
    while depends != b".":
        fast = rts[depends].r * fast
        depends = ht[depends].depends_on

    slow = matrix.col(ht[slow_name].vector)
    depends = ht[slow_name].depends_on
    while depends != b".":
        slow = rts[depends].r * slow
        depends = ht[depends].depends_on

    origin = matrix.col(ht[fast_name].offset)
    depends = ht[fast_name].depends_on
    while depends != b".":
        origin = rts[depends] * origin
        depends = ht[depends].depends_on

    # verify status

    _fast, _slow, _origin = get_canonical_fast_slow_origin(handle)

    assert (origin - _origin).length() < 1e-12
    assert (fast - _fast).length() < 1e-12
    assert (slow - _slow).length() < 1e-12

    # now show something useful - goniometer and detector axis settings and
    # types - this should be sufficient to reconstruct the full geometry

    for s in stack:
        if not s.equipment in (b"goniometer", b"detector"):
            continue

        rt = axis_rt(s, sets.get(s.id, None))

        try:
            # get the axis setting - either as a translation or a rotation
            if s.type == b"translation":
                a = sets[s.id].displacement_start_range_incr[0]
            else:
                a = sets[s.id].angle_start_range_incr[0]
        except KeyError:
            a = 0.0

        print(s.id.decode(), a)


if __name__ == "__main__":
    draw(sys.argv[1])
