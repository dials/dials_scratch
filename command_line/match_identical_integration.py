import sys

from dials.array_family import flex
from annlib_ext import AnnAdaptorSelfInclude


def match(a, b, max_separation=2, key="xyzobs.px.value", scale=(1, 1, 1)):
    """Match reflections from list a and list b, returning tuple of
    flex.size_t indices, optionally specifying the maximum distance and
    key to search on (which is assumed to be a 3-vector column). Can also
    apply relative scales to the vectors before matching in case e.g. very
    wide or very fine slicing."""

    xyz_a = a[key]
    xyz_b = b[key]

    if scale != (1, 1, 1):
        x, y, z = xyz_a.parts()
        x *= scale[0]
        y *= scale[1]
        z *= scale[2]
        xyz_a = flex.vec3_double(x, y, z)

        x, y, z = xyz_b.parts()
        x *= scale[0]
        y *= scale[1]
        z *= scale[2]
        xyz_b = flex.vec3_double(x, y, z)

    a = xyz_a.as_double().as_1d()
    b = xyz_b.as_double().as_1d()
    ann = AnnAdaptorSelfInclude(a, 3)
    ann.query(b)

    mm = flex.size_t(range(xyz_b.size()))
    nn, distance = ann.nn, flex.sqrt(ann.distances)

    sel = distance <= max_separation

    mm = mm.select(sel)
    nn = nn.select(sel)
    distance = distance.select(sel)
    return nn, mm, distance


def compare(file1, file2):
    """Match the reflections between file1, file2, extract intensities and
    compare."""
    data1 = flex.reflection_table.from_file(file1)
    data1 = data1.select(data1.get_flags(data1.flags.integrated))
    data2 = flex.reflection_table.from_file(file2)
    data2 = data2.select(data2.get_flags(data2.flags.integrated))
    nn, mm, distance = match(data1, data2)

    assert flex.max(distance) == 0

    i1, v1 = data1["intensity.prf.value"], data1["intensity.prf.variance"]
    i2, v2 = data2["intensity.prf.value"], data2["intensity.prf.variance"]

    for n, m in zip(nn, mm):
        assert i1[n] - i2[m] == 0
        assert v1[n] - v2[m] == 0


compare(sys.argv[1], sys.argv[2])
