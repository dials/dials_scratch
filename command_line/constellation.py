from __future__ import division, print_function, absolute_import

import numpy
import open3d
import random

from dials.array_family import flex

from dials.util.options import OptionParser
from dials.util.options import flatten_experiments, flatten_reflections


def point_cloud(reflections):
    pc = open3d.geometry.PointCloud()

    points = numpy.zeros((reflections.size(), 3), dtype=numpy.float64)

    xyz = reflections["rlp"].parts()

    for j in range(3):
        points[:, j] = xyz[j].as_numpy_array()

    pc.points = open3d.open3d.utility.Vector3dVector(points)

    pc.paint_uniform_color((0.0, 0.0, 0.0))

    return pc


def unroll_counts(reflections):
    """Really dumb implementation"""
    unrolled = flex.reflection_table()

    flags = flex.int()
    ids = flex.int()
    panels = flex.size_t()
    intensity_value = flex.double()
    intensity_variance = flex.double()
    xyz_value = flex.vec3_double()
    xyz_variance = flex.vec3_double()

    for j, sbox in enumerate(reflections["shoebox"]):
        flag = reflections["flags"][j]
        ident = reflections["id"][j]
        panel = reflections["panel"][j]
        bbox = sbox.bbox
        data = sbox.data
        for _k in range(bbox[4], bbox[5]):
            for _j in range(bbox[2], bbox[3]):
                for _i in range(bbox[0], bbox[1]):
                    counts = data[(_k - bbox[4], _j - bbox[2], _i - bbox[0])]
                    for c in range(int(round(counts))):
                        x = _i + random.random() - 0.5
                        y = _j + random.random() - 0.5
                        z = _k + random.random() - 0.5
                        xyz_value.append((x, y, z))
                        xyz_variance.append((1.0, 1.0, 1.0))
                        intensity_value.append(1)
                        intensity_variance.append(1)
                        flags.append(flag)
                        ids.append(ident)
                        panels.append(panel)
    unrolled["flags"] = flags
    unrolled["id"] = ids
    unrolled["panel"] = panels
    unrolled["intensity.sum.value"] = intensity_value
    unrolled["intensity.sum.variance"] = intensity_variance
    unrolled["xyzobs.px.value"] = xyz_value
    unrolled["xyzobs.px.variance"] = xyz_variance

    return unrolled


def map_to_reciprocal_space(reflections, experiment):
    imageset = experiment.imageset
    reflections.centroid_px_to_mm(imageset.get_detector(), scan=imageset.get_scan())
    reflections.map_centroids_to_reciprocal_space(
        detector=imageset.get_detector(),
        beam=imageset.get_beam(),
        goniometer=imageset.get_goniometer(),
    )


if __name__ == "__main__":

    parser = OptionParser(
        read_experiments=True, read_reflections=True, check_format=False
    )

    params, options = parser.parse_args(show_diff_phil=False)
    experiments = flatten_experiments(params.input.experiments)[0]
    reflections = flatten_reflections(params.input.reflections)[0]
    print(reflections.size())
    reflections = unroll_counts(reflections)
    print(reflections.size())

    map_to_reciprocal_space(reflections, experiments)

    pc = point_cloud(reflections)

    open3d.visualization.draw_geometries([pc])
