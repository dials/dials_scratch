from __future__ import division
from __future__ import print_function

import matplotlib

matplotlib.use("Agg")

from cctbx.array_family import flex

import iotbx.phil

help_message = """

Examples::

"""

phil_scope = iotbx.phil.parse(
    """
height = 50
  .type = float
  .help = "Height of cone around goniometer phi axis in mm"
radius = 18
  .type = float
  .help = "Radius cone around goniometer phi axis in mm"
angle = None
  .type = float
output {
  animation = None
    .type = path
}
"""
)


def run(args):

    from dials.util.options import OptionParser
    from dials.util.options import flatten_datablocks
    from dials.util.masking import GoniometerShadowMaskGenerator
    from libtbx.utils import Sorry
    import libtbx.load_env

    usage = "%s [options] experiments.json" % libtbx.env.dispatcher_name

    parser = OptionParser(
        usage=usage,
        phil=phil_scope,
        read_datablocks=True,
        check_format=False,
        epilog=help_message,
    )

    params, options = parser.parse_args(show_diff_phil=True)
    datablocks = flatten_datablocks(params.input.datablock)
    if len(datablocks) == 0:
        parser.print_help()
        return

    imagesets = []
    for datablock in datablocks:
        imagesets.extend(datablock.extract_imagesets())

    for imageset in imagesets:
        import math

        height = params.height  # mm
        radius = params.radius  # mm

        steps_per_degree = 10
        steps_per_degree = 1

        theta = (
            flex.double([range(360 * steps_per_degree)])
            * math.pi
            / 180
            * 1
            / steps_per_degree
        )
        y = radius * flex.cos(theta)  # x
        z = radius * flex.sin(theta)  # y
        x = flex.double(theta.size(), height)  # z

        coords = flex.vec3_double(zip(x, y, z))
        coords.insert(0, (0, 0, 0))

        gonio = imageset.get_goniometer()
        scan = imageset.get_scan()
        beam = imageset.get_beam()
        detector = imageset.get_detector()

        if params.angle is not None:
            angle = params.angle
        else:
            angle = scan.get_oscillation()[0]
        gonio_masker = GoniometerShadowMaskGenerator(
            gonio, coords, flex.size_t(len(coords), 0)
        )

        from matplotlib import pyplot as plt

        if params.output.animation is not None:

            import matplotlib.animation as manimation

            import os.path

            ext = os.path.splitext(params.output.animation)
            metadata = dict(
                title="Movie Test", artist="Matplotlib", comment="Movie support!"
            )
            if ext[1] == ".mp4":
                FFMpegWriter = manimation.writers["ffmpeg"]
                writer = FFMpegWriter(fps=15, metadata=metadata)
            elif ext[1] == ".gif":
                ImagemagickWriter = manimation.writers["imagemagick_file"]
                writer = ImagemagickWriter(fps=15, metadata=metadata)

            fig = plt.figure()
            l, = plt.plot([], [], c="r", marker=None)
            plt.axes().set_aspect("equal")
            plt.xlim(0, detector[0].get_image_size()[0])
            plt.ylim(0, detector[0].get_image_size()[0])
            plt.gca().invert_yaxis()
            title = plt.axes().set_title("")

            with writer.saving(fig, params.output.animation, 100):
                start, end = scan.get_array_range()
                step_size = 5
                for i in range(start, end, step_size):
                    angle = scan.get_angle_from_array_index(i)
                    shadow_boundary = gonio_masker.project_extrema(detector, angle)
                    x, y = shadow_boundary[0].parts()
                    l.set_data(x.as_numpy_array(), y.as_numpy_array())
                    title.set_text("scan_angle = %.1f degrees" % angle)
                    writer.grab_frame()

            plt.close()

        shadow_boundary = gonio_masker.project_extrema(detector, angle)

        with open("shadow.phil", "wb") as f:
            print("untrusted {", file=f)
            print("  polygon = \\", file=f)
            for c in shadow_boundary[0]:
                print("    %0.f %.0f \\" % (max(c[0], 0), max(c[1], 0)), file=f)
            print("}", file=f)

        import matplotlib.pyplot as plt

        fig = plt.figure()

        x, y, z = coords.parts()
        plt.scatter(x.as_numpy_array(), y.as_numpy_array())
        plt.axes().set_aspect("equal")
        plt.xlabel("x (gonio axis)")
        plt.ylabel("y (perpendicular to beam)")
        plt.savefig("gonio_xy.png")

        plt.scatter(y.as_numpy_array(), z.as_numpy_array())
        plt.axes().set_aspect("equal")
        plt.xlabel("y (perpendicular to beam)")
        plt.ylabel("z (towards beam))")
        plt.savefig("gonio_yz.png")

        plt.scatter(z.as_numpy_array(), x.as_numpy_array())
        plt.axes().set_aspect("equal")
        plt.xlabel("z (towards beam)")
        plt.ylabel("x (gonio axis)")
        plt.savefig("gonio_zx.png")

        for p_id in range(len(detector)):
            x, y = shadow_boundary[p_id].parts()
            fig = plt.figure()
            plt.scatter(x.as_numpy_array(), y.as_numpy_array(), c="r", s=1, marker="x")
            plt.axes().set_aspect("equal")
            plt.xlim(0, detector[p_id].get_image_size()[0])
            plt.ylim(0, detector[p_id].get_image_size()[0])
            plt.gca().invert_yaxis()
            plt.savefig("shadow.png")


if __name__ == "__main__":
    import sys

    run(sys.argv[1:])
