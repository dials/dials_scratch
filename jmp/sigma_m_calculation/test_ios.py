from __future__ import print_function


def get_reflections(experiments):

    from dials.algorithms.profile_model.gaussian_rs import zeta_factor
    from dials.algorithms.spot_prediction import PixelToMillerIndex
    from collections import defaultdict
    from dials.algorithms.spot_prediction import ScanStaticRayPredictor
    from math import floor, sqrt, pi
    from dials.array_family import flex

    reflection_pixels = defaultdict(list)

    transform = PixelToMillerIndex(
        experiments[0].beam,
        experiments[0].detector,
        experiments[0].goniometer,
        experiments[0].scan,
        experiments[0].crystal,
    )

    xsize, ysize = experiments[0].detector[0].get_image_size()
    z0 = experiments[0].scan.get_array_range()[0]
    image = experiments[0].imageset[0][0]
    mask = flex.bool(image.accessor())

    for j in range(ysize):
        print(j)
        for i in range(xsize):
            h = transform.h(0, i + 0.5, j + 0.5, z0 + 0.5)

            hkl = tuple(map(lambda x: int(floor(x + 0.5)), h))

            d = sqrt(sum(map(lambda a: (a[0] - a[1]) ** 2, zip(h, hkl))))

            if d < 0.3:
                foreground = True
            else:
                foreground = False

            mask[j, i] = foreground

            reflection_pixels[hkl].append((j, i, foreground))

    # from matplotlib import pylab
    # pylab.imshow(mask.as_numpy_array())
    # pylab.show()

    I = []
    V = []
    T = []

    predictor = ScanStaticRayPredictor(
        experiments[0].beam.get_s0(),
        experiments[0].goniometer.get_rotation_axis(),
        experiments[0].goniometer.get_fixed_rotation(),
        experiments[0].goniometer.get_setting_rotation(),
        (-2 * pi, 2 * pi),
    )

    UB = experiments[0].crystal.get_A()
    phi0 = experiments[0].scan.get_angle_from_array_index(z0 + 0.5, deg=False)
    m2 = experiments[0].goniometer.get_rotation_axis()
    s0 = experiments[0].beam.get_s0()

    for hkl, pixel_list in reflection_pixels.iteritems():

        rays = predictor(hkl, UB)
        if len(rays) == 0:
            continue
        elif len(rays) == 1:
            dphi = rays[0].angle - phi0
        else:
            dphi0 = ((rays[0].angle - phi0) + pi) % (2 * pi) - pi
            dphi1 = ((rays[1].angle - phi0) + pi) % (2 * pi) - pi
            if abs(dphi0) < abs(dphi1):
                dphi = dphi0
                s1 = rays[0].s1
            else:
                dphi = dphi1
                s1 = rays[1].s1

        if abs(dphi) > 5 * pi / 180.0:
            continue

        try:
            zeta = zeta_factor(m2, s0, s1)
        except Exception:
            continue

        dphi *= zeta

        I_sum = 0
        B_sum = 0
        I_count = 0
        B_count = 0
        for j, i, foreground in pixel_list:
            data = image[j, i]
            if foreground:
                I_sum += data
                I_count += 1
            else:
                B_sum += data
                B_count += 1
        B = B_sum / B_count
        I.append(I_sum - B * I_count)
        V.append(I_sum + B * I_count * (1 + I_count / B_count))
        T.append(dphi * 180 / pi)

    print(min(T), max(T))

    IOS = []
    TN = []
    for i, v, t in zip(I, V, T):
        if i > 0 and v > 0:
            IOS.append(i / sqrt(v))
            TN.append(t)

    from matplotlib import pylab

    pylab.hist(TN, bins=20, weights=IOS)
    pylab.xlabel("DPHI (degrees)")
    pylab.show()


if __name__ == "__main__":

    from dxtbx.model.experiment_list import ExperimentListFactory

    experiments = ExperimentListFactory.from_json_file("experiments.json")

    reflections = get_reflections(experiments)
