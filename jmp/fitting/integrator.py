from __future__ import division
from __future__ import print_function
import dials


def read_experiments(filename):
    from dxtbx.model.experiment_list import ExperimentListFactory

    return ExperimentListFactory.from_json_file(filename)


def read_reference(filename):
    import cPickle as pickle

    return pickle.load(open(filename))


def read_reflections(filename):
    from dials.array_family import flex

    r = flex.reflection_table()
    for f in filename:
        print("Reading %s" % f)
        r.extend(flex.reflection_table.from_pickle(f))
    return r


def integrate(experiments, reflections, reference, grid_size=5, detector_space=False):
    from dials.algorithms.profile_model.modeller import CircleSampler
    from dials.array_family import flex
    from dials.algorithms.profile_model.gaussian_rs.transform import TransformReverse
    from dials.algorithms.profile_model.gaussian_rs.transform import TransformForward
    from dials.algorithms.profile_model.gaussian_rs.transform import (
        TransformReverseNoModel,
    )
    from dials.algorithms.profile_model.gaussian_rs.transform import TransformSpec
    from dials.algorithms.profile_model.gaussian_rs import CoordinateSystem
    from dials.algorithms.integration.fit import ProfileFitter

    selection = reflections.get_flags(reflections.flags.integrated_sum)
    reflections = reflections.select(selection)
    print("Selected %d reflections to integrate" % len(reflections))

    assert len(reference) % 9 == 0
    num_scan_points = len(reference) // 9

    sampler = CircleSampler(
        experiments[0].detector[0].get_image_size(),
        experiments[0].scan.get_array_range(),
        num_scan_points,
    )

    spec = TransformSpec(
        experiments[0].beam,
        experiments[0].detector,
        experiments[0].goniometer,
        experiments[0].scan,
        experiments[0].profile.sigma_b(deg=False),
        experiments[0].profile.sigma_m(deg=False),
        experiments[0].profile.n_sigma() * 1.5,
        grid_size,
    )

    m2 = experiments[0].goniometer.get_rotation_axis()
    s0 = experiments[0].beam.get_s0()

    Iprf = flex.double(len(reflections))
    Vprf = flex.double(len(reflections))
    Cprf = flex.double(len(reflections))
    Fprf = flex.bool(len(reflections))
    Part = reflections["partiality"]

    reflections["intensity.prf_old.value"] = reflections["intensity.prf.value"]
    reflections["intensity.prf_old.variance"] = reflections["intensity.prf.variance"]

    selection = reflections.get_flags(reflections.flags.integrated_prf)

    reflections.unset_flags(~Fprf, reflections.flags.integrated_prf)

    for i, r in enumerate(reflections):

        if selection[i] == False:
            continue

        s1 = r["s1"]
        phi = r["xyzcal.mm"][2]
        xyz = r["xyzcal.px"]
        bbox = r["bbox"]
        panel = r["panel"]
        image = r["shoebox"].data.as_double()
        background = r["shoebox"].background.as_double()
        mask = r["shoebox"].mask.as_1d() == 5  # | (r['shoebox'].mask.as_1d() == 3)
        mask.reshape(image.accessor())
        cs = CoordinateSystem(m2, s0, s1, phi)

        index = sampler.nearest(0, xyz)

        profile, profile_mask = reference[index]

        # print flex.sum(profile)
        # print r['partiality']

        if detector_space:

            transform = TransformReverseNoModel(spec, cs, bbox, panel, profile)
            p = transform.profile()
            d = image
            m = mask
            b = background
            # print flex.sum(p)
            Part[i] = flex.sum(p)
            # ysize, xsize = p.all()[1:3]

            # p1 = flex.double(flex.grid(1, ysize , xsize))
            # d1 = flex.double(flex.grid(1, ysize , xsize))
            # b1 = flex.double(flex.grid(1, ysize , xsize))
            # m1 = flex.double(flex.grid(1, ysize , xsize))
            # for k in range(p.all()[0]):
            #   p1 += p[k:k+1,:,:]
            #   d1 += d[k:k+1,:,:]
            #   b1 += b[k:k+1,:,:]
            #   m1 = m[k:k+1,:,:]

            try:

                fit = ProfileFitter(d, b, m, p, 1e-3, 100)
                assert fit.niter() < 100
                Iprf[i] = fit.intensity()[0]
                Vprf[i] = fit.variance()[0]
                Cprf[i] = fit.correlation()
                Fprf[i] = True
                if r["intensity.sum.value"] > 10 and abs(fit.intensity()) < 1e-3:
                    print(
                        r["miller_index"],
                        i,
                        fit.intensity(),
                        r["intensity.sum.value"],
                        r["intensity.prf_old.value"],
                        Part[i],
                        fit.niter(),
                    )
                # from matplotlib import pylab
                # pylab.imshow(p1.as_numpy_array()[0,:,:], interpolation='none')
                # pylab.show()
            except Exception as e:
                print(e)
                pass

        else:

            try:

                transform = TransformForward(
                    spec, cs, bbox, panel, image, background, mask
                )

                p = profile
                d = transform.profile()
                b = transform.background()
                m = transform.mask() & profile_mask

                # if r['miller_index'] == (9, -25, 74):
                #   print list(p)
                #   print list(m)
                #   print list(b)
                #   print list(d)
                #   exit(0)

                fit = ProfileFitter(d, b, m, p, 1e-3, 100)
                assert fit.niter() < 100
                Iprf[i] = fit.intensity()[0]
                Vprf[i] = fit.variance()[0]
                Cprf[i] = fit.correlation()
                Fprf[i] = True
                print(
                    r["miller_index"], i, fit.intensity(), r["intensity.prf_old.value"]
                )
                # from matplotlib import pylab
                # pylab.imshow(p1.as_numpy_array()[0,:,:], interpolation='none')
                # pylab.show()
            except Exception:
                pass

    reflections["intensity.prf.value"] = Iprf
    reflections["intensity.prf.variance"] = Vprf
    reflections["intensity.prf.correlation"] = Cprf
    reflections.set_flags(Fprf, reflections.flags.integrated_prf)

    return reflections


if __name__ == "__main__":
    from os.path import join
    import sys

    experiments_filename = sys.argv[1]
    reference_filename = sys.argv[2]
    grid_size = int(sys.argv[3])
    if sys.argv[4] == "True":
        detector_space = True
    elif sys.argv[4] == "False":
        detector_space = False
    else:
        raise RuntimeError(sys.argv[4])
    reflections_filename = sorted(sys.argv[5:])

    experiments = read_experiments(experiments_filename)
    reflections = read_reflections(reflections_filename)
    reference = read_reference(reference_filename)

    print("Read %d reflections" % len(reflections))

    reflections = integrate(
        experiments,
        reflections,
        reference[0],
        grid_size=grid_size,
        detector_space=detector_space,
    )

    del reflections["shoebox"]

    reflections.as_pickle("integrated.pickle")
