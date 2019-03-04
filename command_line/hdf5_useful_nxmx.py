from __future__ import print_function
import h5py
import numpy


def get_stupid_value(h5_obj):
    if h5_obj.ndim == 0:
        return h5_obj[()]
    if h5_obj.ndim == 1 and h5_obj.shape == (1,):
        return h5_obj[()][0].item()
    raise RuntimeError("not scalar or 1D length 1 array")


def get_depends_on_stack(h5_file, name):
    from cctbx.sgtbx import rt_mx
    from scitbx.matrix import sqr, col

    obj = h5_file[name]
    depends_on = obj["depends_on"][...]
    actual_depends_on = get_stupid_value(depends_on)

    # now iterate through the stack...
    transformation = actual_depends_on

    obj = h5_file[transformation]

    t_type = obj.attrs["transformation_type"]

    shift = col((0, 0, 0))

    while transformation != ".":
        obj = h5_file[transformation]
        t_type = obj.attrs["transformation_type"]
        if t_type == "translation":
            shift += get_stupid_value(obj) * col(obj.attrs["vector"])
        else:
            raise RuntimeError("panic")
        transformation = obj.attrs["depends_on"]

    return shift


def get_actual_pixel_offset(detector):
    from scitbx.matrix import col

    modules = []
    for thing in detector:
        obj = detector[thing]
        if "NX_class" in obj.attrs:
            if obj.attrs["NX_class"] == "NXdetector_module":
                modules.append(obj)
    if len(modules) == 0:
        raise RuntimeError("no modules found")

    if len(modules) > 1:
        raise RuntimeError("too many modules found")

    module = modules[0]
    offset = module["module_offset"]
    fast = col(module["fast_pixel_direction"].attrs["vector"])
    slow = col(module["slow_pixel_direction"].attrs["vector"])
    actual_offset = col(offset.attrs["offset"])
    return actual_offset, fast, slow


def get_derived_beam_centre(h5_file):
    # first find the detector

    detectors = []

    def visitor(name, obj):
        if isinstance(obj, h5py.Group):
            if "NX_class" in obj.attrs:
                if obj.attrs["NX_class"] == "NXdetector":
                    detectors.append(name)

    h5_file.visititems(visitor)

    detector = detectors[0]

    actual_detector = h5_file[detector]
    offset, fast, slow = get_actual_pixel_offset(actual_detector)
    shift = get_depends_on_stack(h5_file, detector)
    origin = offset + shift
    norm = fast.cross(slow)
    beam = origin - norm * origin.dot(norm)
    return beam.dot(fast), beam.dot(slow)


def hdf5_useful_nxmx(h5_file):

    # wavelength
    wavelength_locations = [
        "entry/instrument/beam/incident_wavelength",
        "entry/sample/beam/incident_wavelength",
    ]
    for wavelength_location in wavelength_locations:
        try:
            wavelength_dataset = h5_file[wavelength_location]
        except KeyError as e:
            continue
        print("Using %s for wavelength" % wavelength_location)
        ndim = wavelength_dataset.ndim
        shape = wavelength_dataset.shape
        print("Wavelength array: %d dimensions; shape = %s" % (ndim, str(shape)))
        wavelengths = wavelength_dataset[()]
        print("====>", wavelengths, wavelength_dataset.attrs["units"])

    # beam centre info
    beam_x_locations = [
        "entry/instrument/detector/beam_center_x",
        "entry/instrument/eiger/beam_centre_x",
    ]
    for beam_x_location in beam_x_locations:
        try:
            beam_x_dataset = h5_file[beam_x_location]
        except KeyError as e:
            continue
        print("Using %s for beam X" % beam_x_location)
        ndim = beam_x_dataset.ndim
        shape = beam_x_dataset.shape
        print(
            "Beam X (i.e. fast) array: %d dimensions; shape = %s" % (ndim, str(shape))
        )
        beam_xs = beam_x_dataset[()]
        print("====>", beam_xs, beam_x_dataset.attrs["units"])

    beam_y_locations = [
        "entry/instrument/detector/beam_center_y",
        "entry/instrument/eiger/beam_centre_y",
    ]
    for beam_y_location in beam_y_locations:
        try:
            beam_y_dataset = h5_file[beam_y_location]
        except KeyError as e:
            continue
        print("Using %s for beam Y" % beam_y_location)
        ndim = beam_y_dataset.ndim
        shape = beam_y_dataset.shape
        print(
            "Beam Y (i.e. slow) array: %d dimensions; shape = %s" % (ndim, str(shape))
        )
        beam_ys = beam_y_dataset[()]
        print("====>", beam_ys, beam_y_dataset.attrs["units"])

    beam_xy = get_derived_beam_centre(h5_file)
    print("Derived beam centre(mm)", beam_xy)

    distance_locations = [
        "entry/instrument/detector/detector_distance",
        "entry/instrument/detector_z/det_z",
    ]
    for distance_location in distance_locations:
        try:
            distance_dataset = h5_file[distance_location]
        except KeyError as e:
            continue
        print("Using %s for distance" % distance_location)
        ndim = distance_dataset.ndim
        shape = distance_dataset.shape
        print("Distance array: %d dimensions; shape = %s" % (ndim, str(shape)))
        distances = distance_dataset[()]
        print("====>", distances, distance_dataset.attrs["units"])


def main(filename):

    f = h5py.File(filename, "r")
    hdf5_useful_nxmx(f)


if __name__ == "__main__":
    import sys

    main(sys.argv[1])
