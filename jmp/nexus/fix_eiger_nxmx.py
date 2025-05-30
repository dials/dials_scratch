from __future__ import print_function


class EigerNXmxFixer(object):
    """
    A hacky class to read an NXmx file

    """

    def __init__(self, input_filename, output_filename):
        import h5py
        import shutil
        from os.path import join
        from scitbx import matrix

        #
        print("Copying %s to %s to fix NXmx file" % (input_filename, output_filename))
        shutil.copy(input_filename, output_filename)

        #
        print("Opening %s to perform fixes" % (output_filename))
        handle = h5py.File(output_filename, "r+")

        # Add some simple datasets
        def create_scalar(handle, path, dtype, value):
            print(
                "Adding dataset %s with value %s"
                % (join(handle.name, path), str(value))
            )
            dataset = handle.create_dataset(path, (), dtype=dtype)
            dataset[()] = value

        # Add NXmx definition
        create_scalar(handle["entry"], "definition", "S4", "NXmx")

        # Add saturation value
        create_scalar(
            handle["entry/instrument/detector"], "saturation_value", "int32", 10000
        )

        # Add detector type
        create_scalar(handle["entry/instrument/detector"], "type", "S4", "PIXEL")

        # Move the beam
        print("Copying /entry/instrument/beam to /entry/sample/beam")
        handle.copy("/entry/instrument/beam", "/entry/sample/beam")

        # Create detector module
        module_path = "/entry/instrument/detector/module"
        print("Creating detector module %s" % (module_path))
        group = handle.create_group(module_path)
        group.attrs["NX_class"] = "NXdetector_module"

        # Add a module index
        create_scalar(group, "module_index", "int64", 0)

        # Create detector data origin
        print(
            "Adding dataset %s with value %s"
            % (join(group.name, "data_origin"), str((0, 0)))
        )
        dataset = group.create_dataset("data_origin", (2,), dtype="int32")
        dataset[0] = 0
        dataset[1] = 0

        # Create detector data size
        print(
            "Adding dataset %s with value %s"
            % (join(group.name, "data_size"), str((1030, 1065)))
        )
        dataset = group.create_dataset("data_size", (2,), dtype="int32")
        dataset[0] = handle["/entry/data/data_000001"].shape[2]
        dataset[1] = handle["/entry/data/data_000001"].shape[1]

        # Add fast_pixel_size dataset
        print(
            "Using /entry/instrument/detector/geometry/orientation/value as fast/slow pixel directions"
        )
        fast_axis = handle["/entry/instrument/detector/geometry/orientation/value"][0:3]
        slow_axis = handle["/entry/instrument/detector/geometry/orientation/value"][3:6]
        create_scalar(
            group,
            "fast_pixel_direction",
            "float32",
            handle["/entry/instrument/detector/x_pixel_size"].value,
        )
        group["fast_pixel_direction"].attrs["transformation_type"] = "translation"
        group["fast_pixel_direction"].attrs["vector"] = fast_axis
        group["fast_pixel_direction"].attrs["offset"] = 0
        group["fast_pixel_direction"].attrs["units"] = "m"
        group["fast_pixel_direction"].attrs[
            "depends_on"
        ] = "/entry/instrument/detector/transformations/translation"

        # Add slow_pixel_size dataset
        create_scalar(
            group,
            "slow_pixel_direction",
            "float32",
            handle["/entry/instrument/detector/y_pixel_size"].value,
        )
        group["slow_pixel_direction"].attrs["transformation_type"] = "translation"
        group["slow_pixel_direction"].attrs["vector"] = slow_axis
        group["slow_pixel_direction"].attrs["offset"] = 0
        group["slow_pixel_direction"].attrs["units"] = "m"
        group["slow_pixel_direction"].attrs[
            "depends_on"
        ] = "/entry/instrument/detector/transformations/translation"

        # Add module offset dataset
        print("Set module offset to be zero relative to detector")
        create_scalar(group, "module_offset", "float32", 0)
        group["module_offset"].attrs["transformation_type"] = "translation"
        group["module_offset"].attrs["vector"] = (0, 0, 0)
        group["module_offset"].attrs["offset"] = 0
        group["module_offset"].attrs["units"] = "m"
        group["module_offset"].attrs[
            "depends_on"
        ] = "/entry/instrument/detector/transformations/translation"

        # Create detector depends_on
        create_scalar(
            handle["/entry/instrument/detector"],
            "depends_on",
            "S1",
            "/entry/instrument/detector/transformations/translation",
        )

        # Add detector position
        print(
            "Using /entry/instrument/detector/geometry/translation/distances as transformation"
        )
        detector_offset_vector = matrix.col(
            handle["/entry/instrument/detector/geometry/translation/distances"][()]
        )
        group = handle.create_group("/entry/instrument/detector/transformations")
        group.attrs["NX_class"] = "NXtransformations"
        create_scalar(group, "translation", "float32", detector_offset_vector.length())
        group["translation"].attrs["transformation_type"] = "translation"
        group["translation"].attrs["vector"] = detector_offset_vector.normalize()
        group["translation"].attrs["offset"] = 0
        group["translation"].attrs["units"] = "m"
        group["translation"].attrs["depends_on"] = "."

        # Create goniometer transformations
        print("Creating group /entry/sample/transformation")
        group = handle.create_group("/entry/sample/transformations")
        group.attrs["NX_class"] = "NXtransformations"

        print("Creating omega transformation:")
        print(" - making up rotation axis to be (1, 0, 0)")
        print(" - making up starting angle to be 0")
        print(
            " - using /entry/sample/goniometer/omega_range_average as oscillation range"
        )

        # Get the number of images
        num_images = 0
        for name in handle["/entry/data"].iterkeys():
            num_images += len(handle[join("/entry/data", name)])
        dataset = group.create_dataset("omega", (num_images,), dtype="float32")
        dataset.attrs["units"] = "degree"
        dataset.attrs["transformation_type"] = "rotation"
        dataset.attrs["vector"] = (1, 0, 0)
        dataset.attrs["offset"] = 0
        dataset.attrs["depends_on"] = "."
        omega_range_average = handle["/entry/sample/goniometer/omega_range_average"][()]
        omega_range_average = int(omega_range_average * 100 + 0.5) / 100.0
        for i in range(num_images):
            angle = omega_range_average * i
            dataset[i] = angle

        # Create sample depends_on
        create_scalar(
            handle["/entry/sample"],
            "depends_on",
            "S%d" % len(dataset.name),
            str(dataset.name),
        )


if __name__ == "__main__":

    import sys

    if len(sys.argv) != 2:
        print("Usage: python fix_eiger_nxmx.py FILE.nxs")
        exit(0)

    input_filename = sys.argv[1]
    output_filename = "fixed.nxs"

    fixer = EigerNXmxFixer(input_filename, output_filename)
