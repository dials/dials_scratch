from __future__ import division
from __future__ import print_function


if __name__ == "__main__":

    import os
    import libtbx.load_env
    from glob import glob
    from dxtbx.imageset import ImageSetFactory

    # Check dials_regression is configured
    try:
        path = libtbx.env.dist_path("dials_regression")
    except KeyError as e:
        print("FAIL: dials_regression not configured")
        raise

    # Find the filenames
    template = os.path.join(path, "centroid_test_data", "centroid_*.cbf")
    filenames = glob(template)

    # Create the sequence
    sequence = ImageSetFactory.new(filenames)
    assert len(sequence) == 1
    sequence = sequence[0]

    # Get the models
    beam = sequence.get_beam()
    detector = sequence.get_detector()
    gonio = sequence.get_goniometer()
    scan = sequence.get_scan()
    print(beam)
    print(detector)
    print(gonio)
    print(scan)

    print("sequence: ", sequence)
    print("sequence indices: ", sequence.indices())
    print("sequence array range: ", sequence.get_array_range())

    # Get a sub sequence
    sub_sequence = sequence[3:6]
    print("sub_sequence: ", sub_sequence)
    print("sub_sequence indices: ", sub_sequence.indices())
    print("sub_sequence array range: ", sub_sequence.get_array_range())

    # Loop through sub sequence
    offset = sub_sequence.get_array_range()[0]
    for i, f in enumerate(sub_sequence):
        print("Image {0} has size {1}".format(i + offset, f.all()))

    # Extract a volume
    volume = sequence.to_array()
    print("sequence volume 1 size: ", volume.all())

    volume = sequence.to_array((2, 7))
    print("sequence volume 2 size: ", volume.all())

    volume = sequence[2:7].to_array()
    print("sequence volume 3 size: ", volume.all())

    volume = sequence.to_array((2, 7, 100, 200, 100, 300))
    print("sequence volume 4 size: ", volume.all())

    print("Format: ", sequence.reader().get_format())
    print("Template: ", sequence.get_template())
