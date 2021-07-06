#!/usr/bin/env python3
# coding: utf-8

import subprocess
import sys

import h5py
import numpy as np


def categorise_datasets(key, value):
    """
    Find all data sets in a HDF5 file that scale with the number of images.

    For our purposes, we would ideally resize these data sets to have length 1.
    Since h5py cannot do that for data sets that are not chunked, we must instead
    delete them and recreate them, using only the first datum.  We cannot do that
    within a visititems call because the keys and values created by h5py's visititems
    are h5py._hl.base.ItemsViewHDF5 objects, which hold a lock on the parent h5py
    group that is being iterated over.  A subgroup or dataset therefore cannot be
    deleted from within a visit or visititems method.

    Instead we can populate a list of the data sets to be truncated.
    """
    if isinstance(value, h5py.Dataset):
        # We'll handle nimages separately.  The datafile links cannot be traversed by
        # visititems (thanks h5py!), so we'll have to deal with those separately too.
        if "nimages" in key or key.startswith("entry/data"):
            pass
        elif value.shape and value.shape[0] == nimages:
            scales_with_nimage.append(key)
        elif value.size * value.dtype.itemsize <= 1024:
            small_datasets.append(key)
        elif value.size * value.dtype.itemsize > 1024:
            large_datasets.append(key)
        else:
            print("Warning:  uncategorised data set '{}'.".format(key))
    elif isinstance(value, h5py.Group):
        groups.append(key)
    else:
        print("Warning:  uncategorised entry '{}'".format(key))


def reduce_size(
    master_file, target_file, large_data_keys, condition, dummy_values=None, dtype=None
):
    selected_keys = filter(condition, large_data_keys)
    # Reduce the size of the first matching data set by setting it to dummy values
    # that compress well and enabling compression (requires chunking).
    first_key = next(selected_keys)
    first_data = master_file[first_key][...]
    target_file.create_dataset(
        first_key,
        data=dummy_values(first_data) if dummy_values else first_data,
        chunks=first_data.shape,
        compression="gzip",
        dtype=dtype,
    )
    # Make all subsequent data sets hard links to the first.
    for key in selected_keys:
        target_file[key] = target_file[first_key]


master_file = sys.argv[1]
data_file = master_file.replace("master", "data_000001")
new_master_file = "dectris_eiger_master.h5"
new_data_file = "dectris_eiger_data_000001.h5"

# Truncate and compress the master file, with some dummy values for ease of compression.
with h5py.File(master_file, "r") as f, h5py.File("dummy_master.h5", "w") as g:
    # Find the number of images.
    nimages_key = "entry/instrument/detector/detectorSpecific/nimages"
    nimages = f[nimages_key][...]

    scales_with_nimage = []
    small_datasets = []
    large_datasets = []
    groups = []

    # Divide the various data sets into various categories.
    f.visititems(categorise_datasets)

    # Copy only the first data file.
    data_key = "entry/data/data_000001"
    g[data_key] = h5py.ExternalLink(new_data_file, "entry/data/data")

    # Truncate those data that scale with number of images to the first entry.
    for key in scales_with_nimage:
        value = f[key][:1]
        g.create_dataset(key, data=value, dtype=value.dtype)

    # Set the image count to 1.
    nimages_key = "entry/instrument/detector/detectorSpecific/nimages"
    g.create_dataset(nimages_key, data=1, dtype="u4")

    # Copy the small data sets over unmodified.
    for key in small_datasets:
        value = f[key]
        g.create_dataset_like(key, value, data=value)

    # Reduce the size of the count rate correction lookup tables.  Most of the values
    # are just (2**16 - 1) so it ought to compress well without dummy values.
    reduce_size(
        f, g, large_datasets, lambda key: "countrate_correction_lookup_table" in key
    )
    # Reduce the size of the countrate correction tables.  These are quite small so
    # we won't resort to dummy values.
    reduce_size(f, g, large_datasets, lambda key: "countrate_correction_table" in key)

    # Reduce the size of the detector module flat-field correction data sets,
    # using ones as dummy values.
    reduce_size(
        f,
        g,
        large_datasets,
        lambda key: "detectorModule" in key and "flatfield" in key,
        np.ones_like,
    )
    # Reduce the size of the detector module pixel masks, using zeroes as dummy values.
    # Dectris only uses five bits for the mask, so the default dtype, uint32,
    # is wasteful.  Use uint8 instead.
    reduce_size(
        f,
        g,
        large_datasets,
        lambda key: "detectorModule" in key and "pixel_mask" in key,
        np.zeros_like,
        np.uint8,
    )
    # Reduce the size of the detector module trimbit data sets (whatever those are).
    # Use 16 as a dummy value throughout (I'm guessing that this is appropriate...)
    reduce_size(
        f,
        g,
        large_datasets,
        lambda key: "trimbit" in key,
        lambda data: np.full_like(data, 16),
    )
    # Reduce the size of the global detector flat-field correction data set,
    # using ones as dummy values.
    reduce_size(
        f,
        g,
        large_datasets,
        lambda key: "detectorModule" not in key and "flatfield" in key,
        np.ones_like,
    )
    # Reduce the size of the global pixel mask, using zeroes as dummy values and a
    # more sensible dtype.
    reduce_size(
        f,
        g,
        large_datasets,
        lambda key: "detectorModule" not in key and "pixel_mask" in key,
        np.zeros_like,
        np.uint8,
    )

    # Copy the NX_class and other attributes from the old file to the new one.
    for key in (
        [nimages_key] + large_datasets + small_datasets + scales_with_nimage + groups
    ):
        for attribute, value in f[key].attrs.items():
            g[key].attrs.create(attribute, value)

# Compress the first data file.
with h5py.File(data_file, "r") as f, h5py.File("dummy_data.h5", "w") as g:
    image_key = "entry/data/data"
    value = f[image_key]
    g.create_dataset(image_key, data=value, chunks=value.shape, compression="gzip")
    # Keep only the first image.
    g[image_key].resize((1, *g[image_key].shape[1:]))

    # Get all the 'NX_class' and other attributes from the old file to the new one.
    data_file_keys = []
    g.visit(lambda key: data_file_keys.append(key))
    for key in data_file_keys:
        for attribute, value in f[key].attrs.items():
            g[key].attrs.create(attribute, value)
    # The number of images is recorded (redundantly) in the attributes of the images
    # data set.   Change its value there too.
    if "image_nr_high" in g[image_key].attrs:
        g[image_key].attrs.modify("image_nr_high", 1)


repack_commands = {
    master_file: "h5repack dummy_master.h5 {} && rm dummy_master.h5".format(
        new_master_file
    ),
    data_file: "h5repack dummy_data.h5 {} && rm dummy_data.h5".format(new_data_file),
}

for filename, command in repack_commands.items():
    print("Repacking {}".format(filename))
    subprocess.call(command, shell=True)
