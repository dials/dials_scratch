# Write a new EIGER HDF5 file with selected frames
#
#  By Takanori Nakane, IPR, Osaka University, 2024-2025
#     tnakane.protein@osaka-u.ac.jp

import h5py

# DO NOT IMPORT hdf5plugin!!!
import numpy as np
import os
import sys

if len(sys.argv) != 4:
    sys.stderr.write(
        "Usage: dials.python select_eiger_frames.py input_master.h5 output_master.h5 hits.lst\n\n"
    )
    sys.stderr.write(
        "hits.lst is the indices of hit images to be selected (1-indexed).\n"
    )
    sys.stderr.write("One way to generate it is as follows:\n")
    sys.stderr.write(
        "    dials.spot_counts_per_image imported.expt strong.refl json=strong.json > count.lst\n"
    )
    sys.stderr.write(
        "    awk 'NF==9 && 0+$4 >= 15 {print $2}' count.lst > hits.lst # 15 is the minimum number of spots\n\n"
    )
    sys.exit(-1)

f = h5py.File(sys.argv[1], "r")
g = h5py.File(sys.argv[2], "w")
hits = np.loadtxt(sys.argv[3], dtype=int).reshape(-1) - 1  # to 0-indexed

# Dr. Naomine Yano (JASRI) pointed out that fixed target serial datasets
# were collected using multiple triggers. Thus, we have to multiply
# the number of frames per trigger with the number of triggers to get the
# total frame number.
n_frames_per_trigger = f["/entry/instrument/detector/detectorSpecific/nimages"][()]
n_triggers = f["/entry/instrument/detector/detectorSpecific/ntrigger"][()]
n_frames = n_frames_per_trigger * n_triggers

frame_to_data = [None] * n_frames
print("Number of frames per trigger", n_frames_per_trigger)
print("Number of triggers", n_triggers)
print("Number of total frames in the input:", n_frames)
print("Number of hits to write:", len(hits))

for dataname in sorted(f["/entry/data"]):
    data = f["/entry/data/" + dataname]
    f_start = (
        int(data.attrs["image_nr_low"]) - 1
    )  # originally 1-indexed in the master file,
    f_end = int(data.attrs["image_nr_high"]) - 1  # so we make them 0-indexed
    #   print(dataname, f_start, f_end)
    for i in range(f_start, f_end + 1):
        frame_to_data[i] = (data, i - f_start)

# Amend and Copy metadata


def copy_visitor(name, node):
    if name == "entry/instrument/detector/detectorSpecific/nimages":
        g[name] = len(hits)
    elif name == "entry/instrument/detector/detectorSpecific/ntrigger":
        g[name] = 1  # put everything in one trigger (thanks to Dr. Yano@JASRI)
    elif (
        name == "entry/sample/goniometer/omega"
        or name == "entry/sample/goniometer/omega_end"
        or name == "entry/sample/transformations/omega"
        or name == "entry/sample/transformations/omega_end"
    ):
        filtered = node[()][hits]
        new_data = g.create_dataset(
            name, data=filtered, compression="gzip", shuffle=True
        )
        for k, v in node.attrs.items():
            new_data.attrs[k] = v
    elif isinstance(node, h5py.Dataset):
        f.copy(name, g, name)
    else:  # group
        new_group = g.require_group(name)
        for k, v in node.attrs.items():
            new_group.attrs[k] = v
        if name == "entry/data":
            return
        # Because old h5py (< 3.11) does not have visititems_links,
        # we have to deal with Soft and External links manually.
        for childname in node:
            childnode = node.get(childname, getlink=True)
            if isinstance(childnode, h5py.SoftLink):
                g[name + "/" + childname] = childnode
            elif isinstance(childnode, h5py.ExternalLink):
                g[name + "/" + childname] = childnode


f.visititems(copy_visitor)
# f.visititems_links(copy_visitor) # only in h5py >= 3.11, while DIALS has 3.8

# Copy images

original_data = f["/entry/data/data_000001"]
_, ny, nx = original_data.shape
dcid = original_data.id.get_create_plist()
assert dcid.get_nfilters() == 1
fid, _, fopt, comment = dcid.get_filter(0)

if h5py.h5z.filter_avail(fid):
    h5py.h5z.unregister_filter(fid)

output_data = g.create_dataset(
    "/entry/data/data_000001",
    shape=(len(hits), ny, nx),
    maxshape=(None, ny, nx),
    dtype=original_data.dtype,
    chunks=original_data.chunks,
    compression=fid,
    compression_opts=fopt,
    allow_unknown_filter=True,
)
for idx, hit_id in enumerate(hits):
    source_ds, source_z = frame_to_data[hit_id]
    print(idx + 1, "/", len(hits), "from", hit_id)  # , source_ds, source_z)
    mask, rawdata = source_ds.id.read_direct_chunk((source_z, 0, 0))
    output_data.id.write_direct_chunk((idx, 0, 0), rawdata, mask)

output_data.attrs["image_nr_low"] = 1  # 1 indexed
output_data.attrs["image_nr_high"] = len(hits)

frame_index = g.create_dataset(
    "/entry/original_frame_index", data=hits + 1, compression="gzip", shuffle=True
)  # 1 indexed

g[
    "/entry/description"
] = "This file was created by select_eiger_frames.py (Takanori Nakane, 2024). The frame indices in the original file are stored in /entry/original_frame_index (1-indexed)."

f.close()
g.close()
