#!/usr/bin/env dials.python

from argparse import ArgumentParser
from typing import Set

import h5py


def find_external_links(root_group: h5py.Group) -> Set[str]:
    """Given an HDF5 Group, returns the filenames of all external links."""

    found_links = set()

    def _visit(name):
        # visit() will only return hardlink Datasets and Groups, non-recursively
        if root_group.get(name, getclass=True) is h5py.Dataset:
            return
        # We need to manually walk this group to check for external links
        group = root_group[name]
        for key in group:
            link = group.get(key, getlink=True)
            if isinstance(link, h5py.ExternalLink):
                found_links.add(link.filename)

    root_group.visit(_visit)
    return found_links


if __name__ == "__main__":
    parser = ArgumentParser(
        description="Reads an H5 file to extract all external links"
    )
    parser.add_argument("file", help="The h5 file to read", metavar="FILE")
    args = parser.parse_args()

    with h5py.File(args.file, "r") as f:
        print("\n".join(sorted(find_external_links(f))))
