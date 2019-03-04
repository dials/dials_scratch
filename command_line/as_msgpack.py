from __future__ import print_function, division

if __name__ == "__main__":
    from dials.array_family import flex
    import sys

    flex.reflection_table.from_pickle(sys.argv[1]).as_msgpack_file(
        sys.argv[1].replace("pickle", "msgpack")
    )
