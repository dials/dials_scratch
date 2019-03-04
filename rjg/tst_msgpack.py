from __future__ import print_function


def run(args):
    assert len(args)
    f = args[0]
    from libtbx import easy_pickle

    refl = easy_pickle.load(f)
    xyzobs_px = refl["xyzobs.px.value"]

    d = {}
    for k in refl.keys():
        if k == "shoebox":
            continue
        d[k] = list(refl[k])

    import msgpack

    with open("reflection.mpk", "wb") as f:
        print(msgpack.packb(d), file=f)


if __name__ == "__main__":
    import sys

    run(sys.argv[1:])
