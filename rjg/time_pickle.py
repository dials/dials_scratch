from __future__ import print_function


def run(file_name):
    import time
    from libtbx import easy_pickle
    from dials.array_family import flex

    t0 = time.time()
    refl = flex.reflection_table.from_pickle(file_name)
    t1 = time.time()
    print("Time reflection_table.from_pickle(): %.3f" % (t1 - t0))

    refl.as_pickle("tmp.pickle")
    t2 = time.time()
    print("Time reflection_table.as_pickle(): %.3f" % (t2 - t1))

    d = dict(((k, refl[k]) for k in refl.keys()))
    t3 = time.time()
    easy_pickle.dump("tmp.pickle", d)
    t4 = time.time()
    print("Time pickle dict: %.3f" % (t4 - t3))

    for k, v in d.iteritems():
        t0 = time.time()
        easy_pickle.dump("tmp.pickle", v)
        t1 = time.time()
        print("Column %s (%s): %.3f" % (k, type(v), t1 - t0))


if __name__ == "__main__":
    import sys

    args = sys.argv[1:]
    assert len(args) == 1
    run(args[0])
