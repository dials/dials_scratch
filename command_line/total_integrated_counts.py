from __future__ import print_function
from dials.array_family import flex
import cPickle as pickle


def total_count(filtered_reflections):
    return (
        flex.sum(filtered_reflections["intensity.sum.value"]),
        filtered_reflections["intensity.sum.value"].size(),
    )


def filter_reflections(reflections):
    selection = reflections.get_flags(reflections.flags.integrated_sum)
    reflections = reflections.select(selection)
    selection = reflections["intensity.sum.variance"] > 0
    reflections = reflections.select(selection)
    return reflections


def main():
    import sys

    data = pickle.load(open(sys.argv[1]))
    print("%e %d" % total_count(filter_reflections(data)))


if __name__ == "__main__":
    main()
