from __future__ import division, print_function


def shelxd_cc_hist(filename):
    """Read the logs from filename (from shelxd) scrape out CC's, compute
    histogram of all, weak, write to stdout"""
    from cctbx.array_family import flex

    all = flex.double()
    weak = flex.double()

    for record in open(filename):
        if not record.startswith(" Try"):
            continue
        a = float(record[31:36])
        w = float(record[38:43])
        all.append(a)
        weak.append(w)

    h_all = flex.histogram(all, n_slots=220, data_min=-10, data_max=100)
    h_weak = flex.histogram(weak, n_slots=220, data_min=-10, data_max=100)

    for b, a, w in zip(h_all.slot_centers(), h_all.slots(), h_weak.slots()):
        print(b, a, w)


if __name__ == "__main__":
    import sys

    shelxd_cc_hist(sys.argv[1])
