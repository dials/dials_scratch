from scitbx.random import variate, normal_distribution
from scitbx.array_family import flex
import math

def centroidify(width, shift, count):
    g = variate(normal_distribution(mean=shift, sigma=width))
    values = flex.double([g.next() for c in range(count)])
    hist = flex.histogram(data=values, n_slots=20,
                          data_min=-10, data_max=10)
    true_mean = flex.sum(values) / values.size()
    true_variance = sum([(v - true_mean) ** 2 for v in values]) / \
      (values.size() - 1)
    total = 1.0 * flex.sum(hist.slots())

    hist_mean = sum([c * v for c, v in zip(hist.slot_centers(),
                                           hist.slots())]) / total

    # equation 6
    hist_var = sum([(v / total) ** 2 * (1.0/12.0) for v in hist.slots()])

    # print input setings
    print '%8.5f %4.1f %4d' % (width ** 2 / count, shift, count),

    # true variance / mean of distribution
    print '%6.3f %8.5f' % (true_mean, true_variance / values.size()),

    # putative values of same derived from histogram
    print '%6.3f %8.5f' % (hist_mean, hist_var)

for width in 0.1, 0.2, 0.5, 1.0, 1.5, 2, 3, 5:
    for shift in 0, 0.1, 0.2, 0.5:
        for count in 10, 20, 50, 100, 200, 500, 1000, 2000, 5000:
            centroidify(width, shift, count)
