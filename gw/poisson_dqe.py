def poisson_source(howmany, counts):
  from scitbx.random import variate, poisson_distribution
  g = variate(poisson_distribution(mean=counts))
  return [g.next() for j in range(howmany)]

def measurement_process(counts, dqe):
  from scitbx.random import variate, uniform_distribution
  g = variate(uniform_distribution(min=0.0, max=1.0))
  result = 0
  for j in range(counts):
    if g.next() < dqe:
      result += 1
  return result

def meanvar(values):
  mean = sum(values) / len(values)
  var = sum([(v - mean) ** 2 for v in values]) / (len(values) - 1)
  return mean, var

#cycles = 10000
#for counts in 1, 10, 100, 1000:
#  for dqe in 0.1, 0.5, 0.9:
#    values = poisson_source(cycles, counts)
#    vmean, vvar = meanvar(values)
#
#    measurements = [measurement_process(v, dqe) for v in values]
#    mmean, mvar = meanvar(measurements)
#    print '%5d %.1f %7.2f %7.2f %7.2f %7.2f' % \
#        (counts, dqe, vmean, vvar, mmean, mvar)

# Take one case: A single pixel with a DQE of 0.5 is exposed 1000 times by a
# stable source with a mean count rate of 100 counts per exposure.
dqe = 0.5
cycles = 1000
counts = 100
exposures = poisson_source(cycles, counts)

# The sample mean and variance should be close to 100
print "Original signal mean: {0:.3f}, var: {1:.3f}".format(*meanvar(exposures))

# The measurement process gives another Poisson distribution, now with mean and
# variance close to 50
measurements = [measurement_process(e, dqe) for e in exposures]
mmean, mvar = meanvar(measurements)
print "Measurements mean: {0:.3f}, var: {1:.3f}".format(mmean, mvar)

# The mean of the reconstructed signal is 1/dqe * mmean, but the variance is
# said to be (1/dqe^2) * mvar
print "Theoretical reconstructed signal mean: {0:.3f}, var: {1:.3f}".format(
  1./dqe * mmean, (1./dqe**2) * mvar)

# This is easy to check: just reconstruct the signal from the measurements
# and look at the variance of the reconstructed signal
inv_dqe = 1./dqe
reconstructed = [inv_dqe * m for m in measurements]
rmean, rvar = meanvar(reconstructed)
print "Observed signal mean: {0:.3f}, var: {1:.3f}".format(rmean, rvar)

# The theoretical and observed signal variances are identical.
from libtbx.test_utils import approx_equal
assert approx_equal((1./dqe**2) * mvar, rvar)
print "OK"
