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

cycles = 10000
for counts in 1, 10, 100, 1000:
  for dqe in 0.1, 0.5, 0.9:
    values = poisson_source(cycles, counts)
    vmean, vvar = meanvar(values)

    measurements = [measurement_process(v, dqe) for v in values]
    mmean, mvar = meanvar(measurements)
    print '%5d %.1f %7.2f %7.2f %7.2f %7.2f' % \
        (counts, dqe, vmean, vvar, mmean, mvar)

