from __future__ import print_function

from math import ceil, floor

for j in range(1, 20):
    for i in range(1, j + 1):
        d = int(floor(float(j) / float(i)))
        r = j % i

        jobs = [d] * i
        k = 0
        while r > 0:
            jobs[k] += 1
            r -= 1
            k += 1

        print("#blocks = %d, #jobs = %d => %s" % (j, i, jobs), sum(jobs) == j)
        # print "#blocks = %d, #jobs = %d, B/J = %d, R = %d" % (j, i, d, r)
