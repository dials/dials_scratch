from __future__ import division
from __future__ import print_function


from numpy.random import poisson

# Compute mean and variance
def mean_and_variance(X):
    N = len(X)
    sum_X = sum(X)
    sum_X2 = sum(x * x for x in X)
    mean = sum_X / N
    var = (sum_X2 - sum_X * sum_X / N) / (N - 1)
    return mean, var


mu = 10
N = 50
ntrials = 10000

# Compute a list of index of dispersions
D_list = []
for i in range(ntrials):
    X = poisson(mu, N)
    mean, var = mean_and_variance(X)
    D = var / mean
    D_list.append(D)

# Compute mean and variance of D(N-1)
DN1_list = [d * (N - 1) for d in D_list]
mean, var = mean_and_variance(DN1_list)

print("Mean Chi(N-1) = ", N - 1)
print("Var Chi(N-1) = ", 2 * (N - 1))
print("Mean D(N-1) = ", mean)
print("Var D(N-1) = ", var)

mean, var = mean_and_variance(D_list)

print("Mean = ", 1)
print("Var = ", 2 / (N - 1))
print("Mean D(N-1) = ", mean)
print("Var D(N-1) = ", var)
