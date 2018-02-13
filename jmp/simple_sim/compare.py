import sys
import cPickle as pickle

im1 = pickle.load(open(sys.argv[1]))
im2 = pickle.load(open(sys.argv[2]))

diff = im1 - im2

print max(diff)
print min(diff)

from matplotlib import pylab
pylab.imshow(diff.as_numpy_array())
pylab.show()
