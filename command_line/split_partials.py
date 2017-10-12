from dials.array_family import flex
import cPickle as pickle
import sys

refl_in = sys.argv[1]
refl_out = sys.argv[2]
part_lim = float(sys.argv[3])

refl = pickle.load(open(refl_in, 'r'))

# remove duff reflections
sel = refl.get_flags(refl.flags.integrated_sum, all=True)
refl = refl.select(sel)
sel = refl.get_flags(refl.flags.integrated, all=True)
refl = refl.select(sel)

# extract partiality, plot histogram, select subset, dump
part = refl['partiality']
parth = flex.histogram(part, n_slots=20, data_min=0.0, data_max=1.0)
parth.show()
edge = part < part_lim
edge_refl = refl.select(edge)

pickle.dump(edge_refl, open(refl_out, 'w'))

print 'Wrote %d reflections with partiality < %.3f to %s' % \
  (len(edge_refl), part_lim, refl_out)
