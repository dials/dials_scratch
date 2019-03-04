from __future__ import print_function, division

from iotbx import mtz
from scitbx.array_family import flex

merged = mtz.object("AUTOMATIC_DEFAULT_scaled.mtz")
unmerged = mtz.object("AUTOMATIC_DEFAULT_scaled_unmerged.mtz")

unique_indices = merged.extract_miller_indices()
all_indices = unmerged.extract_miller_indices()
imean = merged.extract_observations("IMEAN", "SIGIMEAN")
iuniq = unmerged.extract_observations("I", "SIGI")

z = flex.double(iuniq.data.size(), 0.0)

for index, i in zip(unique_indices, imean.data):
    sel = all_indices == index
    z.set_selected(sel, (iuniq.data - i) / iuniq.sigmas)

h = flex.histogram(z, data_min=-5, data_max=5, n_slots=100)

for c, v in zip(h.slot_centers(), h.slots()):
    print(c, v)
