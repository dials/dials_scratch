from __future__ import print_function
from cctbx import uctbx

a = 4.15695 # Angstrom
wavelength = 0.68890

unit_cell = uctbx.unit_cell((a, a, a, 90, 90, 90))
print(unit_cell)
from cctbx.array_family import flex

millers = flex.miller_index()
r = 7
for h in range(-r, r+1):
 for k in range(-r, r+1):
  for l in range(-r, r+1):
   millers.append((h, k, l))
resolution = unit_cell.d(millers)

def sort_uniq(seq):
  set = {}
  map(set.__setitem__, seq, [])
  return sorted(set.keys())
resolution = sort_uniq((round(tt, 4) for tt in resolution))
resolution.reverse()

for i, n in enumerate(resolution):
  if n > 0:
    print("Ring %2d: %7.4f A" % (i + 1, n))
