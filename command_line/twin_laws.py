from __future__ import print_function
from iotbx import mtz
import sys

from iotbx.file_reader import any_file

f = any_file(sys.argv[1], force_type="hkl", raise_sorry_if_errors=True)

# find intensity array - if XDS HKL file need to check label
i = [
    ma
    for ma in f.file_content.as_miller_arrays()
    if (ma.is_xray_intensity_array() or "iobs" in ma.info().label_string())
]

if not i:
    raise RuntimeError("no intensities found")

from mmtbx.scaling.twin_analyses import twin_laws

TL = twin_laws(miller_array=i[0])

print("%d possible operators" % len(TL.operators))
for o in TL.operators:
    print("Le-page delta: %.3f operator: %s" % (o.delta_le_page, o.operator))
