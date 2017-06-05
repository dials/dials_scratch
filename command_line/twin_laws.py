from iotbx import mtz
import sys

m = mtz.object(sys.argv[1])
mad = m.as_miller_arrays_dict()
k = [k for k in mad.keys() if k[2] == 'I'][0]
i = mad[k]

from mmtbx.scaling.twin_analyses import twin_laws
TL = twin_laws(miller_array=i)

print '%d possible operators' % len(TL.operators)
for o in TL.operators:
    print o.delta_le_page
    print o.operator
