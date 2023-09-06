# SPring-8 BL40XU EIGER to NeXus converter

This script adds the detector two theta and the phi axis information by modifying `*_master.h5` from the
EIGER detector at SPring-8 BL40XU and makes it fully NeXus compiant.

## Test data processing

Download [XRDa-159: Cytidine collected at BL40XU, SPring-8](https://xrda.pdbj.org/entry/159).

Working in the same directory:

```
for i in 377 378 379 380;
    do dials.python ~/prog/dials_scratch/tn/SPring-8-BL40XU/sp8bl40xu-to-nxs.py dataset_2_${i}_master.h5 $i.nxs;
done

dials.import *.nxs
dials.find_spots nproc=36 imported.expt

dials.index imported.expt strong.refl space_group=P222
dials.refine indexed.{expt,refl} scan_varying=True

dials.integrate refined.{expt,refl} nproc=16

dials.two_theta_refine integrated.{expt,refl}
dials.scale refined_cell.expt integrated.refl d_min=0.76
dials.export scaled.{expt,refl}
xia2.to_shelx scaled.mtz Cytidine C9N3O5
```

Results:
```
+------------+-------------+---------------+-------------+
|   Imageset |   # indexed |   # unindexed | % indexed   |
|------------+-------------+---------------+-------------|
|          0 |        3066 |            99 | 96.9%       |
|          1 |        3020 |            65 | 97.9%       |
|          2 |        3081 |           133 | 95.9%       |
|          3 |        3113 |           154 | 95.3%       |
+------------+-------------+---------------+-------------+

RMSDs by experiment:
+-------+--------+----------+----------+------------+
|   Exp |   Nref |   RMSD_X |   RMSD_Y |     RMSD_Z |
|    id |        |     (px) |     (px) |   (images) |
|-------+--------+----------+----------+------------|
|     0 |   2577 | 0.12145  |  0.18511 |    0.27112 |
|     1 |   2613 | 0.090719 |  0.17563 |    0.2643  |
|     2 |   2736 | 0.17221  |  0.15853 |    0.26989 |
|     3 |   2594 | 0.1453   |  0.11001 |    0.27678 |
+-------+--------+----------+----------+------------+

  d_max  d_min   #obs  #uniq  mult.   %comp     <I>  <I/sI>  r_mrg  r_meas  r_pim  r_anom  cc1/2   cc_ano
  14.71   2.06   1006     99  10.16  100.00  9232.7    21.8  0.106   0.112  0.034   0.055  0.991*   0.236
   2.06   1.64    856     79  10.84  100.00  5522.4    21.5  0.127   0.134  0.040   0.071  0.994*   0.244
   1.64   1.43    913     83  11.00  100.00  2543.6    17.8  0.135   0.142  0.042   0.065  0.986*  -0.051
   1.43   1.30    744     69  10.78  100.00  2760.3    16.8  0.147   0.156  0.050   0.075  0.987*   0.254
   1.30   1.21    807     80  10.09  100.00  2495.1    14.7  0.150   0.159  0.051   0.122  0.993*   0.251
   1.21   1.14    730     72  10.14  100.00  2344.5    14.2  0.164   0.173  0.054   0.105  0.990*  -0.277
   1.14   1.08    744     76   9.79  100.00  2780.2    14.1  0.157   0.166  0.053   0.071  0.990*  -0.399
   1.08   1.03    636     67   9.49  100.00  2080.2    12.2  0.168   0.179  0.058   0.089  0.980*  -0.173
   1.03   0.99    606     71   8.54  100.00  2312.8    11.1  0.165   0.178  0.063   0.109  0.966*   0.202
   0.99   0.96    723     83   8.71  100.00  1621.4     9.9  0.206   0.219  0.075   0.133  0.975*  -0.475
   0.96   0.93    639     69   9.26  100.00  1135.6     8.9  0.211   0.224  0.073   0.121  0.973*  -0.120
   0.93   0.90    602     73   8.25  100.00  1164.9     7.3  0.205   0.219  0.075   0.141  0.991*  -0.154
   0.90   0.88    539     65   8.29  100.00  1367.3     8.1  0.214   0.228  0.079   0.141  0.977*   0.136
   0.88   0.86    516     67   7.70  100.00  1004.1     6.1  0.233   0.251  0.091   0.186  0.970*   0.160
   0.86   0.84    580     76   7.63  100.00   901.1     5.7  0.246   0.266  0.098   0.174  0.936*  -0.170
   0.84   0.82    522     69   7.57  100.00  1123.0     6.7  0.262   0.282  0.103   0.168  0.936*   0.077
   0.82   0.80    569     83   6.86  100.00  1136.3     5.2  0.231   0.250  0.093   0.179  0.991*  -0.290
   0.80   0.79    357     64   5.58  100.00   969.9     5.0  0.229   0.254  0.107   0.258  0.938*   0.105
   0.79   0.77    188     66   2.85   97.06   832.2     3.5  0.216   0.260  0.141   0.298  0.796*   0.134
   0.77   0.76    115     54   2.13   68.35   918.0     3.5  0.153   0.198  0.125   0.286  0.925*  -0.514
  14.70   0.76  12392   1465   8.46   98.19  2386.5    11.2  0.149   0.159  0.052   0.106  0.989*   0.021
```

The structure can be trivially phased.

However, anisotropic ADP refinement gave non-positive definite (NPD) ADP tensors.
`dials.scale` might have over-sharpened intensities.
Blurring the dataset by B=1.5 helped.

This didn't happen on another (unpublished) dataset, so I guess this is specific to the dataset and
has nothing to do with my converter.

```
servalcat util blur --hklin scaled.mtz -B 1.5
iotbx.reflection_file_converter scaled_blur_1.50.mtz --write_unmerged --shelx=Cytidine_blur_1.50.hkl --label=I,SIGI
# If this crashes saying "assert len(result) == 8", it is too old;
# use the latest version from DIALS.

cp Cytidine.ins Cytidine_blur_1.50.ins
```

## Precision of the goniometer

According to the beam line scientists, the phi axis should be along
(0, -0.70711, 0.70711), i.e. 45 degrees, diagonal.
With this, however, joint indexing in P1 gave skewed angles (90.013(2), 90.1043(19), 90.028(2))
and joint indexing with orthorhombic constraints failed.

So I indexed individual sweeps after setting the phi angles to 0 and compared pairwise
relative orientations with `dials.compare_orientation_matrices`.

This gave (ignore signs):
```
Phi 0 to 90
 (-0.010954586402166367, 0.7092618274876992, -0.7048600266049756) -89.82605289766043
Phi 0 to 180
 (-0.010638440727831471, 0.7093444849680379, -0.704781686215178) -179.85400928198467
Phi 0 to 270
 (-0.010319002254552758, 0.7093494675360645, -0.7047814207956764) 90.24315956411125
Phi 90 to 180
 (-0.010751663136272524, 0.7095661389301968, -0.7045568083721119)  -90.02796284262233
Phi 90 to 270
 (-0.01069443573247834, 0.7095295089450782, -0.7045945678051458) -179.93079951525377
Phi 180 to 270
 (-0.010641742890083093, 0.7095702550980775, -0.7045543317504411) -89.90283701920494
```

Another (unpublished) dataset gave:
```
Phi 0 to 90
 (-0.007710680446344071, 0.7095287730344186, -0.7046342779366661) -90.08314895219911
Phi 0 to 180
 (-0.00908076797878293, 0.7097567763635739, -0.7043882864293691) 179.9283295410904
Phi 0 to 270
 (-0.011141941210044191, 0.7083082581964855, -0.705815321820611) 90.05455263229369
Phi 90 to 180
 (-0.009415374672027079, 0.7087878219796541, -0.7053588973941729) -89.9886355985751
Phi 90 to 270
 (-0.008575861241763653, 0.7077128780865153, -0.706448113306606) -179.86271884322068
Phi 180 to 270
 (-0.007043068295465361, 0.7083135229960625, -0.705862839601216) -89.87425798100747
```

From this I concluded that the phi axis is actually along ~ (0.0106, -0.7094, 0.7047).
With this angles, joint indexing in P1 gave angles closer to 90 degrees (90.0165(7), 90.0112(6), 89.9736(6)).

Note that to perform the above analysis, I hacked the relevant code;
otherwise the command applies rotations allowed in the Laue class and
prints the smallest rotation.

```diff
diff --git a/src/dials/algorithms/indexing/compare_orientation_matrices.py b/src/dials/algorithms/indexing/compare_orientation_matrices.py
index 2ea76ef38..e995b47c2 100644
--- a/src/dials/algorithms/indexing/compare_orientation_matrices.py
+++ b/src/dials/algorithms/indexing/compare_orientation_matrices.py
@@ -33,6 +33,7 @@ def difference_rotation_matrix_axis_angle(crystal_a, crystal_b, target_angle=0):
         axis = axis_angle.axis
         angle = axis_angle.angle() * 180.0 / math.pi
         for sign in (+1, -1):
+            print(op, axis, angle)
             if abs(sign * angle - target_angle) < abs(best_angle - target_angle):
                 best_angle = sign * angle
                 best_axis = tuple(sign * a for a in axis)
```

## Remaining issues

- [ ] NeXus validation
- [ ] Check the absolute hand (should be OK but hasn't been tested)

### NeXus validation

I tried validation using `punx`.

```
pip install pyRestTable pyQt5 punx # I don't know why dependencies are not automatically installed.
punx install
# This installed
# main          user   2023-06-26 08:57:16 d669ffb /home1/XXXX/.config/punx/main
punx validate 127.nxs
```

I don't know why it complains:
```
/entry@None ERROR NXDL field NXmx:title not found 
/entry@None ERROR NXDL field NXmx:end_time not found 
```

The item `/entry/title` has `minOccur="0"` in [NXmx.nxdl.xml](https://github.com/nexusformat/definitions/blob/d669ffb4/applications/NXmx.nxdl.xml#L89),
so it should be optional. The same applies to `end_time`.

I also don't know whether `countrate_correction_applied` is necessary or not.
[NXmx.nxdl.xml](https://github.com/nexusformat/definitions/blob/d669ffb4/applications/NXmx.nxdl.xml#628)
says optional, but [the document](https://manual.nexusformat.org/classes/applications/NXmx.html) says required.
