# Processing SPring-8 BL02B2 datasets

## Problem description

- [DIALS fails to import cbf from Pilatus3 1M at BL02B1, SPRING-8](https://github.com/dials/dials/issues/2488)
- [Crystal structure of a complex of Zinc(II) tetraphenylporphyrin with pyrimidine (XRDa-155)](https://xrda.pdbj.org/entry/155)

The test data in XRDa apparently didn't undergo the swapping of the fast and slow axes unlike the dataset @ndevenish got.

## How to solve the structure

### Sweep 1

```
cd WORKDIR
parallel -P32  dials.python ~/prog/dials_scratch/tn/SPring-8-BL02B1/sp8bl02b1-cbf2smv.py {} {/.}.img ::: ~/data/XRDa-155-BL02B1-ZnTPP/ZNTPP_full_01*.cbf

dials.import ZNTPP_full_01*.img
dials.find_spots nproc=24 imported.expt 
dials.search_beam_position imported.expt strong.refl

dials.index optimised.expt strong.refl
dials.refine indexed.{expt,refl} scan_varying=True

dials.integrate nproc=16 refined.{expt,refl}
dials.two_theta_refine integrated.{expt,refl}
dials.scale refined_cell.expt integrated.refl 

dials.export scaled.{expt,refl}
xia2.to_shelx scaled.mtz ZnTPP ZnC44N4
# trivially solved in Olex2
```

Merging Statistics:

```
d_max  d_min   #obs  #uniq  mult.  %comp     <I>  <I/sI>  r_mrg  r_meas  r_pim  r_anom  cc1/2   cc_ano
10.70   1.14   2431   1372   1.77  95.15  1015.8    36.5  0.065   0.092  0.065   0.137  0.952*   0.000
 1.14   0.90   2700   1414   1.91  97.58   442.2    32.2  0.050   0.071  0.050   0.104  0.987*   0.000
 0.90   0.79   2580   1365   1.89  95.99   223.0    26.1  0.055   0.078  0.055   0.113  0.985*   0.000
 0.79   0.72   2197   1295   1.70  89.25   161.9    21.2  0.060   0.085  0.060   0.124  0.983*   0.000
 0.72   0.67   2437   1340   1.82  92.99   144.6    19.7  0.057   0.081  0.057   0.121  0.986*   0.000
 0.67   0.63   2496   1336   1.87  94.55   105.9    16.6  0.063   0.088  0.063   0.129  0.986*   0.000
 0.63   0.59   2621   1382   1.90  94.21    92.2    15.1  0.068   0.096  0.068   0.142  0.982*   0.000
 0.59   0.57   2604   1364   1.91  94.59    75.9    13.4  0.074   0.105  0.074   0.153  0.980*   0.000
 0.57   0.55   2620   1367   1.92  94.73    56.4    11.0  0.082   0.116  0.082   0.173  0.978*   0.000
 0.55   0.53   2537   1318   1.92  92.23    44.2     9.3  0.095   0.134  0.094   0.187  0.973*   0.000
 0.53   0.51   2206   1229   1.79  85.23    37.8     7.9  0.101   0.143  0.101   0.207  0.966*   0.000
 0.51   0.50   1575    875   1.80  61.15    33.6     7.3  0.105   0.148  0.104   0.221  0.966*   0.000
 0.50   0.48   1274    704   1.81  47.57    30.9     6.9  0.103   0.145  0.103   0.220  0.969*   0.000
 0.48   0.47    903    492   1.84  35.76    25.1     5.9  0.116   0.164  0.116   0.260  0.959*   0.000
 0.47   0.46    861    469   1.84  31.02    22.4     5.5  0.121   0.171  0.120   0.256  0.955*   0.000
 0.46   0.45    556    317   1.75  22.56    20.1     4.5  0.156   0.221  0.156   0.339  0.952*   0.000
 0.45   0.44    460    266   1.73  18.77    17.4     4.3  0.145   0.204  0.145   0.269  0.953*   0.000
 0.44   0.43    290    170   1.71  11.63    15.8     3.7  0.155   0.219  0.155   0.321  0.941*   0.000
 0.43   0.43    161    112   1.44   7.81    12.7     2.9  0.195   0.274  0.193   0.489  0.910*   0.000
 0.43   0.42     60     51   1.18   3.59    10.1     2.1  0.161   0.227  0.161   0.300  0.947*   0.000
10.69   0.42  33569  18238   1.84  63.30   184.8    16.7  0.064   0.090  0.063   0.131  0.970*   0.000
```

## Sweep 4

This worked in the same way as the sweep 1.

## Multi-sweep indexing

The CHI axis is along (0, 0, 1), i.e., parallel to the beam towards the source, not (0, -1, 0)!!

After fixing this, multi-sweep indexing worked.

## Issues

- [ ] This is PAD, not CCD.
- [ ] Test sweeps with non-zero two theta angles (XRDa-155 does not have such sweeps, though).
- [X] Multi-sweep indexing.
- [ ] The beam center is off; is the header wrong or is my interpretation wrong?
- [ ] I don't know if the hand is correct; the test data is P-1 so I cannot check.
- [ ] Ideally SPring-8 people should use full CBF. If they don't, our converter should write full CBF.
