# Processing SPring-8 BL02B2 datasets

## Problem description

- [DIALS fails to import cbf from Pilatus3 1M at BL02B1, SPRING-8](https://github.com/dials/dials/issues/2488)
- [Crystal structure of a complex of Zinc(II) tetraphenylporphyrin with pyrimidine (XRDa-155)](https://xrda.pdbj.org/entry/155)

The test data in XRDa apparently didn't undergo the swapping of the fast and slow axes unlike the dataset @ndevenish got.

The deposited dataset contains CBF files written by PILATUS and metadata in INF files.
The INF files are similar to the RIGAKU SMV header.

## Conversion to full CBF

I wrote a converter from CBF+INF to full CBF.

```
parallel -P36 dials.python ~/prog/dials_scratch/tn/SPring-8-BL02B1/sp8bl02b1-fullcbf.py {} {/.}.cbf :::  ~/data/XRDa-155-BL02B1-ZnTPP/ZNTPP_full_*.cbf

mkdir process
cd process
dials.import ../ZNTPP_full_*.cbf

dials.find_spots nproc=36 imported.expt

dials.index optimised.expt strong.refl
dials.refine indexed.{expt,refl} scan_varying=True

dials.integrate refined.* nproc=36

dials.two_theta_refine integrated.{expt,refl}
dials.scale refined_cell.expt integrated.refl d_min=0.40

dials.export scaled.{expt,refl}
xia2.to_shelx scaled.mtz ZnTPP ZnC44N4
# Trivially solved and refined in Olex2 and Servalcat
# Bonding electrons were visible :)
```

```
Shift: 0.01, 2.07 mm (0.1, 12.0 px)

+------------+-------------+---------------+-------------+
|   Imageset |   # indexed |   # unindexed | % indexed   |
|------------+-------------+---------------+-------------|
|          0 |       20137 |           137 | 99.3%       |
|          1 |       11350 |            37 | 99.7%       |
|          2 |       10750 |            82 | 99.2%       |
|          3 |       15769 |           133 | 99.2%       |
|          4 |        8894 |            54 | 99.4%       |
|          5 |        5811 |            19 | 99.7%       |
|          6 |        5691 |            75 | 98.7%       |
|          7 |        8397 |           132 | 98.5%       |
+------------+-------------+---------------+-------------+

RMSDs by experiment:
+-------+--------+----------+----------+------------+
|   Exp |   Nref |   RMSD_X |   RMSD_Y |     RMSD_Z |
|    id |        |     (px) |     (px) |   (images) |
|-------+--------+----------+----------+------------|
|     0 |  17890 |  0.165   |  0.16745 |    0.21241 |
|     1 |   9093 |  0.15843 |  0.15011 |    0.22843 |
|     2 |   9045 |  0.14796 |  0.1843  |    0.20834 |
|     3 |  13107 |  0.15457 |  0.21879 |    0.22071 |
|     4 |   7841 |  0.13846 |  0.16304 |    0.21424 |
|     5 |   4760 |  0.14229 |  0.13068 |    0.22278 |
|     6 |   4904 |  0.12599 |  0.15705 |    0.20286 |
|     7 |   7058 |  0.11939 |  0.17651 |    0.22187 |
+-------+--------+----------+----------+------------+

  d_max  d_min    #obs  #uniq  mult.   %comp    <I>  <I/sI>  r_mrg  r_meas  r_pim  r_anom  cc1/2  cc_ano
  10.70   1.09   17814   1662  10.72  100.00  514.1    44.6  0.056   0.058  0.018   0.039  0.999*  0.118*
   1.09   0.86   18773   1666  11.27  100.00  181.2    37.2  0.066   0.070  0.021   0.045  0.998*  0.090*
   0.86   0.75   16055   1669   9.62  100.00   94.5    28.5  0.080   0.085  0.027   0.061  0.996*  0.050
   0.75   0.68   16280   1633   9.97  100.00   73.4    24.4  0.091   0.096  0.030   0.063  0.997*  0.023
   0.68   0.63   17722   1684  10.52  100.00   52.8    21.6  0.111   0.116  0.036   0.073  0.996* -0.071
   0.63   0.60   17557   1642  10.69  100.00   44.0    18.9  0.124   0.131  0.040   0.081  0.996*  0.014
   0.60   0.57   17784   1652  10.77  100.00   34.4    16.4  0.146   0.154  0.047   0.092  0.995*  0.046
   0.57   0.54   17849   1666  10.71  100.00   24.3    13.5  0.184   0.194  0.059   0.108  0.993*  0.010
   0.54   0.52   17427   1673  10.42  100.00   17.7    10.8  0.237   0.250  0.077   0.139  0.989*  0.001
   0.52   0.50   14336   1667   8.60  100.00   15.6     8.7  0.263   0.280  0.095   0.175  0.979* -0.036
   0.50   0.49   11534   1647   7.00   99.94   13.5     6.8  0.296   0.320  0.120   0.229  0.967*  0.003
   0.49   0.47    8731   1610   5.42   98.77   12.6     5.3  0.309   0.343  0.144   0.296  0.942*  0.001
   0.47   0.46    8121   1696   4.79   97.86   10.5     4.4  0.344   0.385  0.170   0.368  0.919*  0.008
   0.46   0.45    6963   1607   4.33   98.65    8.5     3.2  0.437   0.497  0.230   0.477  0.880*  0.015
   0.45   0.44    6497   1669   3.89   99.23    6.9     2.6  0.522   0.604  0.295   0.626  0.794*  0.018
   0.44   0.43    5493   1616   3.40   99.38    6.1     2.1  0.572   0.679  0.356   0.751  0.758* -0.092
   0.43   0.42    5003   1646   3.04   99.58    5.6     1.8  0.638   0.777  0.435   0.934  0.717* -0.003
   0.42   0.41    4552   1659   2.74   99.88    5.3     1.5  0.689   0.860  0.508   1.024  0.642*  0.054
   0.41   0.41    4449   1674   2.66   99.76    4.3     1.2  0.865   1.088  0.653   1.211  0.542* -0.118
   0.41   0.40    4401   1671   2.63   99.70    3.5     1.0  0.966   1.216  0.732   1.327  0.487*  0.045
  10.70   0.40  237341  33109   7.17   99.63   56.6    12.7  0.089   0.095  0.032   0.086  0.999*  0.121*
```

## How to solve the structure (old way)

This is my old approach.

Here I converted CBF+INF files into the RIGAKU Saturn SMV format.
This worked but I don't like this approach because the format is proprietary and not well-defined.

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

The CHI axis is along (0, 0, -1), i.e., parallel to the beam towards the detector, not (0, -1, 0)!!

Sweeps 5 to 7 are from the same crystal but with the detector two theta at -30.
This information is stored in `CCD_DETECTOR_VECTORS` not `CCD_GONIO_VALUES` (why?).

```
parallel -P32 dials.python ~/prog/dials_scratch/tn/SPring-8-BL02B1/sp8bl02b1-cbf2smv.py {} {/.}.img  ::: ~/data/XRDa-155-BL02B1-ZnTPP/ZNTPP_full_0*.cbf
dials.import ZNTPP_full_0*.img
dials.find_spots imported.expt nproc=36
dials.search_beam_position imported.expt strong.refl
dials.index optimised.expt strong.refl
dials.refine indexed.{expt,refl} scan_varying=True
dials.integrate refined.* nproc=36
dials.two_theta_refine integrated.{expt,refl}
dials.scale refined_cell.expt integrated.refl d_min=0.40
```

```
+------------+-------------+---------------+-------------+
|   Imageset |   # indexed |   # unindexed | % indexed   |
|------------+-------------+---------------+-------------|
|          0 |       13954 |            53 | 99.6%       |
|          1 |        6128 |            22 | 99.6%       |
|          2 |        5906 |            31 | 99.5%       |
|          3 |        9938 |            45 | 99.5%       |
|          4 |        5269 |            20 | 99.6%       |
|          5 |        3086 |            13 | 99.6%       |
|          6 |        3001 |            38 | 98.7%       |
|          7 |        4658 |            50 | 98.9%       |
+------------+-------------+---------------+-------------+

d_max  d_min    #obs  #uniq  mult.   %comp    <I>  <I/sI>  r_mrg  r_meas  r_pim  r_anom  cc1/2  cc_ano
10.70   1.09   17789   1662  10.70  100.00  507.5    36.3  0.061   0.064  0.020   0.042  0.998*  0.166*
 1.09   0.86   18790   1668  11.26  100.00  183.6    27.3  0.088   0.092  0.027   0.058  0.997* -0.072
 0.86   0.75   16045   1668   9.62  100.00   94.8    19.5  0.118   0.125  0.040   0.090  0.991*  0.055
 0.75   0.68   16312   1635   9.98  100.00   74.5    16.0  0.145   0.153  0.048   0.106  0.991*  0.003
 0.68   0.63   17706   1683  10.52  100.00   53.7    13.6  0.180   0.189  0.058   0.125  0.987* -0.052
 0.63   0.60   17572   1644  10.69  100.00   45.3    11.5  0.199   0.209  0.064   0.133  0.986*  0.039
 0.60   0.57   17798   1652  10.77  100.00   35.1     9.7  0.235   0.247  0.075   0.155  0.978* -0.155
 0.57   0.54   17866   1666  10.72  100.00   24.9     7.6  0.280   0.294  0.089   0.186  0.974* -0.065
 0.54   0.52   17461   1677  10.41  100.00   18.5     5.9  0.331   0.348  0.107   0.217  0.963*  0.001
 0.52   0.50   14321   1666   8.60  100.00   16.4     4.7  0.343   0.365  0.124   0.247  0.948* -0.091
 0.50   0.49   11535   1651   6.99   99.94   14.4     3.6  0.363   0.393  0.147   0.297  0.934* -0.103
 0.49   0.47    8712   1604   5.43   98.71   13.4     2.8  0.362   0.402  0.169   0.346  0.917*  0.004 
 0.47   0.46    8139   1701   4.78   97.76   11.3     2.3  0.414   0.465  0.207   0.447  0.861* -0.046 
 0.46   0.45    6932   1602   4.33   98.71    9.8     1.6  0.526   0.601  0.283   0.602  0.807* -0.033 
 0.45   0.44    6482   1672   3.88   99.17    8.2     1.3  0.624   0.728  0.362   0.783  0.690* -0.024 
 0.44   0.43    5482   1613   3.40   99.45    7.4     1.0  0.716   0.856  0.457   0.962  0.649* -0.050 
 0.43   0.42    5020   1658   3.03   99.58    6.9     0.9  0.824   1.010  0.573   1.244  0.587* -0.048 
 0.42   0.41    4509   1648   2.74   99.88    6.7     0.7  0.868   1.086  0.646   1.320  0.495* -0.016 
 0.41   0.41    4459   1679   2.66   99.76    5.8     0.6  1.092   1.377  0.829   1.538  0.412* -0.042 
 0.41   0.40    4401   1674   2.63   99.70    4.7     0.5  1.253   1.585  0.961   1.717  0.340* -0.128 
10.70   0.40  237331  33123   7.17   99.63   57.3     8.4  0.124   0.133  0.045   0.124  0.998*  0.071*
```

## Issues

- [X] This is PAD, not CCD. => Solved with the full CBF.
- [X] Sweeps with non-zero two theta angles.
- [X] Multi-sweep indexing.
- [ ] dials.image_viewer shows the origin at upper left when fs=(1,0,0) and ss=(0,1-,0).  
      However, ours have fs=(-1,0,0) and ss=(0,1,0) so the origin is at the lower right. Is this correct?
- [X] The beam center is off => I guess the header is not accurate.
- [ ] I don't know if the hand is correct; the test data is P-1 so I cannot check.
- [X] Ideally SPring-8 people should use full CBF. If they don't, our converter should write full CBF.
