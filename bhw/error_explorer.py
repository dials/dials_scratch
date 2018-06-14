#!/usr/bin/env python
#-*- coding:utf-8 -*-

"""Learning a little about the distribution of spot intensities and their
errors."""

import sys
import numpy as np
from scipy import stats as ss
from iotbx import mtz
from cctbx.array_family import flex
from matplotlib import colors, pyplot as plt


def data_from_unmerged_mtz(filename):
  m = mtz.object(filename)  #Parse MTZ, with lots of useful methods.
  ind = m.extract_miller_indices()  #A flex array of Miller indices.
  cols = m.columns()  #Generates columns (augmented flex arrays).
  col_dict = { c.label() : c for c in cols }  #A dict of all the columns.
  
  I, sigI, x, y = (
    col_dict[label].extract_values().as_double()
    for label in ('I', 'SIGI', 'XDET', 'YDET')
  )
  
  return ind, I, sigI, x, y

def mean_error_stddev(ind, I, sigI):
  #Record the multiplicities.
  multis = flex.int([])
  
  #For weighted averaging.
  weights = 1 / flex.pow2(sigI)
  
  #Find unique Miller indices.
  ind_unique = flex.miller_index(np.unique(ind, axis=0))
  
  #Calculate the weighted mean intensities and the standard errors on them.
  sImeans = flex.double([])
  Imeans = flex.double([])
  
  #Calculate the standard deviations from unbiased weighted variances.
  stddevs = flex.double([])
  
  for hkl in ind_unique:
    sel = (ind == hkl).iselection()
    multi = sel.size()
    multis.extend( multi * flex.int(np.ones(multi).astype(np.uint32)) )
    
    selI = I.select(sel)
    
    #Deal with multiplicity == 1 cases, for which variance is undefined.
    if multi == 1:
      Imeans.extend(selI)
      stddevs.extend(sigI.select(sel))
    
    else:
      ones = flex.double(np.ones(multi).astype(np.float64))
      weight = weights.select(sel)
      sum_weight = flex.sum(weight)
      
      Imean = flex.sum( weight * selI ) / sum_weight
      Imeans.extend( Imean * ones )
      
      variance = (
        flex.sum( weight * flex.pow2( selI - Imean ) ) /
        ( sum_weight - flex.sum(flex.pow2(weight)) / sum_weight )
      )
      stddev = flex.sqrt(variance * ones)
      stddevs.extend(stddev)
  
  return ind_unique, multis, Imeans, sImeans, stddevs

def plot_time_series(ind, ind_unique, I, sigI, overlay=False):
  for hkl in ind_unique:
    sel = (ind == hkl).iselection()
    plt.errorbar(
      sel,
      I.select(sel),
      yerr=sigI.select(sel),
      ls="--"
    )
    if overlay:
      plt.errorbar(
        sel,
        overlay[0].select(sel),
        yerr = overlay[1].select(sel),
        ls = "-",
        color = "k",
        lw = .5
      )

def probplots(I, Imeans, stddevs, multis, x, y):
  scaled = ( (I - Imeans) / stddevs )
  order = flex.sort_permutation(scaled)
  osm, osr = ss.probplot(scaled, fit=False)
  
  fig = plt.figure()
  
  extreme = np.ceil(np.abs(osr - osm).max())
  norm = colors.SymLogNorm(
    vmin = -extreme, vmax = extreme,
    linthresh = .02, linscale=1,
  )
  cmap_kws = {'cmap' : 'coolwarm_r', 'norm' : norm}
  
  ax = fig.add_subplot(221)
  ax.set_title('Normal probability plot')
  ax.set_xlabel('Order statistic medians, $m$')
  ax.set_ylabel(
    'Ordered responses, ' +
    r'$\frac{I_{\mathbf{h}, i} - \bar{I}_\mathbf{h}}{s_\mathbf{h}}$'
  )
  ax.plot(osm, osr, '.b')
  ax.plot([-5,5], [-5,5], '-g')
  
  ax = fig.add_subplot(223)
  ax.set_title(
    'Difference between ordered responses and order statistic medians'
  )
  ax.set_xlabel('Multiplicity')
  ax.set_ylabel(
    r'$\frac{I_{\mathbf{h}, i} - \bar{I}_\mathbf{h}}{s_\mathbf{h}} - m$'
  )
  ax.plot(multis.select(order), osr - osm, '.')

  ax = fig.add_subplot(122)
  ax.set_title(
    'Difference between ordered responses and order statistic medians'
  )
  ax.set_xlabel('Detector x position (pixels)')
  ax.set_ylabel('Detector y position (pixels)')
  ax.set_aspect('equal', 'box')
  det_map = ax.scatter(
    x.select(order),
    y.select(order),
    c = osr - osm,
    marker = ',',
    s = 0.5,
    **cmap_kws
  )
  cbar = fig.colorbar(det_map, ax=ax, **cmap_kws)
  
  return osm, osr, order


if __name__ == "__main__":
  ind, I, sigI, x, y = data_from_unmerged_mtz(
    sys.argv[1] #Give an unmerged MTZ file as an argument.
  )
  
  #plot_time_series(ind, ind_unique, I, sigI)
  #plt.show()
  
  ind_unique, multis, Imeans, sImeans, stddevs = mean_error_stddev(
    ind, I, sigI
  )
  
  #plot_time_series(ind, ind_unique, I, sigI, overlay=(Imeans, stddevs))
  #plt.show()
  
  sel = (multis > 1).iselection()
  I, Imeans, stddevs, multis, x, y = map(
    lambda x: x.select(sel),
    (I, Imeans, stddevs, multis, x, y)
  )
  
  osm, osr, order = probplots(I, Imeans, stddevs, multis, x, y)
  I, Imeans, stddevs, multis, x, y = map(
    lambda x: x.select(order),
    (I, Imeans, stddevs, multis, x, y)
  )
  
  plt.show()
