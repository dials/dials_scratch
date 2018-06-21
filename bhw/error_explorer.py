#!/usr/bin/env python
#-*- coding:utf-8 -*-

"""Learning a little about the distribution of diffraction spot intensities
and their errors, reading from an unmerged MTZ file."""

import sys
import numpy as np
from scipy import stats as ss
from iotbx import mtz
from cctbx.array_family import flex
from matplotlib import colors, pyplot as plt


class DataDist:
  def __init__(self, filename, keep_singles=False, propagate_errors=False):
    self.filename = filename
    self.data_from_unmerged_mtz(filename)
    self.select_by_multiplicity(keep_singles)
    self.kept_singles = keep_singles
    self.mean_error_stddev(propagate_errors)
    self.make_z()
    

  def data_from_unmerged_mtz(self, filename):
    m = mtz.object(filename)  #Parse MTZ, with lots of useful methods.
    
    self.ind = m.extract_miller_indices()  #A flex array of Miller indices.
    
    cols = m.columns()  #Generates columns (augmented flex arrays).
    col_dict = { c.label() : c for c in cols }  #A dict of all the columns.
    self.I, self.sigI, self.x, self.y = (
      col_dict[label].extract_values().as_double()
      for label in ('I', 'SIGI', 'XDET', 'YDET')
    )
  
  def select_by_multiplicity(self, keep_singles=False):    
    # Find unique Miller indices.
    ind_unique = flex.miller_index(np.unique(self.ind, axis=0))

    # Record the multiplicities.
    multis = flex.int([])
    
    for hkl in ind_unique:
      sel = (self.ind == hkl).iselection()
      multi = sel.size()
      multis.extend(multi * flex.int(np.ones(multi).astype(np.uint32)))
    
    if keep_singles:
      self.multis = multis
      self.ind_unique = ind_unique
    else:
      sel = (multis != 1).iselection()
      self.multis = multis.select(sel)
      self.ind, self.I, self.sigI, self.x, self.y = map(
        lambda x: x.select(sel),
        (self.ind, self.I, self.sigI, self.x, self.y)
      )
      self.ind_unique = flex.miller_index(np.unique(self.ind, axis=0))
  
  def mean_error_stddev(self, propagate_errors=False):
    # Calculate the weighted mean intensities.
    self.Imeans = flex.double([])
    
    if propagate_errors:
      # Calculate the standard errors on the means.
      self.sigImeans = flex.double([])
      
      # Calculate the standard deviations from unbiased weighted variances.
      self.stddevs = flex.double([])
    
    # For weighted averaging.
    weights = 1 / flex.pow2(self.sigI)
    
    for hkl in self.ind_unique:
      sel = (self.ind == hkl).iselection()
      multi = sel.size()
      ones = flex.double(np.ones(multi).astype(np.float64))
      weight = weights.select(sel)
      sum_weight = flex.sum(weight)
      selI = self.I.select(sel)
      
      Imean = flex.sum( weight * selI ) / sum_weight
      self.Imeans.extend( Imean * ones )
      
      if propagate_errors:
        sigImean = flex.sqrt(1 / sum_weight)
        self.sigImeans.extend(sigImean)
        
        if multi == 1:
          self.stddevs.extend(self.sigI.select(sel))
        else:
          variance = (
            flex.sum(weight * flex.pow2(selI - Imean)) /
            (sum_weight - flex.sum(flex.pow2(weight)) / sum_weight)
          )
          stddev = flex.sqrt(variance * ones)
          self.stddevs.extend(stddev)
  
  def make_z(self, error='sigma'):
    if error in ('stddev', 'sigImean'):
      error = {'stddev':self.stddevs, 'sigImean':self.sigImeans}[error]
    else:
      error = self.sigI
    
    self.z = (self.I - self.Imeans) / error
    self.order = flex.sort_permutation(self.z)
  
  def plot_time_series(self, overlay_mean=False, overlay_error=False):
    for hkl in self.ind_unique:
      sel = (self.ind == hkl).iselection()
      
      if overlay_error == 'sigImean':
        yerr = self.sigImeans.select(sel)
      elif overlay_error:
        yerr = self.stddevs.select(sel)
      else:
        yerr = None
      
      plt.errorbar(
        sel,
        self.I.select(sel),
        yerr=self.sigI.select(sel),
        ls="--"
      )
      if overlay_mean:
        plt.errorbar(
          sel,
          self.Imeans.select(sel),
          yerr = yerr,
          ls = "-",
          color = "k",
          lw = .5
        )

  def probplot(self):
    self.osm, self.osr = ss.probplot(self.z, fit=False)
    
    fig, ax = plt.subplots()
    
    ax.set_title('Normal probability plot')
    ax.set_xlabel('Order statistic medians, $m$')
    ax.set_ylabel(
      r'Ordered responses, $z$'
    )
    ax.plot(self.osm, self.osr, '.b')
    ax.plot([-5,5], [-5,5], '-g')
    
    fig.savefig(self.filename.split('.')[0] + '_probplot')
    plt.close()
  
  def deviation_vs_multiplicity(self):
    if not (self.osm.any() and self.osr.any()):
      self.probplot()
    
    fig, ax = plt.subplots()
    
    ax.set_title(
      r'Difference between ordered responses, $z$, '
      + r'and order statistic medians, $m$'
    )
    ax.set_xlabel('Multiplicity')
    ax.set_ylabel(r'$z - m$')
    ax.plot(self.multis.select(self.order), self.osr - self.osm, '.')
    
    fig.savefig(self.filename.split('.')[0] + '_deviation_vs_multiplicity')
    plt.close()
  
  def deviation_map(self):
    if not (self.osm.any() and self.osr.any()):
      self.probplot()
    
    extreme = np.ceil(np.abs(self.osr - self.osm).max())
    norm = colors.SymLogNorm(
      vmin = -extreme, vmax = extreme,
      linthresh = .02, linscale=1,
    )
    cmap_kws = {'cmap' : 'coolwarm_r', 'norm' : norm}
    
    fig, ax = plt.subplots()
    
    ax.set_title(
      r'Difference between ordered responses, $z$, '
      + r'and order statistic medians, $m$'
    )
    ax.set_xlabel('Detector x position (pixels)')
    ax.set_ylabel('Detector y position (pixels)')
    ax.set_aspect('equal', 'box')
    det_map = ax.scatter(
      self.x.select(self.order),
      self.y.select(self.order),
      c = self.osr - self.osm,
      marker = ',',
      s = 0.5,
      **cmap_kws
    )
    cbar = fig.colorbar(det_map, ax=ax, **cmap_kws)
    cbar.set_label(r'$z - m$')
    
    fig.savefig(self.filename.split('.')[0] + '_deviation_detector_map')
    plt.close()

  # TODO Add time series plot (v easy)
  
if __name__ == "__main__":
  # TODO handle multiple input MTZ files.
  data = DataDist(sys.argv[1]) #Give an unmerged MTZ file as an argument.
  
  data.probplot()
  data.deviation_vs_multiplicity()
  data.deviation_map()
