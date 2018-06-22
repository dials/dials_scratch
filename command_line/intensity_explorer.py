#!/usr/bin/env python
#-*- coding:utf-8 -*-

"""Learning a little about the distribution of diffraction spot intensities
and their errors, reading from an unmerged MTZ file."""

# TODO Proper documentation.

import os
import sys
import math
import numpy as np
from scipy import stats as ss
from iotbx import mtz
from cctbx.array_family import flex
from matplotlib import colors, pyplot as plt


class DataDist:
  def __init__(self,
      filename,
      outfile=False,
      keep_singles=False,
      propagate_errors=True
  ):
    if outfile:
      self.outfile = outfile
    else:
      self.outfile = filename
    self.data_from_unmerged_mtz(filename)
    (self.multis, self.ind, self.I, self.sigI,
      self.x, self.y, self.image, self.ind_unique,
      self.kept_singles) = self.select_by_multiplicity(keep_singles)
    self.kept_singles = keep_singles
    self.mean_error_stddev(propagate_errors)
    self.make_z()

  def data_from_unmerged_mtz(self, filename):
    m = mtz.object(filename)  #Parse MTZ, with lots of useful methods.

    self.ind = m.extract_miller_indices()  #A flex array of Miller indices.

    cols = m.columns()  #Generates columns (augmented flex arrays).
    col_dict = { c.label() : c for c in cols }  #A dict of all the columns.
    self.I, self.sigI, self.x, self.y, self.image = (
      col_dict[label].extract_values().as_double()
      for label in ('I', 'SIGI', 'XDET', 'YDET', 'BATCH')
    )

  def select_by_multiplicity(self, keep_singles=False):
    # Find unique Miller indices.
    ind_unique = set(self.ind)

    # Record the multiplicities.
    multis = flex.int(self.ind.size(), 0)

    for hkl in ind_unique:
      sel = (self.ind == hkl).iselection()
      multis.set_selection(sel, sel.size())

    # Drop multiplicity-1 data unless instructed otherwise.
    if not keep_singles:
      sel = (multis != 1).iselection()
      multis, ind, I, sigI, x, y, image = map(
        lambda col: col.select(sel),
        (multis, self.ind, self.I, self.sigI, self.x, self.y, self.image)
      )
      ind_unique = set(self.ind)

    return multis, ind, I, sigI, x, y, image, ind_unique, keep_singles

# FIXME Dodgy order of data in Imeans, sigImeans, stddevs.
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
        sigImean = math.sqrt(1 / sum_weight)
        self.sigImeans.extend(sigImean * ones)

        if multi == 1:
          self.stddevs.extend(self.sigI.select(sel))
        else:
          variance = (
            flex.sum(weight * flex.pow2(selI - Imean)) /
            (sum_weight - flex.sum(flex.pow2(weight)) / sum_weight)
          )
          stddev = flex.sqrt(variance * ones)
          self.stddevs.extend(stddev)
    #self.Imeans = self.gw_imeans()

  def make_z(self, error='sigma'):
    if error in ('stddev', 'sigImean'):
      error = {'stddev':self.stddevs, 'sigImean':self.sigImeans}[error]
    else:
      error = self.sigI

    self.z = (self.I - self.Imeans) / error

    self.order = flex.sort_permutation(self.z)

  def gw_imeans(self):
    gw_imean = flex.double(self.I.size(), 0.0)
    for hkl in self.ind_unique:
      sel = (self.ind == hkl)
      _I = self.I.select(sel)
      _sigI = self.sigI.select(sel)
      _weights = 1.0 / (_sigI * _sigI)
      _Imean = flex.sum(_weights * _I) / flex.sum(_weights)
      gw_imean.set_selected(sel, _Imean)
    return gw_imean

  def plot_z_histogram(self):
    fig, ax = plt.subplots()


    ax.set_title('Z Histogram')
    ax.set_xlabel('Z')
    ax.set_ylabel('N')
    ax.hist(self.z, label='Z', bins=100, range=(-10, 10))
    fig.savefig(
      os.path.splitext(os.path.basename(self.outfile))[0] + '_zhistogram'
    )
    plt.close()

  def plot_symmetry_equivalents(self,
    overlay_mean=False, overlay_error=False
  ):
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

    fig.savefig(
      os.path.splitext(os.path.basename(self.outfile))[0]
      + '_probplot'
    )
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

    fig.savefig(
      os.path.splitext(os.path.basename(self.outfile))[0]
      + '_deviation_vs_multiplicity'
    )
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

    fig.savefig(
      os.path.splitext(os.path.basename(self.outfile))[0]
      + '_deviation_detector_map'
    )
    plt.close()

  def time_series(self):
    fig, ax = plt.subplots()

    ax.set_title(
      r'Difference between ordered responses, $z$, '
      + r'and order statistic medians, $m$'
    )
    ax.set_xlabel('Approximate chronology (image number)')
    ax.set_ylabel(r'$z - m$')

    ax.plot(self.image, self.z, '.')

    fig.savefig(
      os.path.splitext(os.path.basename(self.outfile))[0]
      + '_deviation_time_series'
    )
    plt.close()

  def deviation_vs_IsigI(self):
    if not (self.osm.any() and self.osr.any()):
      self.probplot()

    fig, ax = plt.subplots()

    ax.set_title(
      r'Difference between ordered responses, $z$, '
      + r'and order statistic medians, $m$'
    )
    ax.set_xlabel(r'\bar{I}_\mathbf{h} / \sigma_\mathbf{h}')
    ax.set_ylabel(r'$z$')
    ax.plot(
      self.Imeans/self.sigImeans,
      #self.I/self.sigI, 
      self.z,
      #self.osr - self.osm,
      '.')

    fig.savefig(
      os.path.splitext(os.path.basename(self.outfile))[0]
      + '_deviation_vs_IsigI'
    )
    plt.close()

if __name__ == "__main__":
  # TODO Handle multiple input MTZ files.
  # TODO Allow determination of output filename root.
  data = DataDist(sys.argv[1]) #Give an unmerged MTZ file as an argument.

  data.probplot()
  data.deviation_vs_multiplicity()
  data.deviation_map()
  data.time_series()
  data.plot_z_histogram()
  data.deviation_vs_IsigI()
