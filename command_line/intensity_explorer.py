#!/usr/bin/env python
#-*- coding:utf-8 -*-

"""Learning a little about the distribution of diffraction spot intensities
and their errors, reading from an unmerged MTZ file."""

# TODO Proper documentation.
# TODO Sprinkle some sensible tests around the place.

import os
import sys
import math
from scipy import stats as ss
from iotbx import mtz
from cctbx.array_family import flex
from matplotlib import colors, pyplot as plt


class DataDist:
  def __init__(self, filename, outfile=False, keep_singles=False, error=None):
    if outfile:
      self.outfile = outfile
    else:
      self.outfile = filename
    (self.ind,
      self.I,
      self.sigI,
      self.x,
      self.y,
      self.image) = self._data_from_unmerged_mtz(filename)
    (self.ind,
      self.I,
      self.sigI,
      self.x,
      self.y,
      self.image,
      self.multis,
      self.ind_unique,
      self.kept_singles) = self._select_by_multiplicity(keep_singles)
    self.Imeans, self.sigImeans, self.stddevs = self._mean_error_stddev()
    self.z, self.order = self._make_z(error)
    self.osm = self._probplot_data()


  def _data_from_unmerged_mtz(self, filename):
    m = mtz.object(filename)  #Parse MTZ, with lots of useful methods.

    ind = m.extract_miller_indices()  #A flex array of Miller indices.

    cols = m.columns()  #Generates columns (augmented flex arrays).
    col_dict = { c.label() : c for c in cols }  #A dict of all the columns.
    I, sigI, x, y, image = (
      col_dict[label].extract_values().as_double()
      for label in ('I', 'SIGI', 'XDET', 'YDET', 'BATCH')
    )

    return ind, I, sigI, x, y, image


  def _select_by_multiplicity(self, keep_singles=False):
    # Find unique Miller indices.
    ind_unique = set(self.ind)
    # Record the multiplicities.
    multis = flex.int(self.ind.size(), 0)

    for hkl in ind_unique:
      sel = (self.ind == hkl).iselection()
      multis.set_selected(sel, sel.size())

    # Drop multiplicity-1 data unless instructed otherwise.
    if not keep_singles:
      sel = (multis != 1).iselection()
      multis = multis.select(sel)
      ind = self.ind.select(sel)
      I = self.I.select(sel)
      sigI = self.sigI.select(sel)
      x = self.x.select(sel)
      y = self.y.select(sel)
      image = self.image.select(sel)
      ind_unique = set(ind)

    return ind, I, sigI, x, y, image, multis, ind_unique, keep_singles


  def _mean_error_stddev(self):
    # Calculate the weighted mean intensities.
    Imeans = flex.double(self.ind.size(), 0)
    # Calculate the standard errors on the means.
    sigImeans = flex.double(self.ind.size(), 0)
    # Calculate the standard deviations from unbiased weighted variances.
    stddevs = flex.double(self.ind.size(), 0)
    # For weighted averaging.
    weights = 1 / flex.pow2(self.sigI)

    for hkl in self.ind_unique:
      sel = (self.ind == hkl).iselection()
      multi = self.multis.select(sel)[0]
      weight = weights.select(sel)
      sum_weights = flex.sum(weight)
      I = self.I.select(sel)
      Imean = flex.sum( weight * I ) / sum_weights
      Imeans.set_selected(sel, Imean)
      sigImean = math.sqrt(1 / sum_weights)
      sigImeans.set_selected(sel, sigImean)
      if multi == 1:
        stddevs.set_selected(sel, self.sigI.select(sel))
      else:
        variance = (
          flex.sum(weight * flex.pow2(I - Imean)) /
          (sum_weights - flex.sum(flex.pow2(weight)) / sum_weights)
        )
        stddev = math.sqrt(variance)
        stddevs.set_selected(sel, stddev)

    return Imeans, sigImeans, stddevs


  def _make_z(self, error='sigma'):
    if error in ('stddev', 'sigImean'):
      error = {'stddev':self.stddevs, 'sigImean':self.sigImeans}[error]
    else:
      error = self.sigI

    z = (self.I - self.Imeans) / error
    order = flex.sort_permutation(z)

    return z, order


  def _probplot_data(self):
    osm = flex.double(self.z.size(), 0)
    osm.set_selected(
      self.order, flex.double(ss.probplot(self.z, fit=False)[0]) )

    return osm


  def plot_z_histogram(self):
    fig, ax = plt.subplots()

    ax.set_title(r'$z$ histogram')
    ax.set_xlabel(r'$z$')
    ax.set_ylabel(r'$N$')
    ax.hist(self.z, label='$z$', bins=100, range=(-5, 5))
    fig.savefig(
      os.path.splitext(os.path.basename(self.outfile))[0] + '_zhistogram',
      transparent = True
    )
    plt.close()


  def _plot_symmetry_equivalents(self,
    overlay_mean=False, overlay_error=False
  ):
    '''Really just a test function.  Slow.  You probably don't want to use.'''
    fig, ax = plt.subplots()
    
    for hkl in self.ind_unique:
      sel = (self.ind == hkl).iselection()

      if overlay_error == 'sigImean':
        yerr = self.sigImeans.select(sel)
      elif overlay_error:
        yerr = self.stddevs.select(sel)
      else:
        yerr = None

      plt.errorbar(
        sel,#self.image.select(sel),
        self.I.select(sel),
        yerr=self.sigI.select(sel),
        ls="--"
      )
      if overlay_mean:
        plt.errorbar(
          sel,#self.image.select(sel),
          self.Imeans.select(sel),
          yerr = yerr,
          ls = "-",
          color = "k",
          lw = .5
        )
      
    #fig.savefig(
    #  os.path.splitext(os.path.basename(self.outfile))[0] + 'testfig'
    #)
    plt.show()


  def probplot(self):
    fig, ax = plt.subplots()

    ax.set_title('Normal probability plot')
    ax.set_xlabel('Order statistic medians, $m$')
    ax.set_ylabel(r'Ordered responses, $z$')
    ax.set_ylim(-10,10)
    ax.plot(self.osm, self.z, '.b')
    ax.plot([-5,5], [-5,5], '-g')

    fig.savefig(
      os.path.splitext(os.path.basename(self.outfile))[0]
      + '_probplot',
      transparent = True
    )
    plt.close()


  def z_vs_multiplicity(self):
    fig, ax = plt.subplots()

    ax.set_title(
      r'$z$-scores versus multiplicity'
    )
    ax.set_xlabel('Multiplicity')
    ax.set_ylabel(r'$z$')
    ax.plot(self.multis, self.z, '.')

    fig.savefig(
      os.path.splitext(os.path.basename(self.outfile))[0]
      + '_z_vs_multiplicity',
      transparent = True
    )
    plt.close()


  def z_map(self, minimum=0):
    sel = (abs(self.z) >= minimum).iselection()

    extreme = math.ceil(flex.max(flex.abs(self.z)))
    norm = colors.SymLogNorm(
      vmin = -extreme, vmax = extreme,
      linthresh = .02, linscale=1,
    )
    cmap_kws = {'cmap' : 'coolwarm_r', 'norm' : norm}

    fig, ax = plt.subplots()

    ax.set_title(
      r'$z$-scores versus detector position'
    )
    ax.set_xlabel('Detector x position (pixels)')
    ax.set_ylabel('Detector y position (pixels)')
    ax.set_aspect('equal', 'box')
    ax.set_xlim(min(self.x)-5, max(self.x)+5)
    ax.set_ylim(min(self.y)-5, max(self.y)+5)
    det_map = ax.scatter(
      self.x.select(sel),
      self.y.select(sel),
      c = self.z.select(sel),
      marker = ',',
      s = 0.5,
      **cmap_kws
    )
    cbar = fig.colorbar(det_map, ax=ax, **cmap_kws)
    cbar.set_label(r'$z$')

    fig.savefig(
      os.path.splitext(os.path.basename(self.outfile))[0]
      + '_z_detector_map',
      transparent = True
    )
    plt.close()


  def time_series(self):
    fig, ax = plt.subplots()

    ax.set_title(
      r'Time series of $z$-scores'
    )
    ax.set_xlabel('Approximate chronology (image number)')
    ax.set_ylabel(r'$z$')

    ax.plot(self.image, self.z, '.')

    fig.savefig(
      os.path.splitext(os.path.basename(self.outfile))[0]
      + '_z_time_series',
      transparent = True
    )
    plt.close()


  def z_vs_IsigI(self):
    fig, ax = plt.subplots()

    ax.set_title(
      r'$z$-scores versus spot intensity'
    )
    ax.set_xlabel(r'$\bar{I}_\mathbf{h} / \sigma_\mathbf{h}$')
    ax.set_ylabel(r'$z$')
    ax.set_ylim(-10,10)
    ax.set_xscale('log')
    ax.plot(
      abs(self.Imeans/self.sigImeans),
      self.z,
      '.')

    fig.savefig(
      os.path.splitext(os.path.basename(self.outfile))[0]
      + '_z_vs_I_over_sigma',
      transparent = True
    )
    plt.close()


  def z_vs_I(self):
    fig, ax = plt.subplots()

    ax.set_title(
      r'$z$-scores versus spot intensity'
    )
    ax.set_xlabel(r'$\bar{I}_\mathbf{h}$')
    ax.set_ylabel(r'$z$')
    ax.set_ylim(-10,10)
    ax.set_xscale('log')
    ax.plot(
      abs(self.Imeans),
      self.z,
      '.')

    fig.savefig(
      os.path.splitext(os.path.basename(self.outfile))[0]
      + '_z_vs_I',
      transparent = True
    )
    plt.close()


if __name__ == "__main__":
  # TODO Handle multiple input MTZ files.
  # TODO Allow determination of output filename root.
  data = DataDist(sys.argv[1]) #Give an unmerged MTZ file as an argument.

  data.probplot()
  data.z_vs_multiplicity()
  data.z_map()
  data.time_series()
  data.plot_z_histogram()
  data.z_vs_IsigI()
  data.z_vs_I()
