#!/usr/bin/env python
#-*- coding:utf-8 -*-
"""
Examine the distribution of diffraction spot intensities.

This module defines a class IntensityDist, with several methods for exploring
the distribution of measured spot intensities in an X-ray diffraction
experiment.  The user may wish to use this information to inform decisions
regarding the error model employed in analysing the data.  Data are passed in
as an unmerged MTZ file (see http://www.ccp4.ac.uk/html/mtzformat.html) and the
resulting IntensityDist instance contains the pertinent columns of data, along
with normal order statistic medians of the z-scores of the intensities, for constructing a normal probability plot (See
https://www.itl.nist.gov/div898/handbook/eda/section3/normprpl.htm).

If called as a script, read data from an unmerged MTZ file; generate a
histogram and a normal probability plot of the z-scores of the intensity data,
along with plots of z as a function of batch number, of multiplicity, of
detector position, of measured multiplicity, of absolute intensity and of
I/sigma.

Example:
  $ dials.python intensity_explorer.py <unmerged MTZ file>
"""

# FIXME Docstrings are in Google-ish format — move to Sphinx-ish.
# TODO Once ∃ a dials tool for (unmerged MTZ) –> (exp list, refl table), use it

from __future__ import absolute_import, division, print_function
import os
import sys
import math
from scipy import stats as ss
from iotbx import mtz
from cctbx.array_family import flex
from matplotlib import colors, pyplot as plt


class IntensityDist(object):
  """
  Store intensity data and generate normal order statistic medians.
  
  Attributes:
    ind
      (cctbx_array_family_flex_ext.miller_index): Miller indices.
    I (cctbx_array_family_flex_ext.double):       Measured intensity data.
    sigI (cctbx_array_family_flex_ext.double):    Measured intensity standard
                                                    deviations.
    x (cctbx_array_family_flex_ext.double):       Detector position, x (fast)
                                                    axis component.
    y (cctbx_array_family_flex_ext.double):       Detector position, y (slow)
                                                    axis component.
    image (cctbx_array_family_flex_ext.int):      Batch (image) number.
    multis (cctbx_array_family_flex_ext.int):     Measured multiplicity of
                                                    symmetry-equivalent spots.
    Imeans (cctbx_array_family_flex_ext.double):  Weighted means of symmetry-
                                                    equivalent reflection
                                                    intensities.
    sigImeans
      (cctbx_array_family_flex_ext.double):   Standard deviation on the 
                                                weighted mean intensities.
    stddevs
      (cctbx_array_family_flex_ext.double):   Sample standard deviations,
                                                calculated as square-root of
                                                unbiased weighted sample
                                                variances of symmetry-
                                                equivalent reflection
                                                intensities.
    z (cctbx_array_family_flex_ext.double):   z-scores of weighted mean
                                                intensities.
    order
      (scitbx_array_family_flex_ext.size_t):  Index with which to sort the
                                                z-scores in ascending order.
                                                Useful for making a normal
                                                probability plot.
    osm (cctbx_array_family_flex_ext.double): Normal order statistic medians of
                                                the z-scores.
    ind_unique (set):     Set of observed symmetry-inequivalent Miller indices.
    kept_singles (bool):  Indicates whether multiplicity-1 reflections were
                            retained.
                            Defaults to False.
    outfile (str):        File root for generated plots.
                            Defaults to MTZ input file root.
  """

  def __init__(self, filename, outfile=None, keep_singles=False, 
              error='sigma'):
    """
    Generate z-scores and normal probability plot from an unmerged MTZ file
    
    Args:
      filename (str):       Unmerged MTZ input file.
      outfile (str):        File root for output PNG plots.  If None, the root
                              of the input filename is used.
                              Defaults to None.
      keep_singles (bool):  Choose whether to keep multiplicity-1 reflections.
                              Defaults to False.
      error (str):          Measure of spread to use in normalising the
                              z-scores, i.e. z = (I - <I>) / error.
        Possible values for error:
        'sigma':    Use measured sigma values;
        'stddev':   Use sample standard deviations calculated as square-root of
          unbiased weighted sample variances of symmetry-equivalent reflection
          intensities;
        'sigImean': Use standard deviation on the weighted mean intensities.
          Mathematically meaningless, this is just for debugging.
        Defaults to 'sigma'.
    """
    self.ind = None
    self.I = None
    self.sigI =None
    self.x = None
    self.y = None
    self.image = None
    self.multis = None
    self.Imeans = None
    self.sigImeans = None
    self.stddevs = None
    self.z = None
    self.order = None
    self.osm = None
    self.ind_unique = None
    self.kept_singles = None
    self.outfile = None

    if outfile:
      self.outfile = os.path.splitext(os.path.basename(outfile))[0]
    else:
      self.outfile = os.path.splitext(os.path.basename(filename))[0]
    (self.ind, self.I, self.sigI, self.x, self.y,
      self.image) = self._data_from_unmerged_mtz(filename)
    (self.ind, self.I, self.sigI,
      self.x, self.y, self.image,
      self.multis, self.ind_unique,
      self.kept_singles) = self._select_by_multiplicity(keep_singles)
    self.Imeans, self.sigImeans, self.stddevs = self._mean_error_stddev()
    self.z, self.order = self._make_z(error)
    self.osm = self._probplot_data()


  def _data_from_unmerged_mtz(self, filename):
    m = mtz.object(filename)  #Parse MTZ, with lots of useful methods.

    ind = m.extract_miller_indices()  #A flex array of Miller indices.

    cols = m.columns()  #Generates columns (augmented flex arrays).
    col_dict = { c.label() : c for c in cols }  #A dict of all the columns.
    I, sigI, x, y = (
      col_dict[label].extract_values().as_double()
      for label in ('I', 'SIGI', 'XDET', 'YDET')
    )
    image = col_dict['BATCH'].extract_values().as_double().iround()

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
    """Plot a hitogram of the z-scores of the weighted mean intensities."""
    fig, ax = plt.subplots()

    ax.set_title(r'$z$ histogram')
    ax.set_xlabel(r'$z$')
    ax.set_ylabel(r'$N$')
    ax.hist(self.z, label='$z$', bins=100, range=(-5, 5))
    fig.savefig(self.outfile + '_zhistogram', transparent = True)
    plt.close()


  def _plot_symmetry_equivalents(self,
    overlay_mean=False, overlay_error=False
  ):
    """Really just a test plot.  Slow.  You probably don't want to use."""
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
      
    #fig.savefig(self.outfile + 'testfig')
    plt.show()


  def probplot(self, **kwargs):
    """Create a normal probability plot from the z-scores."""
    fig, ax = plt.subplots()

    ax.set_title('Normal probability plot')
    ax.set_xlabel('Order statistic medians, $m$')
    ax.set_ylabel(r'Ordered responses, $z$')
    ax.set_ylim(-10,10)
    ax.plot(self.osm, self.z, '.b', **kwargs)
    ax.plot([-5,5], [-5,5], '-g')

    fig.savefig(self.outfile + '_probplot', transparent = True)
    plt.close()


  def plot_z_vs_multiplicity(self, **kwargs):
    "Plot intensity z-scores versus multiplicity."
    fig, ax = plt.subplots()

    ax.set_title(
      r'$z$-scores versus multiplicity'
    )
    ax.set_xlabel('Multiplicity')
    ax.set_ylabel(r'$z$')
    ax.plot(self.multis, self.z, '.', **kwargs)

    fig.savefig(self.outfile + '_z_vs_multiplicity', transparent = True)
    plt.close()


  def plot_z_map(self, minimum=0):
    """Plot a z-score heatmap of the detector.
    
    Beware, this is only meaningful if the data have a single geometry model.
    """
    sel = (flex.abs(self.z) >= minimum).iselection()

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
    ax.set_xlim(flex.min(self.x)-5, flex.max(self.x)+5)
    ax.set_ylim(flex.min(self.y)-5, flex.max(self.y)+5)
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

    fig.savefig(self.outfile + '_z_detector_map', transparent = True)
    plt.close()


  def plot_time_series(self, **kwargs):
    """Plot a crude time series of z-scores.
    
    Batch (image) number is used as a proxy for time."""
    fig, ax = plt.subplots()

    ax.set_title(
      r'Time series of $z$-scores'
    )
    ax.set_xlabel('Approximate chronology (image number)')
    ax.set_ylabel(r'$z$')

    ax.plot(self.image, self.z, '.', **kwargs)

    fig.savefig(self.outfile + '_z_time_series', transparent = True)
    plt.close()


  def plot_z_vs_IsigI(self, **kwargs):
    """Plot z-scores versus I/sigma."""
    fig, ax = plt.subplots()

    ax.set_title(
      r'$z$-scores versus spot intensity'
    )
    ax.set_xlabel(r'$\bar{I}_\mathbf{h} / \sigma_\mathbf{h}$')
    ax.set_ylabel(r'$z$')
    ax.set_ylim(-10,10)
    ax.set_xscale('log')
    ax.plot(flex.abs(self.Imeans/self.sigImeans), self.z, '.', **kwargs)

    fig.savefig(self.outfile + '_z_vs_I_over_sigma', transparent = True)
    plt.close()


  def plot_z_vs_I(self, **kwargs):
    """Plot z-scores versus absolute intensity."""
    fig, ax = plt.subplots()

    ax.set_title(
      r'$z$-scores versus spot intensity'
    )
    ax.set_xlabel(r'$\bar{I}_\mathbf{h}$')
    ax.set_ylabel(r'$z$')
    ax.set_ylim(-10,10)
    ax.set_xscale('log')
    ax.plot(flex.abs(self.Imeans), self.z, '.', **kwargs)

    fig.savefig(self.outfile + '_z_vs_I', transparent = True)
    plt.close()


if __name__ == "__main__":
  # TODO Handle multiple input MTZ files.
  # TODO Allow determination of output filename root.
  data = IntensityDist(sys.argv[1]) #Give an unmerged MTZ file as an argument.

  data.plot_z_histogram()
  data.probplot()
  data.time_series()
  data.z_map()
  data.z_vs_multiplicity()
  data.z_vs_I()
  data.z_vs_IsigI()