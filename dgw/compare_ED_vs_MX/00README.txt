ED vs MX geometry simulated refinement comparison
=================================================

Introduction
------------

Here, investigate differences in behaviour of geometry refinement for typical
electron diffraction geometry compared with X-ray MX. Start where the previous
investigation, performed for the CCP-EM Spring Symposium, left off. From that
take two files:

* ``observed.pickle`` - 55.8 degree sweep of indexed reflections for a lysozyme
  electron diffraction dataset. The original images were provided by Taimin
  Yang at Stockholm University.
* ``experiments_ED.json`` - the associated experimental geometry for the real
  ED experiment. The geometry is complex, with the beam forming an angle of 7.5
  degrees with the detector normal.

I want to produce regularised experimental geometries for both typical MX
geometry and ED geometry, and then produce indexed centroids for both
geometries by prediction and adding some noise.

Regularise input
----------------

The script ``regularise_experiments.py`` reads in ``experiments_ED.json`` and
simplifies the geometry for the ED experiment as follows:

* Force beam direction along -Z axis.
* Create a new detector of the same overall size as the original Timepix
  detector, and at the same distance, but use a single panel instead of 4
  panels, and make it orthogonal to the beam, intersecting in the centre. We
  ignore thickness and material, so assume there is no parallax correction.
* Increase the scan range to a full turn. The simulation of new 'observations'
  from the updated geometry may cause reflections near the edges of the
  original scan to appear outside that scan range (particularly an issue for
  the conversion to MX geometry). During refinement, reflections outside the
  scan range will be discarded, so this ensures that all the original
  observations will be within the new scan range.

The updated geometry is then written to ``experiments_ED_regularised.json``.

Now a simple MX geometry is constructed:

* Change wavelength to 12 keV.
* Construct a single panel detector with dimensions of a Pilatus 6M at a
  distance of 200 mm, with the beam centre in the middle of the panel.

In either case, the rotation axis is not changed from its original direction,
which for ``experiments_ED.json`` is unusual, being about (-0.755, -0.656, 0).
The reason to not change this to something conventional such as (1, 0, 0) is
that the crystal model orientation would have to be changed so that the same
reflections as contained in ``observed.pickle`` would be observed by rotation
around the new axis. This way is simpler and should still demonstrate the
differences between MX and ED refinement that I'm interested in.

Simulate observations
---------------------

The script ``create_indexed.py`` takes the original observations in
``observed.pickle`` along with experimental descriptions in .json files to
generate new centroids for those observations. Here we run it like this::

  dials.python create_indexed.py observed.pickle \
    experiments_ED_regularised.json experiments_MX_regularised.json

This creates two new files, ``experiments_ED_regularised.pickle`` and
``experiments_MX_regularised.pickle`` containing the simulated observations.
These are created by predicting reflections using the regularised geometry in
each case and then adding error to the predicted centroids to form
observations. In either case the error vector is the same in pixels/images,
but this is scaled appropriately into mm/rad for each experiment. The centroid
variances from spot-finding are left untouched for the centroids in
pixels/images, but again these are scaled appropriately for the centroids in
mm/rad. Only reflections that can be predicted with both geometries are written
to the output files. This ensures that refinement can be performed using the
equivalent set of reflections in each case, to help comparison.

Refinement
----------

**TODO**

