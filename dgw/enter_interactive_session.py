#!/usr/bin/env cctbx.python

#
#  Copyright (C) (2013) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

"""
Enter an interactive (ipython) session with useful models constructed for
testing things.

"""

# Python and cctbx imports
from __future__ import division
import sys
from math import pi
from scitbx import matrix
from libtbx.phil import parse
from libtbx.test_utils import approx_equal
from scitbx.array_family import flex

# Get modules to build models and minimiser using PHIL
from dials.test.algorithms.refinement import setup_geometry
from dials.test.algorithms.refinement import setup_minimiser

# We will set up a mock scan and a mock experiment list
from dxtbx.model.scan import scan_factory
from dxtbx.model.experiment.experiment_list import ExperimentList, Experiment

# Model parameterisations
from dials.algorithms.refinement.parameterisation.detector_parameters import \
    DetectorParameterisationSinglePanel
from dials.algorithms.refinement.parameterisation.beam_parameters import \
    BeamParameterisation
from dials.algorithms.refinement.parameterisation.crystal_parameters import \
    CrystalOrientationParameterisation, CrystalUnitCellParameterisation

# Symmetry constrained parameterisation for the unit cell
from cctbx.uctbx import unit_cell
from rstbx.symmetry.constraints.parameter_reduction import \
    symmetrize_reduce_enlarge

# Reflection prediction
from dials.algorithms.spot_prediction import IndexGenerator
from dials.algorithms.refinement.prediction import ScansRayPredictor, \
  ExperimentsPredictor
from dials.algorithms.spot_prediction import ray_intersection
from cctbx.sgtbx import space_group, space_group_symbols

# Parameterisation of the prediction equation
from dials.algorithms.refinement.parameterisation.prediction_parameters import \
    XYPhiPredictionParameterisation

# Imports for the target function
from dials.algorithms.refinement.target import \
    LeastSquaresPositionalResidualWithRmsdCutoff

from dials.algorithms.refinement.reflection_manager import ReflectionManager

# Import helper functions
from dials.algorithms.refinement.refinement_helpers import print_model_geometry

#############################
# Setup experimental models #
#############################

args = sys.argv[1:]
master_phil = parse("""
    include scope dials.test.algorithms.refinement.geometry_phil
    include scope dials.test.algorithms.refinement.minimiser_phil
    """, process_includes=True)

models = setup_geometry.Extract(master_phil, cmdline_args = args)

mydetector = models.detector
mygonio = models.goniometer
mycrystal = models.crystal
mybeam = models.beam

# Build a mock scan for a 180 degree sweep
sf = scan_factory()
myscan = sf.make_scan(image_range = (1,1800),
                      exposure_times = 0.1,
                      oscillation = (0, 0.1),
                      epochs = range(1800),
                      deg = True)
sweep_range = myscan.get_oscillation_range(deg=False)
temp = myscan.get_oscillation(deg=False)
im_width = temp[1] - temp[0]
assert sweep_range == (0., pi)
assert approx_equal(im_width, 0.1 * pi / 180.)

# Build an experiment list
experiments = ExperimentList()
experiments.append(Experiment(
      beam=mybeam, detector=mydetector, goniometer=mygonio,
      scan=myscan, crystal=mycrystal, imageset=None))

###########################
# Parameterise the models #
###########################

det_param = DetectorParameterisationSinglePanel(mydetector)
s0_param = BeamParameterisation(mybeam, mygonio)
xlo_param = CrystalOrientationParameterisation(mycrystal)
xluc_param = CrystalUnitCellParameterisation(mycrystal)

# Fix beam to the X-Z plane (imgCIF geometry), fix wavelength
s0_param.set_fixed([True, False, True])

# Fix crystal parameters
#xluc_param.set_fixed([True, True, True, True, True, True])

########################################################################
# Link model parameterisations together into a parameterisation of the #
# prediction equation                                                  #
########################################################################

pred_param = XYPhiPredictionParameterisation(experiments,
  [det_param], [s0_param], [xlo_param], [xluc_param])

################################
# Apply known parameter shifts #
################################

# shift detector by 1.0 mm each translation and 2 mrad each rotation
det_p_vals = det_param.get_param_vals()
p_vals = [a + b for a, b in zip(det_p_vals,
                                [1.0, 1.0, 1.0, 2., 2., 2.])]
det_param.set_param_vals(p_vals)

# shift beam by 2 mrad in free axis
s0_p_vals = s0_param.get_param_vals()
p_vals = list(s0_p_vals)

p_vals[0] += 2.
s0_param.set_param_vals(p_vals)

# rotate crystal a bit (=2 mrad each rotation)
xlo_p_vals = xlo_param.get_param_vals()
p_vals = [a + b for a, b in zip(xlo_p_vals, [2., 2., 2.])]
xlo_param.set_param_vals(p_vals)

# change unit cell a bit (=0.1 Angstrom length upsets, 0.1 degree of
# gamma angle)
xluc_p_vals = xluc_param.get_param_vals()
cell_params = mycrystal.get_unit_cell().parameters()
cell_params = [a + b for a, b in zip(cell_params, [0.1, 0.1, 0.1, 0.0,
                                                   0.0, 0.1])]
new_uc = unit_cell(cell_params)
newB = matrix.sqr(new_uc.fractionalization_matrix()).transpose()
S = symmetrize_reduce_enlarge(mycrystal.get_space_group())
S.set_orientation(orientation=newB)
X = tuple([e * 1.e5 for e in S.forward_independent_parameters()])
xluc_param.set_param_vals(X)

#############################
# Generate some reflections #
#############################

print "Reflections will be generated with the following geometry:"
print_model_geometry(mybeam, mydetector, mycrystal)
print "Target values of parameters are"
msg = "Parameters: " + "%.5f " * len(pred_param)
print msg % tuple(pred_param.get_param_vals())
print

# All indices in a 2.0 Angstrom sphere
resolution = 2.0
index_generator = IndexGenerator(mycrystal.get_unit_cell(),
                space_group(space_group_symbols(1).hall()).type(), resolution)
indices = index_generator.to_array()

# Predict rays within the sweep range
ray_predictor = ScansRayPredictor(experiments, sweep_range)
obs_refs = ray_predictor.predict(indices)

print "Total number of reflections excited", len(obs_refs)

# Take only those rays that intersect the detector
intersects = ray_intersection(mydetector, obs_refs)
obs_refs = obs_refs.select(intersects)

# Make a reflection predictor and re-predict for all these reflections. The
# result is the same, but we gain also the flags and xyzcal.px columns
ref_predictor = ExperimentsPredictor(experiments)
obs_refs['id'] = flex.size_t(len(obs_refs), 0)
obs_refs = ref_predictor.predict(obs_refs)

# Set 'observed' centroids from the predicted ones
obs_refs['xyzobs.mm.value'] = obs_refs['xyzcal.mm']

# Invent some variances for the centroid positions of the simulated data
im_width = 0.1 * pi / 180.
px_size = mydetector[0].get_pixel_size()
var_x = flex.double(len(obs_refs), (px_size[0] / 2.)**2)
var_y = flex.double(len(obs_refs), (px_size[1] / 2.)**2)
var_phi = flex.double(len(obs_refs), (im_width / 2.)**2)
obs_refs['xyzobs.mm.variance'] = flex.vec3_double(var_x, var_y, var_phi)

print "Total number of observations made", len(obs_refs)

###############################
# Undo known parameter shifts #
###############################

s0_param.set_param_vals(s0_p_vals)
det_param.set_param_vals(det_p_vals)
xlo_param.set_param_vals(xlo_p_vals)
xluc_param.set_param_vals(xluc_p_vals)

print "Initial values of parameters are"
msg = "Parameters: " + "%.5f " * len(pred_param)
print msg % tuple(pred_param.get_param_vals())
print

#####################################
# Select reflections for refinement #
#####################################

refman = ReflectionManager(obs_refs, experiments)

##############################
# Set up the target function #
##############################

# The current 'achieved' criterion compares RMSD against 1/3 the pixel size and
# 1/3 the image width in radians. For the simulated data, these are just made up

mytarget = LeastSquaresPositionalResidualWithRmsdCutoff(experiments,
    ref_predictor, refman, pred_param)

################################
# Set up the refinement engine #
################################

refiner = setup_minimiser.Extract(master_phil,
                                  mytarget,
                                  pred_param,
                                  cmdline_args = args).refiner

print "Prior to refinement the experimental model is:"
print_model_geometry(mybeam, mydetector, mycrystal)

# get a CS-PAD detector for testing
import os
import libtbx.load_env
dials_regression = libtbx.env.find_in_repositories(
    relative_path="dials_regression",
    test=os.path.isdir)
data_dir = os.path.join(dials_regression, "refinement_test_data",
                        "hierarchy_test")
datablock_path = os.path.join(data_dir, "datablock.json")
assert os.path.exists(datablock_path)

# load models
from dxtbx.datablock import DataBlockFactory
datablock = DataBlockFactory.from_serialized_format(datablock_path, check_format=False)
im_set = datablock[0].extract_imagesets()[0]
from copy import deepcopy
cspad = deepcopy(im_set.get_detector())

# get also a hierarchical type P6M detector
data_dir = os.path.join(dials_regression, "refinement_test_data",
                        "metrology", "i03-2.5A-thaumatin-20140514-split")
datablock_path = os.path.join(data_dir, "datablock.json")
datablock = DataBlockFactory.from_serialized_format(datablock_path, check_format=False)
im_set = datablock[0].extract_imagesets()[0]
from copy import deepcopy
p6m = deepcopy(im_set.get_detector())
print p6m[0]

# get a P12M (cropped to middle 18 modules)
data_dir = os.path.join(dials_regression, "indexing_test_data",
                        "i23-multiple-germanate")
datablock_path = os.path.join(data_dir, "datablock.json")
datablock = DataBlockFactory.from_serialized_format(datablock_path, check_format=False)
im_set = datablock[0].extract_imagesets()[0]
from copy import deepcopy
p12m = deepcopy(im_set.get_detector())
print p12m[0]

# hierarchical parameterisation
from dials.algorithms.refinement.parameterisation.detector_parameters import \
    DetectorParameterisationHierarchical, get_panel_groups_at_depth, \
    get_panel_ids_at_root

# parameterise
cspad_param = DetectorParameterisationHierarchical(cspad,level=0)

# enter interactive console
from dials.util.command_line import interactive_console; interactive_console()
#refiner.run()

#print
#print "Refinement has completed with the following geometry:"
#print_model_geometry(mybeam, mydetector, mycrystal)
