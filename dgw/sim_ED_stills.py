"""Script for simulating still images with electron diffraction geometry
and a tetragonal lyzozyme cell"""

from scitbx.array_family import flex
from simtbx.nanoBragg import shapetype
from simtbx.nanoBragg import nanoBragg

from scitbx import matrix
from cctbx.eltbx import henke
from dxtbx.model.detector import DetectorFactory
from dxtbx.model.beam import BeamFactory
from dxtbx.model import Crystal
from dxtbx.model.experiment_list import Experiment
from dxtbx.model.experiment_list import ExperimentList
from multiprocessing import Pool
import argparse
from random import uniform, seed

pdb_lines = """HEADER TEST
CRYST1   78.840   78.840   38.290  90.00  90.00  90.00 P 43 21 2
ATOM      1  O   HOH A   1      56.829   2.920  55.702  1.00 20.00           O
ATOM      2  O   HOH A   2      49.515  35.149  37.665  1.00 20.00           O
ATOM      3  O   HOH A   3      52.667  17.794  69.925  1.00 20.00           O
ATOM      4  O   HOH A   4      40.986  20.409  18.309  1.00 20.00           O
ATOM      5  O   HOH A   5      46.896  37.790  41.629  1.00 20.00           O
ATOM      6 SED  MSE A   6       1.000   2.000   3.000  1.00 20.00          SE
END
"""


class Simulation(object):
    def __init__(self, start, num):

        self.start = start
        self.num = num
        assert self.start > 0 and self.num > 1

        self._setup_models()

    def _setup_models(self):

        # Set up detector
        image_size = (1024, 1024)
        beam_centre_mm = (
            self.pixel_size * image_size[0] / 2,
            self.pixel_size * image_size[1] / 2,
        )
        self.detector = DetectorFactory().simple(
            "PAD",
            self.distance,
            beam_centre_mm,
            "+x",
            "-y",
            (self.pixel_size, self.pixel_size),
            image_size,
            trusted_range=(-1, 1000000),
        )

        self.beam = BeamFactory().make_polarized_beam(
            sample_to_source=(0.0, 0.0, 1.0),
            wavelength=self.wavelength,
            polarization=(0, 1, 0),
            polarization_fraction=self.polarization_fraction,
        )

        # Set up simulated structure factors
        self.sfall = self.fcalc_from_pdb(resolution=2.0)

        self._unit_cell = self.cell_from_pdb()
        a, b, c, aa, bb, cc = self._unit_cell.parameters()
        self.experiments = ExperimentList()

        for _ in range(self.num):
            # Set up crystal - does not need to be correct, it is overwritten anyway
            crystal = Crystal(
                real_space_a=(a, 0, 0),
                real_space_b=(0, b, 0),
                real_space_c=(0, 0, c),
                space_group_symbol="P 43 21 2",
            )

            self.experiments.append(
                Experiment(beam=self.beam, detector=self.detector, crystal=crystal)
            )

    def set_imageset(self, filename, expr_no):
        from dxtbx.format.FormatSMVJHSim import FormatSMVJHSim

        exp = self.experiments[expr_no]
        imset = FormatSMVJHSim.get_imageset(
            [filename],
            beam=exp.beam,
            detector=exp.detector,
            goniometer=exp.goniometer,
            scan=exp.scan,
        )
        exp.imageset = imset

    def cell_from_pdb(self):
        from iotbx import pdb

        pdb_inp = pdb.input(source_info=None, lines=pdb_lines)
        return pdb_inp.crystal_symmetry().unit_cell()

    def fcalc_from_pdb(self, resolution):
        from iotbx import pdb

        pdb_inp = pdb.input(source_info=None, lines=pdb_lines)
        xray_structure = pdb_inp.xray_structure_simple()
        wavelength = self.wavelength

        # Assuming X-ray scattering here - does not matter for the geometry. How
        # would one use electron scattering instead?
        scatterers = xray_structure.scatterers()
        for sc in scatterers:
            expected_henke = henke.table(sc.element_symbol()).at_angstrom(wavelength)
            sc.fp = expected_henke.fp()
            sc.fdp = expected_henke.fdp()

        primitive_xray_structure = xray_structure.primitive_setting()
        P1_primitive_xray_structure = primitive_xray_structure.expand_to_p1()
        fcalc = P1_primitive_xray_structure.structure_factors(
            d_min=resolution, anomalous_flag=True, algorithm="direct"
        ).f_calc()

        return fcalc.amplitudes()

    def generate_image(self, image_no=1):

        expr_no = image_no - self.start
        crystal = self.experiments[int(expr_no)].crystal

        # Set the beam for this image
        beam = self.beam

        print("Generating image {0}".format(image_no))
        print(
            "DIALS beam centre will be", self.detector[0].get_beam_centre(beam.get_s0())
        )

        # Construct simulation. Ugh, it prints a load of stuff from C++ that is not
        # easy to suppress.
        print("Ignore the following output from simtbx")
        print("#######################################")
        SIM = nanoBragg(self.detector, beam, verbose=1)
        print("#######################################")

        # Set Ncells to give approx 200nm cube
        a, b, c, aa, bb, cc = self._unit_cell.parameters()
        Na = int(round(2000 / a))
        Nb = int(round(2000 / b))
        Nc = int(round(2000 / c))
        print("setting Ncells:", (Na, Nb, Nc))
        SIM.Ncells_abc = (Na, Nb, Nc)
        # set different random number seed for noise generation for each image
        SIM.seed = image_no
        SIM.oversample = 1
        SIM.progress_meter = False  # only really useful for long runs
        # SIM.default_F=100 # this will become F000, marking the beam center
        # use crystal structure to initialize Fhkl array
        SIM.Fhkl = self.sfall

        # This does not 'stick': "ERROR: cannot initialize without a cell"
        # SIM.Amatrix = self.crystal.get_A()

        # WORKAROUND: Instead, use nanoBragg to set the A matrix by missets, then
        # update the dxtbx crystal to keep track
        seed(image_no)
        SIM.missets_deg = (uniform(0, 360), uniform(0, 360), uniform(0, 360))
        # Apparently this needs the transpose
        crystal.set_A(matrix.sqr(SIM.Amatrix).transpose())

        SIM.xtal_shape = shapetype.Gauss  # fastest option, least realistic
        SIM.flux = 1e12  # quanta/s
        SIM.beamsize_mm = self.beamsize_mm  # assumes round beam
        SIM.exposure_s = 0.1

        SIM.divergence_hv_mrad = self.divergence_hv_mrad
        SIM.divsteps_hv = 6, 6

        # Set mosaic spread _before_ setting the number of domains.  If the
        # mosaicity is zero, the domain count is always reset to 1.
        SIM.mosaic_spread_deg = 0.1
        SIM.mosaic_domains = 10

        # Set detector noise and offset parameters to zero
        SIM.adc_offset_adu = 0
        SIM.readout_noise_adu = 0
        SIM.detector_calibration_noise_pct = 0

        # Now actually burn up some CPU
        SIM.add_nanoBragg_spots()

        # Amplify spot signal
        SIM.raw_pixels *= self.amplifying_factor

        # Write out the noise-free image with a pedestal matching that reported in
        # the header
        SIM.raw_pixels += SIM.adc_offset_adu
        fileout = "intimage_{0:05d}.img".format(image_no)
        SIM.to_smv_format(fileout=fileout, intfile_scale=1)
        SIM.raw_pixels -= SIM.adc_offset_adu

        # Add background scatter: interpolation points for sin(theta/lambda)
        # vs structure factor. Model ED images approximately using an
        # exponential fall off
        stol = flex.double_range(0, 51, 2) / 100
        scatt = 70 * flex.exp(-7 * stol)
        bg = flex.vec2_double(stol, scatt)
        SIM.Fbg_vs_stol = bg
        SIM.amorphous_sample_thick_mm = self.amorphous_sample_thick_mm
        SIM.amorphous_density_gcm3 = 1
        SIM.amorphous_molecular_weight_Da = 18
        SIM.add_background()

        # Add scatter from rough approximation to air
        # bg = flex.vec2_double([(0,14.1),(0.045,13.5),(0.174,8.35),(0.35,4.78),
        #     (0.5,4.22)])
        # SIM.Fbg_vs_stol = bg
        # SIM.amorphous_sample_thick_mm = 35 # between beamstop and collimator
        # SIM.amorphous_density_gcm3 = 1.2e-3
        # SIM.amorphous_sample_molecular_weight_Da = 28 # nitrogen = N2
        # SIM.add_background()

        # detector PSF params for apply_psf(), called within add_noise()
        SIM.detector_psf_kernel_radius_pixels = 0
        SIM.detector_psf_fwhm_mm = 0.05
        SIM.detector_psf_type = shapetype.Gauss

        # Add the various noise sources
        SIM.add_noise()

        # Write out image with noise
        fileout = "noiseimage_{0:05d}.img".format(image_no)
        SIM.to_smv_format(fileout=fileout, intfile_scale=1)

        SIM.free_all()

        # Set an imageset in the experiment list using the noiseimage
        self.set_imageset(fileout, expr_no)

        self.experiments[expr_no : expr_no + 1].as_file(
            "experiments_%05d.json" % image_no
        )

    def generate_all_images(self):
        for image in range(self.start, self.start + self.num):
            self.generate_image(image)

    def generate_all_in_parallel(self, proc):
        with Pool(proc) as pool:
            pool.map(self.generate_image, range(self.start, self.start + self.num))


class SimulationED(Simulation):
    def __init__(self, start, num):

        self.start = start
        self.num = num
        assert self.start > 0 and self.num > 1

        self.distance = 1590.00
        self.pixel_size = 0.055

        # unpolarized 200 keV beam
        self.wavelength = 0.02508
        self.polarization_fraction = 0.5

        # Values found by trial-and-error to produce reasonable images
        self.beamsize_mm = 0.01
        self.divergence_hv_mrad = 0.07, 0.07
        self.amplifying_factor = 100
        self.amorphous_sample_thick_mm = 0.3

        self._setup_models()


class SimulationMX(Simulation):
    def __init__(self, start, num):

        self.start = start
        self.num = num
        assert self.start > 0 and self.num > 1

        self.distance = 100.0
        self.pixel_size = 0.172

        # Unpolarized 12 keV beam
        self.wavelength = 1.0332
        self.polarization_fraction = 0.5

        # Values found by trial-and-error to produce reasonable images
        self.beamsize_mm = 0.1
        self.divergence_hv_mrad = 3.2, 3.2
        self.amplifying_factor = 1e7
        self.amorphous_sample_thick_mm = 0.0001

        self._setup_models()


def run():
    p = argparse.ArgumentParser()
    p.add_argument("start", type=int)
    p.add_argument("num", default=1, type=int)
    p.add_argument("proc", type=int)
    args = p.parse_args()

    # Edit here to SimulationMX if required
    sim = SimulationED(int(args.start), int(args.num))
    sim.generate_all_in_parallel(args.proc)


if __name__ == "__main__":
    run()
