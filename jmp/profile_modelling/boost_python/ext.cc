#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dials/error.h>
#include <dials_scratch/jmp/profile_modelling/target.h>
#include <dials_scratch/jmp/profile_modelling/target_forward.h>

using namespace boost::python;

namespace dials {
namespace algorithms {
namespace boost_python {

void CovarianceMatrix_compose(
    CovarianceMatrix &self,
    boost::python::object reciprocal_lattice_point_spread,
    boost::python::object mosaic_block_angular_spread,
    boost::python::object wavelength_spread) {

  // Check the input
  DIALS_ASSERT(boost::python::len(reciprocal_lattice_point_spread) == 6);

  // Construct parameter array
  scitbx::af::shared<double> parameters;

  // Add reciprocal lattice point spread parameters
  for (std::size_t i = 0; i < 6; ++i) {
    parameters.push_back(
        boost::python::extract<double>(reciprocal_lattice_point_spread[i])());
  }

  // Add the mosaic block angular spread parameters
  if (self.use_mosaic_block_angular_spread()) {
    DIALS_ASSERT(mosaic_block_angular_spread != boost::python::object());
    parameters.push_back(
        boost::python::extract<double>(mosaic_block_angular_spread)());
  }

  // Add the wavelength spread parameters
  if (self.use_wavelength_spread()) {
    DIALS_ASSERT(wavelength_spread != boost::python::object());
    parameters.push_back(boost::python::extract<double>(wavelength_spread)());
  }

  // Set the parameters
  self.set_parameters(parameters.const_ref());
}

static scitbx::af::shared<std::string>
MLTarget3D_parameter_names(bool use_mosaic_block_angular_spread,
                           bool use_wavelength_spread) {
  scitbx::af::shared<std::string> result;
  result.push_back("RLP: L00");
  result.push_back("RLP: L10");
  result.push_back("RLP: L11");
  result.push_back("RLP: L20");
  result.push_back("RLP: L21");
  result.push_back("RLP: L22");
  if (use_mosaic_block_angular_spread) {
    result.push_back("Mosaic Block Spread");
  }
  if (use_wavelength_spread) {
    result.push_back("Wavelength Spread");
  }
  return result;
}

BOOST_PYTHON_MODULE(dials_scratch_jmp_profile_modelling_ext) {
  class_<ReciprocalLatticePointSpread>("ReciprocalLatticePointSpread", no_init)
      .def(init<mat3<double>>((arg("A"))))
      .def(init<mat3<double>, double6>((arg("A"), arg("parameters"))))
      .add_property("A", &ReciprocalLatticePointSpread::get_A)
      .add_property("parameters", &ReciprocalLatticePointSpread::get_parameters,
                    &ReciprocalLatticePointSpread::set_parameters)
      .add_property("covariance", &ReciprocalLatticePointSpread::get_covariance,
                    &ReciprocalLatticePointSpread::set_covariance);

  class_<MosaicBlockAngularSpread>("MosaicBlockAngularSpread", no_init)
      .def(init<vec3<double>>((arg("r"))))
      .def(init<vec3<double>, double>((arg("r"), arg("parameter"))))
      .add_property("r", &MosaicBlockAngularSpread::get_r)
      .add_property("parameter", &MosaicBlockAngularSpread::get_parameter,
                    &MosaicBlockAngularSpread::set_parameter)
      .add_property("covariance", &MosaicBlockAngularSpread::get_covariance,
                    &MosaicBlockAngularSpread::set_covariance);

  class_<WavelengthSpread>("WavelengthSpread", no_init)
      .def(init<vec3<double>, vec3<double>>((arg("s0"), arg("r"))))
      .def(init<vec3<double>, vec3<double>, double>(
          (arg("s0"), arg("r"), arg("parameter"))))
      .add_property("s0", &WavelengthSpread::get_s0)
      .add_property("r", &WavelengthSpread::get_r)
      .add_property("parameter", &WavelengthSpread::get_parameter,
                    &WavelengthSpread::set_parameter)
      .add_property("covariance", &WavelengthSpread::get_covariance,
                    &WavelengthSpread::set_covariance);

  class_<CovarianceMatrix>("CovarianceMatrix", no_init)
      .def(init<mat3<double>, vec3<double>, vec3<double>, bool, bool>(
          (arg("A"), arg("s0"), arg("r"),
           arg("use_mosaic_block_angular_spread") = true,
           arg("use_wavelength_spread") = true)))
      .add_property("A", &CovarianceMatrix::get_A)
      .add_property("s0", &CovarianceMatrix::get_s0)
      .add_property("r", &CovarianceMatrix::get_r)
      .add_property("use_mosaic_block_angular_spread",
                    &CovarianceMatrix::use_mosaic_block_angular_spread)
      .add_property("use_wavelength_spread",
                    &CovarianceMatrix::use_wavelength_spread)
      .def("num_parameters", &CovarianceMatrix::num_parameters)
      .add_property("parameters", &CovarianceMatrix::get_parameters,
                    &CovarianceMatrix::set_parameters)
      .add_property("covariance", &CovarianceMatrix::get_covariance)
      .def("reciprocal_lattice_point_spread",
           &CovarianceMatrix::get_reciprocal_lattice_point_spread,
           return_value_policy<copy_const_reference>())
      .def("mosaic_block_angular_spread",
           &CovarianceMatrix::get_mosaic_block_angular_spread,
           return_value_policy<copy_const_reference>())
      .def("wavelength_spread", &CovarianceMatrix::get_wavelength_spread,
           return_value_policy<copy_const_reference>())
      .def("compose", &CovarianceMatrix_compose,
           (arg("reciprocal_lattice_point_spread") = boost::python::object(),
            arg("mosaic_block_angular_spread") = boost::python::object(),
            arg("wavelength_spread") = boost::python::object()));

  class_<Model2D>("Model2D", no_init)
      .def(init<const Beam &, const Panel &, mat3<double>, vec3<double>>(
          (arg("beam"), arg("panel"), arg("sigma"), arg("r0"))))
      .def("f", &Model2D::f)
      .def("J", &Model2D::J)
      .def("coord", &Model2D::coord);

  class_<Model3D>("Model3D", no_init)
      .def(init<const Beam &, const Panel &, const Goniometer &, const Scan &,
                mat3<double>, vec3<double>>((arg("beam"), arg("panel"),
                                             arg("goniometer"), arg("scan"),
                                             arg("sigma"), arg("r0"))))
      .def("f", &Model3D::f)
      .def("J", &Model3D::J)
      .def("coord", &Model3D::coord);

  class_<MLTarget2DSingle>("MLTarget2DSingle", no_init)
      .def(init<const Beam &, const Panel &, mat3<double>, vec3<double>, int6,
                std::size_t>((arg("beam"), arg("panel"), arg("sigma"),
                              arg("r0"), arg("bbox"), arg("num_integral"))))
      .def("log_likelihood", &MLTarget2DSingle::log_likelihood,
           (arg("data"), arg("mask")))
      .def("compute_P", &MLTarget2DSingle::compute_P);

  class_<MLTarget3DSingle>("MLTarget3DSingle", no_init)
      .def(init<const Beam &, const Panel &, const Goniometer &, const Scan &,
                mat3<double>, vec3<double>, int6, std::size_t>(
          (arg("beam"), arg("panel"), arg("goniometer"), arg("scan"),
           arg("sigma"), arg("r0"), arg("bbox"), arg("num_integral"))))
      .def("log_likelihood", &MLTarget3DSingle::log_likelihood,
           (arg("data"), arg("mask")))
      .def("compute_P", &MLTarget3DSingle::compute_P)
      .def("simulate", &MLTarget3DSingle::simulate);

  class_<MLTarget3D>("MLTarget3D", no_init)
      .def(init<const Experiment &, const reflection_table &, std::size_t, bool,
                bool>((arg("experiment"), arg("reflections"),
                       arg("num_integral") = 5,
                       arg("use_mosaic_block_angular_spread") = false,
                       arg("use_wavelength_spread") = false)))
      .def("log_likelihood", &MLTarget3D::log_likelihood, (arg("parameters")))
      .def("covariance", &MLTarget3D::covariance, (arg("parameters")))
      .def("simulate", &MLTarget3D::simulate, (arg("index"), arg("parameters")))
      .def("num_reflections", &MLTarget3D::num_reflections)
      .def("parameter_names", &MLTarget3D_parameter_names,
           (arg("use_mosaic_block_angular_spread") = false,
            arg("use_wavelength_spread") = false))
      .staticmethod("parameter_names");

  class_<ReciprocalLatticePointDistribution>(
      "ReciprocalLatticePointDistribution", no_init)
      .def(init<miller_index, double6>((arg("h0"), arg("parameters"))))
      .add_property("h0", &ReciprocalLatticePointDistribution::h0)
      .add_property("parameters",
                    &ReciprocalLatticePointDistribution::get_parameters,
                    &ReciprocalLatticePointDistribution::set_parameters)
      .add_property("covariance",
                    &ReciprocalLatticePointDistribution::get_covariance,
                    &ReciprocalLatticePointDistribution::set_covariance)
      .def("sample", &ReciprocalLatticePointDistribution::sample);

  class_<WavelengthDistribution>("WavelengthDistribution", no_init)
      .def(init<double, double>((arg("wavelength"), arg("parameters"))))
      .add_property("wavelength", &WavelengthDistribution::wavelength)
      .add_property("parameters", &WavelengthDistribution::get_parameters,
                    &WavelengthDistribution::set_parameters)
      .def("sample", &WavelengthDistribution::sample);

  class_<BeamVectorDistribution>("BeamVectorDistribution", no_init)
      .def(init<vec3<double>, double>((arg("s0"), arg("parameters"))))
      .add_property("s0", &BeamVectorDistribution::s0)
      .add_property("parameters", &BeamVectorDistribution::get_parameters,
                    &BeamVectorDistribution::set_parameters)
      .def("sample", &BeamVectorDistribution::sample);

  class_<MLTargetForward3DSingle>("MLTargetForward3DSingle", no_init)
      .def(init<const Panel &, const Goniometer &, const Scan &, bool, int6,
                ReciprocalLatticePointDistribution, BeamVectorDistribution,
                std::size_t>((arg("beam"), arg("panel"), arg("goniometer"),
                              arg("scan"), arg("entering"), arg("bbox"),
                              arg("reciprocal_lattice_point_distribution"),
                              arg("beam_vector_distribution"),
                              arg("num_integral"))))
      .def("log_likelihood", &MLTargetForward3DSingle::log_likelihood,
           (arg("data"), arg("mask")))
      .def("simulate", &MLTargetForward3DSingle::simulate);
}

} // namespace boost_python
} // namespace algorithms
} // namespace dials
