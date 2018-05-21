#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dials_scratch/jbe/scaling_code/scaling_helper.h>

namespace dials_scratch { namespace scaling { namespace boost_python {

  using namespace boost::python;

  using scitbx::sparse::matrix;

    void export_elementwise_square()
    {
    def("elementwise_square", &elementwise_square, (
      arg("m")));
    }

    void export_sph_harm_table()
    {
    def("create_sph_harm_table", &create_sph_harm_table, (
      arg("s0_theta_phi"),
      arg("s1_theta_phi"),
      arg("lmax")));
    }

    void export_rotate_vectors_about_axis()
    {
    def("rotate_vectors_about_axis", &rotate_vectors_about_axis, (
      arg("rot_axis"),
      arg("vectors"),
      arg("angles")));
    }

    void export_calc_theta_phi()
    {
    def("calc_theta_phi", &calc_theta_phi,(
      arg("xyz")));
    }

    void export_calc_sigmasq()
    {
    def("calc_sigmasq", &calc_sigmasq, (
      arg("jacobian_transpose"),
      arg("var_cov")));
    }

    void export_row_multiply()
    {
    def("row_multiply", &row_multiply, (
      arg("m"),
      arg("v")));
    }

}}} // namespace = dials_scratch::scaling::boost_python
