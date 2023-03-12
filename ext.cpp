#include <boost/python.hpp>
#include <cctype>
#include <dxtbx/model/beam.h>
#include <dxtbx/model/detector.h>
#include <dxtbx/model/panel.h>
#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/mat3.h>
#include <scitbx/vec3.h>

namespace dials_scratch {
namespace ext {

// make a python list

static boost::python::list make_list(size_t n) {
  boost::python::list result;
  for (size_t i = 0; i < n; i++) {
    result.append(i);
  }
  return result;
}

// make a flex array (much more flexible)

static scitbx::af::shared<int> make_flex(size_t n) {
  scitbx::af::shared<int> result;
  for (size_t i = 0; i < n; i++) {
    result.push_back(i);
  }
  return result;
}

// using flex arrays

static int sum(scitbx::af::shared<int> array) {
  int result = 0;
  for (size_t i = 0; i < array.size(); i++) {
    result += array[i];
  }
  return result;
}

// use dxtbx things

std::string detector_as_string(const dxtbx::model::Detector &detector) {
  std::stringstream ss;
  ss << detector;
  return ss.str();
}

// TODO in here implement a mosaicity tensor rather than an isotropic
// r factor => constrain the mosaicity by the crystal symmetry (external
// to this routine)

scitbx::af::versa<double, scitbx::af::c_grid<2>>
x_map(const dxtbx::model::Panel &panel, const dxtbx::model::Beam &beam,
      const scitbx::mat3<double> &UB_inv, int oversample, double r,
      double d_min) {
  size_t width = panel.get_image_size()[0];
  size_t height = panel.get_image_size()[1];

  scitbx::af::versa<double, scitbx::af::c_grid<2>> map;
  map.resize(scitbx::af::c_grid<2>(height, width));

  scitbx::af::tiny<double, 2> xy;
  scitbx::vec3<double> s0(beam.get_s0());
  double winv = 1.0 / beam.get_wavelength();

  double r2 = r * r;

  size_t offset = 0;
  for (size_t j = 0; j < height; j++) {
    for (size_t i = 0; i < width; i++) {

      double value = 0.0;

      if (d_min > 0) {
        xy[0] = i + 0.5;
        xy[1] = j + 0.5;
        double d = panel.get_resolution_at_pixel(s0, xy);
        if (d < d_min) {
          map[offset] = 0.0;
          offset++;
          continue;
        }
      }

      for (size_t _j = 0; _j < oversample; _j++) {
        for (size_t _i = 0; _i < oversample; _i++) {

          xy[0] = i + ((_i + 0.5) / oversample);
          xy[1] = j + ((_j + 0.5) / oversample);

          scitbx::vec3<double> p(panel.get_pixel_lab_coord(xy));
          scitbx::vec3<double> q = p.normalize() * winv - s0;
          scitbx::vec3<double> hkl = UB_inv * q;
          double d2 = (hkl[0] - round(hkl[0])) * (hkl[0] - round(hkl[0])) +
                      (hkl[1] - round(hkl[1])) * (hkl[1] - round(hkl[1])) +
                      (hkl[2] - round(hkl[2])) * (hkl[2] - round(hkl[2]));

          // TODO here figure out the tensor description of this i.e.
          // ([R]d) ** 2
          value += exp(-d2 / r2);
        }
      }

      map[offset] = value / (oversample * oversample);
      offset++;
    }
  }

  return map;
}

scitbx::af::versa<double, scitbx::af::c_grid<2>>
q_map(const dxtbx::model::Panel &panel, const dxtbx::model::Beam &beam,
      const scitbx::mat3<double> &UB, int oversample) {
  size_t width = panel.get_image_size()[0];
  size_t height = panel.get_image_size()[1];

  scitbx::mat3<double> UB_inv = UB.inverse();

  scitbx::af::versa<double, scitbx::af::c_grid<2>> map;
  map.resize(scitbx::af::c_grid<2>(height, width));

  scitbx::af::tiny<double, 2> xy;
  scitbx::vec3<double> s0(beam.get_s0());
  double winv = 1.0 / beam.get_wavelength();

  size_t offset = 0;
  for (size_t j = 0; j < height; j++) {
    for (size_t i = 0; i < width; i++) {

      double value = 0.0;

      for (size_t _j = 0; _j < oversample; _j++) {
        for (size_t _i = 0; _i < oversample; _i++) {

          xy[0] = i + ((_i + 0.5) / oversample);
          xy[1] = j + ((_j + 0.5) / oversample);

          scitbx::vec3<double> p(panel.get_pixel_lab_coord(xy));
          scitbx::vec3<double> q = p.normalize() * winv - s0;
          scitbx::vec3<double> hkl = UB_inv * q;
          scitbx::vec3<int> ihkl;
          ihkl[0] = round(hkl[0]);
          ihkl[1] = round(hkl[1]);
          ihkl[2] = round(hkl[2]);
          scitbx::vec3<double> mq = UB * ihkl;
          value += (q - mq).length();
        }
      }

      map[offset] = value / (oversample * oversample);
      offset++;
    }
  }

  return map;
}

void init_module() {
  using namespace boost::python;
  def("make_list", make_list, (arg("size")));
  def("make_flex", make_flex, (arg("size")));
  def("sum", sum, (arg("array")));
  def("detector_as_string", detector_as_string, (arg("detector")));
  def("x_map", x_map,
      (arg("panel"), arg("beam"), arg("UB_inv"), arg("oversample"), arg("r")));
  def("q_map", q_map,
      (arg("panel"), arg("beam"), arg("UB"), arg("oversample")));
}

} // namespace ext
} // namespace dials_scratch

BOOST_PYTHON_MODULE(dials_scratch_ext) { dials_scratch::ext::init_module(); }
