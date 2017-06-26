#include <random>
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <scitbx/constants.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/array_family/reflection_table.h>
#include <dxtbx/model/detector.h>
#include <dxtbx/model/scan.h>
#include <dxtbx/model/crystal.h>
#include <dials/algorithms/spot_prediction/pixel_to_miller_index.h>
#include <dials/algorithms/background/glm/robust_poisson_mean.h>
#include <dials/algorithms/integration/sum/summation.h>
#include <dials/model/data/mask_code.h>
#include <dials/error.h>

using namespace boost::python;

namespace dials { namespace algorithms { namespace boost_python {

  using scitbx::constants::pi;
  using dials::model::Valid;
  using dials::model::Background;
  using dials::model::Foreground;

  namespace detail {

    double min4(double a, double b, double c, double d) {
      return std::min(
          std::min(a, b),
          std::min(c, d));
    }

    double max4(double a, double b, double c, double d) {
      return std::max(
          std::max(a, b),
          std::max(c, d));
    }

    template <typename T>
    T median(const af::const_ref<T> &x) {
      af::shared<T> temp(x.begin(), x.end());
      std::nth_element(temp.begin(), temp.begin() + temp.size() / 2, temp.end());
      return temp[temp.size() / 2];
    }

  }


  class Simulator {
  public:

    Simulator(
        const Beam &beam,
        const Detector &detector,
        const Crystal &crystal,
        double mosaicity,
        std::size_t num_sample,
        const af::const_ref<int,af::c_grid<2> > &data,
        const af::const_ref<bool,af::c_grid<2> > &input_mask) {

      DIALS_ASSERT(mosaicity > 0);
      DIALS_ASSERT(num_sample > 0);

      std::size_t xsize = detector[0].get_image_size()[0];
      std::size_t ysize = detector[0].get_image_size()[1];

      model_ = af::versa< double, af::c_grid<2> >(af::c_grid<2>(ysize, xsize), 0);
      mask_ = af::versa< int, af::c_grid<2> >(af::c_grid<2>(ysize, xsize), 0);
      labels_ = af::versa< int, af::c_grid<2> >(af::c_grid<2>(ysize, xsize), 0);

      std::map< cctbx::miller::index<>, std::vector<std::size_t> > index_pixel_map;

      // Initialise the transform
      PixelToMillerIndex transform(beam, detector, crystal);

      // Compute the normal constant
      double K = std::pow((1.0 / std::sqrt(2*pi) * mosaicity), 3);

      for (std::size_t j = 0; j < ysize; ++j) {
        for (std::size_t i = 0; i < xsize; ++i) {

          // Map the corners of the pixel
          vec3<double> A = transform.h(0, i, j);
          vec3<double> B = transform.h(0, i+1, j);
          vec3<double> C = transform.h(0, i, j+1);
          vec3<double> D = transform.h(0, i+1, j+1);
          vec3<double> E = transform.h(0, i+0.5, j+0.5);

          // The integer miller index
          vec3<double> h0(
              std::floor(E[0] + 0.5),
              std::floor(E[1] + 0.5),
              std::floor(E[2] + 0.5));

          // The distance from all corners
          double distance_E = (E - h0).length();

          // Compute the approximate area of the pixel in reciprocal space
          double area_abc = 0.5 * (B-A).cross(C-A).length();
          double area_bcd = 0.5 * (B-D).cross(C-D).length();
          double area = area_abc + area_bcd;

          // Integrate the 3D normal distribution over the pixel
          double sum_f = 0.0;
          for (std::size_t jj = 0; jj < num_sample; ++jj) {
            for (std::size_t ii = 0; ii < num_sample; ++ii) {
              vec3<double> h = transform.h(0,
                  i + ii / (double)num_sample,
                  j + jj / (double)num_sample);

              double distance = (h - h0).length();

              sum_f += K * std::exp(-0.5 * std::pow(distance / mosaicity, 2));
            }
          }

          if (input_mask(j,i)) {
            mask_(j,i) = Valid;
          }

          if (distance_E < 3 * mosaicity) {
            mask_(j,i) |= Foreground;
          } else if (distance_E < 6 * mosaicity) {
            mask_(j,i) |= Background;
          }

          model_(j,i) = area * sum_f / (double)(num_sample * num_sample);
          index_pixel_map[h0].push_back(i + j * xsize);
        }
      }

      typedef std::map< cctbx::miller::index<>, std::vector<std::size_t> >::iterator iterator;
      /* std::size_t count = 0; */
      /* for (iterator it = index_pixel_map.begin(); it != index_pixel_map.end(); ++it) { */
      /*   indices_.push_back(it->first); */
      /*   for (std::vector<std::size_t>::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2) { */
      /*     labels_[*it2] = count; */
      /*   } */

      /*   count += 1; */
      /* } */

      af::shared< cctbx::miller::index<> > miller_index;
      af::shared< double> intensity;
      af::shared< double> variance;
      af::shared< double> scale;
      af::shared< int> flags;
      af::shared< vec3<double> > xyzcal;
      af::shared< vec3<double> > xyzobs;


      for (iterator it = index_pixel_map.begin(); it != index_pixel_map.end(); ++it) {

        cctbx::miller::index<> h = it->first;
        const std::vector<std::size_t> &pixels = it->second;

        // Extract pixels
        af::shared<double> reflection_data(pixels.size());
        af::shared<double> reflection_pred(pixels.size());
        af::shared<int> reflection_mask(pixels.size());
        af::shared< vec3<double> > reflection_pnts(pixels.size());
        for (std::size_t i = 0; i < pixels.size(); ++i) {
          int index = pixels[i];
          double x = (index % xsize) + 0.5;
          double y = ((int)(index / xsize)) + 0.5;
          reflection_data[i] = data[pixels[i]];
          reflection_mask[i] = mask_[pixels[i]];
          reflection_pred[i] = model_[pixels[i]];
          reflection_pnts[i] = vec3<double>(x, y, 0);
        }

        try {

          // Compute the background
          double B = compute_background(
              reflection_data.const_ref(),
              reflection_mask.const_ref());

          // Compute the intensity
          vec2<double> IV = compute_intensity(
              reflection_data.const_ref(),
              reflection_mask.const_ref(),
              B);

          double I = IV[0];
          double V = IV[1];
          double S = compute_scale(
              reflection_pred.const_ref(),
              reflection_mask.const_ref());

          vec3<double> xcal = compute_centroid(
              reflection_pred.const_ref(),
              reflection_mask.const_ref(),
              reflection_pnts.const_ref());

          vec3<double> xobs = compute_centroid(
              reflection_data.const_ref(),
              reflection_mask.const_ref(),
              reflection_pnts.const_ref());


          miller_index.push_back(h);
          intensity.push_back(I);
          variance.push_back(V);
          scale.push_back(S);
          flags.push_back(af::IntegratedSum);
          xyzcal.push_back(xcal);
          xyzobs.push_back(xobs);
        } catch (dials::error e) {
          std::cout << e.what() << std::endl;
          continue;
        }
      }

      DIALS_ASSERT(miller_index.size() == intensity.size());
      DIALS_ASSERT(miller_index.size() == variance.size());
      DIALS_ASSERT(miller_index.size() == scale.size());
      DIALS_ASSERT(miller_index.size() == flags.size());
      DIALS_ASSERT(miller_index.size() == xyzcal.size());
      DIALS_ASSERT(miller_index.size() == xyzobs.size());
      DIALS_ASSERT(reflections_.size() == 0);

      reflections_.resize(miller_index.size());
      reflections_["id"] = af::shared<std::size_t>(miller_index.size(), 0);
      reflections_["panel"] = af::shared<std::size_t>(miller_index.size(), 0);
      reflections_["miller_index"] = miller_index;
      reflections_["intensity.sum.value"] = intensity;
      reflections_["intensity.sum.variance"] = variance;
      reflections_["intensity.scale"] = scale;
      reflections_["flags"] = flags;
      reflections_["xyzcal.px"] = xyzcal;
      reflections_["xyzobs.px.value"] = xyzobs;
    }

    vec3<double> compute_centroid(
        af::const_ref<double> data,
        af::const_ref<int> mask,
        af::const_ref< vec3<double> > points) {
      af::shared<double> v;
      af::shared< vec3<double> > p;
      int mask_code = Valid | Foreground;
      for (std::size_t i = 0; i < mask.size(); ++i) {
        if ((mask[i] & mask_code) == mask_code) {
          v.push_back(data[i]);
          p.push_back(points[i]);
        }
      }

      CentroidPoints<double, vec3<double> > centroid(v.const_ref(), p.const_ref());
      return centroid.mean();
    }

    double compute_background(
        af::const_ref<double> data,
        af::const_ref<int> mask) const {

      std::size_t min_pixels = 10;
      double tuning_constant = 1.345;
      std::size_t max_iter = 100;

      // Compute number of background pixels
      std::size_t num_background = 0;
      int mask_code = Valid | Background;
      for (std::size_t i = 0; i < mask.size(); ++i) {
        if ((mask[i] & mask_code) == mask_code) {
          num_background++;
        }
      }
      DIALS_ASSERT(num_background >= min_pixels);

      // Allocate some arrays
      af::shared<double> Y(num_background, 0);
      std::size_t j = 0;
      for (std::size_t i = 0; i < mask.size(); ++i) {
        if ((mask[i] & mask_code) == mask_code) {
          DIALS_ASSERT(j < Y.size());
          DIALS_ASSERT(data[i] >= 0);
          Y[j++] = data[i];
        }
      }
      DIALS_ASSERT(j == Y.size());

      // Compute the median value for the starting value
      double median = detail::median(Y.const_ref());
      if (median == 0) {
        median = 1.0;
      }

      // Compute the result
      RobustPoissonMean result(
          Y.const_ref(),
          median,
          tuning_constant,
          1e-3,
          max_iter);
      DIALS_ASSERT(result.converged());

      // Compute the background
      return result.mean();
    }

    vec2<double> compute_intensity(
        af::const_ref<double> data,
        af::const_ref<int> mask,
        double B) const {

      std::size_t min_pixels = 1;

      af::shared<double> background(data.size());
      for (std::size_t i = 0; i < background.size(); ++i) {
        background[i] = B;
      }

      Summation<double> summation(data, background.const_ref(), mask);


      DIALS_ASSERT(summation.n_signal() >= min_pixels);

      return vec2<double>(summation.intensity(), summation.variance());
    }

    double compute_scale(
        af::const_ref<double> pred,
        af::const_ref<int> mask) const {
      double result = 0;
      int mask_code = Valid | Foreground;
      for (std::size_t i = 0; i < pred.size(); ++i) {
        if ((mask[i] & mask_code) == mask_code) {
          result += pred[i];
        }
      }
      return result;
    }

    af::versa< double, af::c_grid<2> > model() const {
      return model_;
    }

    af::versa< int, af::c_grid<2> > mask() const {
      return mask_;
    }

    af::versa< int, af::c_grid<2> > labels() const {
      return labels_;
    }

    af::shared< cctbx::miller::index<> > indices() const {
      return indices_;
    }

    af::reflection_table reflections() const {
      return reflections_;
    }

  protected:

    af::versa< double, af::c_grid<2> > model_;
    af::versa< int, af::c_grid<2> > mask_;
    af::versa< int, af::c_grid<2> > labels_;
    af::shared< cctbx::miller::index<> > indices_;
    af::reflection_table reflections_;
  };


  class ModelOld {
  public:

    ModelOld(
          const Beam &beam,
          const Detector &detector,
          const Crystal &crystal,
          double mosaicity,
          std::size_t num_sample,
          const af::const_ref<int,af::c_grid<2> > &input_data,
          const af::const_ref<bool,af::c_grid<2> > &input_mask)
      : image_data_(input_data.accessor()),
        image_mask_(input_mask.accessor()),
        image_pred_(input_data.accessor()) {

      // Check the input
      DIALS_ASSERT(mosaicity > 0);
      DIALS_ASSERT(num_sample > 0);

      // Copy the mask and data arrays
      std::copy(input_data.begin(), input_data.end(), image_data_.begin());
      std::copy(input_mask.begin(), input_mask.end(), image_mask_.begin());

      // Generate the image model
      generate_image_model(beam, detector, crystal, mosaicity, num_sample);

      // integrate the data
      integrate();
    }

    void generate_image_model(
          const Beam &beam,
          const Detector &detector,
          const Crystal &crystal,
          double mosaicity,
          std::size_t num_sample) {

      std::size_t xsize = detector[0].get_image_size()[0];
      std::size_t ysize = detector[0].get_image_size()[1];

      // Initialise the transform
      PixelToMillerIndex transform(beam, detector, crystal);

      // Compute the normal constant
      double K = std::pow(1.0 / (std::sqrt(2*pi) * mosaicity), 3);

      for (std::size_t j = 0; j < ysize; ++j) {
        for (std::size_t i = 0; i < xsize; ++i) {

          // Map the corners of the pixel
          vec3<double> A = transform.h(0, i, j);
          vec3<double> B = transform.h(0, i+1, j);
          vec3<double> C = transform.h(0, i, j+1);
          vec3<double> D = transform.h(0, i+1, j+1);
          vec3<double> E = transform.h(0, i+0.5, j+0.5);

          // The integer miller index
          vec3<double> h0(
              std::floor(E[0] + 0.5),
              std::floor(E[1] + 0.5),
              std::floor(E[2] + 0.5));

          // The distance from all corners
          double distance_E = (E - h0).length();

          // Compute the approximate area of the pixel in reciprocal space
          double area_abc = 0.5 * (B-A).cross(C-A).length();
          double area_bcd = 0.5 * (B-D).cross(C-D).length();
          double area = area_abc + area_bcd;

          // Integrate the 3D normal distribution over the pixel
          double sum_f = 0.0;
          for (std::size_t jj = 0; jj < num_sample; ++jj) {
            for (std::size_t ii = 0; ii < num_sample; ++ii) {
              vec3<double> h = transform.h(0,
                  i + ii / (double)num_sample,
                  j + jj / (double)num_sample);

              double distance = (h - h0).length();

              sum_f += K * std::exp(-0.5 * std::pow(distance / mosaicity, 2));
            }
          }

          // Assign pixel as foreground or background depending on mosaicity
          if (distance_E < 3 * mosaicity) {
            image_mask_(j,i) |= Foreground;
          } else if (distance_E < 6 * mosaicity) {
            image_mask_(j,i) |= Background;
          }

          // Compute the predicted integrated intensity on the pixel
          image_pred_(j,i) = area * sum_f / (double)(num_sample * num_sample);

          // Add a pixel to the lookup
          lookup_[h0].push_back(i + j * xsize);
        }
      }
    }

    void integrate() {

      typedef LookupMap::iterator iterator;

      std::size_t xsize = image_data_.accessor()[1];
      std::size_t ysize = image_data_.accessor()[0];

      for (iterator it = lookup_.begin(); it != lookup_.end(); ++it) {

        cctbx::miller::index<> h = it->first;
        const std::vector<std::size_t> &pixels = it->second;

        // Extract pixels
        af::shared<double> reflection_data(pixels.size());
        af::shared<double> reflection_pred(pixels.size());
        af::shared<int> reflection_mask(pixels.size());
        af::shared< vec3<double> > reflection_pnts(pixels.size());
        for (std::size_t i = 0; i < pixels.size(); ++i) {
          int index = pixels[i];
          double x = (index % xsize) + 0.5;
          double y = ((int)(index / xsize)) + 0.5;
          reflection_data[i] = image_data_[pixels[i]];
          reflection_mask[i] = image_data_[pixels[i]];
          reflection_pred[i] = image_pred_[pixels[i]];
          reflection_pnts[i] = vec3<double>(x, y, 0);
        }

        reflection_data_.push_back(reflection_data);
        reflection_mask_.push_back(reflection_mask);
        reflection_pred_.push_back(reflection_pred);

        try {

          // Compute the background
          double B = compute_background(
              reflection_data.const_ref(),
              reflection_mask.const_ref());

          // Compute the intensity
          vec2<double> IV = compute_intensity(
              reflection_data.const_ref(),
              reflection_mask.const_ref(),
              B);

          double I = IV[0];
          double V = IV[1];
          double S = compute_scale(
              reflection_pred.const_ref(),
              reflection_mask.const_ref());

          vec3<double> xcal = compute_centroid(
              reflection_pred.const_ref(),
              reflection_mask.const_ref(),
              reflection_pnts.const_ref());

          vec3<double> xobs = compute_centroid(
              reflection_data.const_ref(),
              reflection_mask.const_ref(),
              reflection_pnts.const_ref());

          background_.push_back(B);
          intensity_.push_back(I);
          variance_.push_back(V);
          scale_.push_back(S);
          success_.push_back(true);
          //std::cout << "Success" << std::endl;
        } catch (dials::error e) {
          background_.push_back(-1);
          intensity_.push_back(-1);
          variance_.push_back(-1);
          scale_.push_back(-1);
          success_.push_back(false);
          //std::cout << e.what() << std::endl;
          continue;
        }
      }
    }

    vec3<double> compute_centroid(
        af::const_ref<double> data,
        af::const_ref<int> mask,
        af::const_ref< vec3<double> > points) {
      af::shared<double> v;
      af::shared< vec3<double> > p;
      int mask_code = Valid | Foreground;
      for (std::size_t i = 0; i < mask.size(); ++i) {
        if ((mask[i] & mask_code) == mask_code) {
          v.push_back(data[i]);
          p.push_back(points[i]);
        }
      }

      CentroidPoints<double, vec3<double> > centroid(v.const_ref(), p.const_ref());
      return centroid.mean();
    }

    double compute_background(
        af::const_ref<double> data,
        af::const_ref<int> mask) const {

      std::size_t min_pixels = 10;
      double tuning_constant = 1.345;
      std::size_t max_iter = 100;

      // Compute number of background pixels
      std::size_t num_background = 0;
      int mask_code = Valid | Background;
      for (std::size_t i = 0; i < mask.size(); ++i) {
        if ((mask[i] & mask_code) == mask_code) {
          num_background++;
        }
      }
      DIALS_ASSERT(num_background >= min_pixels);

      // Allocate some arrays
      af::shared<double> Y(num_background, 0);
      std::size_t j = 0;
      for (std::size_t i = 0; i < mask.size(); ++i) {
        if ((mask[i] & mask_code) == mask_code) {
          DIALS_ASSERT(j < Y.size());
          DIALS_ASSERT(data[i] >= 0);
          Y[j++] = data[i];
        }
      }
      DIALS_ASSERT(j == Y.size());

      // Compute the median value for the starting value
      double median = detail::median(Y.const_ref());
      if (median == 0) {
        median = 1.0;
      }

      // Compute the result
      RobustPoissonMean result(
          Y.const_ref(),
          median,
          tuning_constant,
          1e-3,
          max_iter);
      DIALS_ASSERT(result.converged());

      // Compute the background
      return result.mean();
    }

    vec2<double> compute_intensity(
        af::const_ref<double> data,
        af::const_ref<int> mask,
        double B) const {

      std::size_t min_pixels = 1;

      af::shared<double> background(data.size());
      for (std::size_t i = 0; i < background.size(); ++i) {
        background[i] = B;
      }

      Summation<double> summation(data, background.const_ref(), mask);


      DIALS_ASSERT(summation.n_signal() >= min_pixels);

      return vec2<double>(summation.intensity(), summation.variance());
    }

    double compute_scale(
        af::const_ref<double> pred,
        af::const_ref<int> mask) const {
      double result = 0;
      int mask_code = Valid | Foreground;
      for (std::size_t i = 0; i < pred.size(); ++i) {
        if ((mask[i] & mask_code) == mask_code) {
          result += pred[i];
        }
      }
      return result;
    }

    af::versa< double, af::c_grid<2> > image_pred() const {
      return image_pred_;
    }

    af::versa< double, af::c_grid<2> > image_data() const {
      return image_data_;
    }

    af::versa< int, af::c_grid<2> > image_mask() const {
      return image_mask_;
    }

    af::shared< double > pred(std::size_t index) const {
      DIALS_ASSERT(index < reflection_pred_.size());
      return reflection_pred_[index];
    }

    af::shared< double > data(std::size_t index) const {
      DIALS_ASSERT(index < reflection_data_.size());
      return reflection_data_[index];
    }

    af::shared< int > mask(std::size_t index) const {
      DIALS_ASSERT(index < reflection_mask_.size());
      return reflection_mask_[index];
    }

    double background(std::size_t index) const {
      DIALS_ASSERT(index < background_.size());
      return background_[index];
    }

    double intensity(std::size_t index) const {
      DIALS_ASSERT(index < intensity_.size());
      return intensity_[index];
    }

    double variance(std::size_t index) const {
      DIALS_ASSERT(index < variance_.size());
      return variance_[index];
    }

    double scale(std::size_t index) const {
      DIALS_ASSERT(index < scale_.size());
      return scale_[index];
    }

    double success(std::size_t index) const {
      DIALS_ASSERT(index < success_.size());
      return success_[index];
    }

    std::size_t size() const {
      return intensity_.size();
    }

  protected:

    typedef std::map< cctbx::miller::index<>, std::vector<std::size_t> > LookupMap;
    LookupMap lookup_;

    af::versa< double, af::c_grid<2> > image_data_;
    af::versa< int,    af::c_grid<2> > image_mask_;
    af::versa< double, af::c_grid<2> > image_pred_;

    af::shared< af::shared< double > > reflection_data_;
    af::shared< af::shared< int    > > reflection_mask_;
    af::shared< af::shared< double > > reflection_pred_;

    af::shared< double > intensity_;
    af::shared< double > variance_;
    af::shared< double > background_;
    af::shared< double > scale_;
    af::shared< bool   > success_;
  };


  class Model {
  public:

    Model(
          const Beam &beam,
          const Detector &detector,
          const Crystal &crystal,
          double mosaicity,
          std::size_t num_sample,
          const af::const_ref<int,af::c_grid<2> > &input_data,
          const af::const_ref<bool,af::c_grid<2> > &input_mask,
          af::reflection_table reflections)
      : image_data_(input_data.accessor()),
        image_mask_(input_mask.accessor()),
        image_pred_(input_data.accessor()) {

      // Check the input
      DIALS_ASSERT(mosaicity > 0);
      DIALS_ASSERT(num_sample > 0);

      // Copy the mask and data arrays
      std::copy(input_data.begin(), input_data.end(), image_data_.begin());
      std::copy(input_mask.begin(), input_mask.end(), image_mask_.begin());

      // Generate the pixel index lookup
      generate_pixel_index_lookup(beam, detector, crystal);

      // Generate reflection pixel model
      generate_image_model(beam, detector, crystal, reflections, mosaicity, num_sample);

      // Integrate reflections
      integrate_reflections();

      /* // Generate the image model */

      /* // integrate the data */
      /* integrate(); */
    }

    void generate_pixel_index_lookup(

          const Beam &beam,
          const Detector &detector,
          const Crystal &crystal) {

      std::size_t xsize = detector[0].get_image_size()[0];
      std::size_t ysize = detector[0].get_image_size()[1];

      // Initialise the transform
      PixelToMillerIndex transform(beam, detector, crystal);

      for (std::size_t j = 0; j < ysize; ++j) {
        for (std::size_t i = 0; i < xsize; ++i) {

          // Map the corners of the pixel
          vec3<double> h = transform.h(0, i+0.5, j+0.5);

          // The integer miller index
          vec3<double> h0(
              std::floor(h[0] + 0.5),
              std::floor(h[1] + 0.5),
              std::floor(h[2] + 0.5));

          // Add a pixel to the lookup
          lookup_[h0].push_back(i + j * xsize);
        }
      }
    }

    void generate_image_model(
          const Beam &beam,
          const Detector &detector,
          const Crystal &crystal,
          af::reflection_table reflections,
          double mosaicity,
          std::size_t num_sample) {

      af::const_ref< cctbx::miller::index<> > miller_index = reflections["miller_index"];

      std::size_t xsize = detector[0].get_image_size()[0];
      std::size_t ysize = detector[0].get_image_size()[1];

      // Initialise the transform
      PixelToMillerIndex transform(beam, detector, crystal);

      // Compute the normal constant
      double K = std::pow(1.0 / (std::sqrt(2*pi) * mosaicity), 3);

      for (std::size_t i = 0; i < miller_index.size(); ++i) {

        vec3<double> h0(
            miller_index[i][0],
            miller_index[i][1],
            miller_index[i][2]);
        const std::vector<std::size_t> &pixels = lookup_[h0];

        // Extract pixels
        af::shared<double> reflection_data(pixels.size());
        af::shared<double> reflection_pred(pixels.size());
        af::shared<int> reflection_mask(pixels.size());
        af::shared< vec3<double> > reflection_pnts(pixels.size());
        for (std::size_t k = 0; k < pixels.size(); ++k) {
          int index = pixels[k];
          double i = (index % xsize) + 0.5;
          double j = ((int)(index / xsize)) + 0.5;

          // Map the corners of the pixel
          vec3<double> A = transform.h(0, i, j);
          vec3<double> B = transform.h(0, i+1, j);
          vec3<double> C = transform.h(0, i, j+1);
          vec3<double> D = transform.h(0, i+1, j+1);
          vec3<double> E = transform.h(0, i+0.5, j+0.5);

          // The distance from all corners
          double distance_E = (E - h0).length();

          // Compute the approximate area of the pixel in reciprocal space
          double area_abc = 0.5 * (B-A).cross(C-A).length();
          double area_bcd = 0.5 * (B-D).cross(C-D).length();
          double area = area_abc + area_bcd;

          // Integrate the 3D normal distribution over the pixel
          double sum_f = 0.0;
          for (std::size_t jj = 0; jj < num_sample; ++jj) {
            for (std::size_t ii = 0; ii < num_sample; ++ii) {
              vec3<double> h = transform.h(0,
                  i + ii / (double)num_sample,
                  j + jj / (double)num_sample);

              double distance = (h - h0).length();

              sum_f += K * std::exp(-0.5 * std::pow(distance / mosaicity, 2));
            }
          }

          reflection_pnts[k] = vec3<double>(i, j, 0);

          // Set pixel data
          reflection_data[k] = image_data_[index];

          // Assign pixel as foreground or background depending on mosaicity
          if (distance_E < 0.3) {//3 * mosaicity) {
            reflection_mask[k] = image_mask_[index] | Foreground;
          } else if (distance_E < 0.5) {// * mosaicity) {
            reflection_mask[k] = image_mask_[index] | Background;
          }

          // Compute the predicted integrated intensity on the pixel
          reflection_pred[k] = area * sum_f / (double)(num_sample * num_sample);
        }

        // Add to arrays
        reflection_data_.push_back(reflection_data);
        reflection_mask_.push_back(reflection_mask);
        reflection_pred_.push_back(reflection_pred);
        reflection_pnts_.push_back(reflection_pnts);
      }
    }

    void integrate_reflections() {

      typedef LookupMap::iterator iterator;

      std::size_t xsize = image_data_.accessor()[1];
      std::size_t ysize = image_data_.accessor()[0];

      for (std::size_t i = 0; i < reflection_data_.size(); ++i) {

        af::shared<double> reflection_data = reflection_data_[i];
        af::shared<double> reflection_pred = reflection_pred_[i];
        af::shared<int> reflection_mask = reflection_mask_[i];
        af::shared< vec3<double> > reflection_pnts = reflection_pnts_[i];

        try {

          // Compute the background
          double B = compute_background(
              reflection_data.const_ref(),
              reflection_mask.const_ref());

          // Compute the intensity
          vec2<double> IV = compute_intensity(
              reflection_data.const_ref(),
              reflection_mask.const_ref(),
              B);

          double I = IV[0];
          double V = IV[1];
          double S = compute_scale(
              reflection_pred.const_ref(),
              reflection_mask.const_ref());

          vec3<double> xcal = compute_centroid(
              reflection_pred.const_ref(),
              reflection_mask.const_ref(),
              reflection_pnts.const_ref());

          vec3<double> xobs = compute_centroid(
              reflection_data.const_ref(),
              reflection_mask.const_ref(),
              reflection_pnts.const_ref());

          observed_.push_back(xobs);
          predicted_.push_back(xcal);
          background_.push_back(B);
          intensity_.push_back(I);
          variance_.push_back(V);
          scale_.push_back(S);
          success_.push_back(true);
          //std::cout << "Success" << std::endl;
        } catch (dials::error e) {
          observed_.push_back(vec3<double>(-1,-1,-1));
          predicted_.push_back(vec3<double>(-1,-1,-1));
          background_.push_back(-1);
          intensity_.push_back(-1);
          variance_.push_back(-1);
          scale_.push_back(-1);
          success_.push_back(false);
          //std::cout << e.what() << std::endl;
          continue;
        }
      }
    }

    vec3<double> compute_centroid(
        af::const_ref<double> data,
        af::const_ref<int> mask,
        af::const_ref< vec3<double> > points) {
      af::shared<double> v;
      af::shared< vec3<double> > p;
      int mask_code = Valid | Foreground;
      for (std::size_t i = 0; i < mask.size(); ++i) {
        if ((mask[i] & mask_code) == mask_code) {
          v.push_back(data[i]);
          p.push_back(points[i]);
        }
      }

      CentroidPoints<double, vec3<double> > centroid(v.const_ref(), p.const_ref());
      return centroid.mean();
    }

    double compute_background(
        af::const_ref<double> data,
        af::const_ref<int> mask) const {

      std::size_t min_pixels = 10;
      double tuning_constant = 1.345;
      std::size_t max_iter = 100;

      // Compute number of background pixels
      std::size_t num_background = 0;
      int mask_code = Valid | Background;
      for (std::size_t i = 0; i < mask.size(); ++i) {
        if ((mask[i] & mask_code) == mask_code) {
          num_background++;
        }
      }
      DIALS_ASSERT(num_background >= min_pixels);

      // Allocate some arrays
      af::shared<double> Y(num_background, 0);
      std::size_t j = 0;
      for (std::size_t i = 0; i < mask.size(); ++i) {
        if ((mask[i] & mask_code) == mask_code) {
          DIALS_ASSERT(j < Y.size());
          DIALS_ASSERT(data[i] >= 0);
          Y[j++] = data[i];
        }
      }
      DIALS_ASSERT(j == Y.size());

      // Compute the median value for the starting value
      double median = detail::median(Y.const_ref());
      if (median == 0) {
        median = 1.0;
      }

      // Compute the result
      RobustPoissonMean result(
          Y.const_ref(),
          median,
          tuning_constant,
          1e-3,
          max_iter);
      DIALS_ASSERT(result.converged());

      // Compute the background
      return result.mean();
    }

    vec2<double> compute_intensity(
        af::const_ref<double> data,
        af::const_ref<int> mask,
        double B) const {

      std::size_t min_pixels = 1;

      af::shared<double> background(data.size());
      for (std::size_t i = 0; i < background.size(); ++i) {
        background[i] = B;
      }

      Summation<double> summation(data, background.const_ref(), mask);


      DIALS_ASSERT(summation.n_signal() >= min_pixels);

      return vec2<double>(summation.intensity(), summation.variance());
    }

    double compute_scale(
        af::const_ref<double> pred,
        af::const_ref<int> mask) const {
      double result = 0;
      int mask_code = Valid | Foreground;
      for (std::size_t i = 0; i < pred.size(); ++i) {
        if ((mask[i] & mask_code) == mask_code) {
          result += pred[i];
        }
      }
      return result;
    }

    af::versa< double, af::c_grid<2> > image_pred() const {
      return image_pred_;
    }

    af::versa< double, af::c_grid<2> > image_data() const {
      return image_data_;
    }

    af::versa< int, af::c_grid<2> > image_mask() const {
      return image_mask_;
    }

    af::shared< double > pred(std::size_t index) const {
      DIALS_ASSERT(index < reflection_pred_.size());
      return reflection_pred_[index];
    }

    af::shared< double > data(std::size_t index) const {
      DIALS_ASSERT(index < reflection_data_.size());
      return reflection_data_[index];
    }

    af::shared< int > mask(std::size_t index) const {
      DIALS_ASSERT(index < reflection_mask_.size());
      return reflection_mask_[index];
    }

    vec3<double> observed(std::size_t index) const {
      DIALS_ASSERT(index < observed_.size());
      return observed_[index];
    }

    vec3<double> predicted(std::size_t index) const {
      DIALS_ASSERT(index < predicted_.size());
      return predicted_[index];
    }

    double background(std::size_t index) const {
      DIALS_ASSERT(index < background_.size());
      return background_[index];
    }

    double intensity(std::size_t index) const {
      DIALS_ASSERT(index < intensity_.size());
      return intensity_[index];
    }

    double variance(std::size_t index) const {
      DIALS_ASSERT(index < variance_.size());
      return variance_[index];
    }

    double scale(std::size_t index) const {
      DIALS_ASSERT(index < scale_.size());
      return scale_[index];
    }

    double success(std::size_t index) const {
      DIALS_ASSERT(index < success_.size());
      return success_[index];
    }

    std::size_t size() const {
      return intensity_.size();
    }

  protected:

    typedef std::map< cctbx::miller::index<>, std::vector<std::size_t> > LookupMap;
    LookupMap lookup_;

    af::versa< double, af::c_grid<2> > image_data_;
    af::versa< int,    af::c_grid<2> > image_mask_;
    af::versa< double, af::c_grid<2> > image_pred_;

    af::shared< af::shared< double > > reflection_data_;
    af::shared< af::shared< int    > > reflection_mask_;
    af::shared< af::shared< double > > reflection_pred_;
    af::shared< af::shared< vec3<double> > > reflection_pnts_;

    af::shared< vec3<double> > observed_;
    af::shared< vec3<double> > predicted_;
    af::shared< double > intensity_;
    af::shared< double > variance_;
    af::shared< double > background_;
    af::shared< double > scale_;
    af::shared< bool   > success_;
  };


  BOOST_PYTHON_MODULE(dials_scratch_jmp_stills_ext)
  {
    class_<Model>("Model", no_init)
      .def(init<
        const Beam &,
        const Detector &,
        const Crystal &,
        double,
        std::size_t,
        const af::const_ref<int,  af::c_grid<2> >,
        const af::const_ref<bool, af::c_grid<2> >&,
        af::reflection_table>())
      .def("image_data", &Model::image_data)
      .def("image_mask", &Model::image_mask)
      .def("image_pred", &Model::image_pred)
      .def("data", &Model::data)
      .def("mask", &Model::mask)
      .def("pred", &Model::pred)
      .def("observed", &Model::observed)
      .def("predicted", &Model::predicted)
      .def("background", &Model::background)
      .def("intensity", &Model::intensity)
      .def("variance", &Model::variance)
      .def("scale", &Model::scale)
      .def("success", &Model::success)
      .def("__len__", &Model::size)
      ;
  }

}}} // namespace dials_scratch::examples::boost_python
