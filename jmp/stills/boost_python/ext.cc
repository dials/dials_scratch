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
#include <dials/model/data/shoebox.h>
#include <dials/error.h>

using namespace boost::python;

namespace dials { namespace algorithms { namespace boost_python {

  using scitbx::constants::pi;
  using dials::model::Valid;
  using dials::model::Background;
  using dials::model::Foreground;
  using dials::model::Shoebox;
  using dxtbx::model::Beam;
  using dxtbx::model::Detector;
  using dxtbx::model::Crystal;

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


  class Model {
  public:

    typedef std::map< cctbx::miller::index<>, std::vector<std::size_t> > LookupMap;

    Model(
          const Beam &beam,
          const Detector &detector,
          const Crystal &crystal,
          af::reflection_table reflections,
          const af::const_ref<int,af::c_grid<2> > &input_data,
          const af::const_ref<bool,af::c_grid<2> > &input_mask,
          double mosaicity,
          double bandpass,
          double foreground_limit,
          double background_limit,
          std::size_t num_samples,
          bool predict_all,
          bool use_mosaicity_for_mask)
      : beam_(beam),
        detector_(detector),
        crystal_(crystal),
        reflections_(reflections),
        image_data_(input_data.accessor()),
        image_mask_(input_mask.accessor()),
        mosaicity_(mosaicity),
        bandpass_(bandpass),
        foreground_limit_(foreground_limit),
        background_limit_(background_limit),
        num_samples_(num_samples),
        predict_all_(predict_all),
        use_mosaicity_for_mask_(use_mosaicity_for_mask),
        update_pixel_lookup_(true),
        update_image_model_(true),
        update_reflection_data_(true),
        first_(true),
        wavelength_values_(5),
        wavelength_weights_(5) {

      // Check the input
      if (predict_all) {
        DIALS_ASSERT(reflections.size() == 0);
      } else {
        DIALS_ASSERT(reflections.size() > 0);
      }

      DIALS_ASSERT(mosaicity > 0);
      DIALS_ASSERT(bandpass >= 0);
      DIALS_ASSERT(foreground_limit_ > 0);
      DIALS_ASSERT(background_limit_ > foreground_limit_);
      DIALS_ASSERT(num_samples > 0);
      DIALS_ASSERT(input_data.accessor().all_eq(input_mask.accessor()));
      DIALS_ASSERT(detector.size() == 1);
      DIALS_ASSERT(detector[0].get_image_size()[0] == input_data.accessor()[1]);
      DIALS_ASSERT(detector[0].get_image_size()[1] == input_data.accessor()[0]);

      // Copy the mask and data arrays
      std::copy(input_data.begin(), input_data.end(), image_data_.begin());
      std::copy(input_mask.begin(), input_mask.end(), image_mask_.begin());

      // We use guass-hermite polynomials to approximate a normal distribution
      // of wavelengths. We use the 5th degree polynomial and hard code the
      // roots below. The actual wavelength values are given by:
      // sqrt(2) * sigma * x + mu
      //
      // The weights are then
      //
      // 2^(n-1) n! sqrt(pi) / (n^2 H_(n-1)(x)^2)
      wavelength_values_[0] = -std::sqrt(5.0/2.0 + std::sqrt(5.0/2.0));
      wavelength_values_[1] = -std::sqrt(5.0/2.0 - std::sqrt(5.0/2.0));
      wavelength_values_[2] = 0;
      wavelength_values_[3] = std::sqrt(5.0/2.0 - std::sqrt(5.0/2.0));
      wavelength_values_[4] = std::sqrt(5.0/2.0 + std::sqrt(5.0/2.0));
      for (std::size_t i = 0; i < 5; ++i) {
        double x = wavelength_values_[i];
        double H4 = 16.0*x*x*x*x - 48.0*x*x + 12.0;
        wavelength_weights_[i] = 16 * 120 * std::sqrt(pi) / (25 * H4);
      }

    }

    Beam get_beam() const {
      return beam_;
    }

    void set_beam(const Beam &beam) {
      beam_ = beam;
      update_pixel_lookup_ = true;
      update_image_model_ = true;
      update_reflection_data_ = true;
    }

    Detector get_detector() const {
      return detector_;
    }

    void set_detector(const Detector &detector) {
      DIALS_ASSERT(detector[0].get_image_size().all_eq(detector_[0].get_image_size()));
      detector_ = detector;
      update_pixel_lookup_ = true;
      update_image_model_ = true;
      update_reflection_data_ = true;
    }

    Crystal get_crystal() const {
      return crystal_;
    }

    void set_crystal(const Crystal &crystal) {
      crystal_ = crystal;
      update_pixel_lookup_ = true;
      update_image_model_ = true;
      update_reflection_data_ = true;
    }

    af::reflection_table get_reflections() const {
      return reflections_;
    }

    void set_reflections(const af::reflection_table &reflections) {
      reflections_ = reflections;
      update_image_model_ = true;
      update_reflection_data_ = true;
    }

    af::versa< int, af::c_grid<2> > get_image_data() const {
      return image_data_;
    }

    void set_image_data(const af::const_ref<int, af::c_grid<2> > &image_data) {
      DIALS_ASSERT(image_data.accessor().all_eq(image_data_.accessor()));
      std::copy(image_data.begin(), image_data.end(), image_data_.begin());
      update_image_model_ = true;
      update_reflection_data_ = true;
    }

    af::versa< int, af::c_grid<2> > get_image_mask() const {
      return image_mask_;
    }

    void set_image_mask(const af::const_ref<int, af::c_grid<2> > &image_mask) {
      DIALS_ASSERT(image_mask.accessor().all_eq(image_mask_.accessor()));
      std::copy(image_mask.begin(), image_mask.end(), image_mask_.begin());
      update_image_model_ = true;
      update_reflection_data_ = true;
    }

    af::versa< double, af::c_grid<2> > get_image_pred() const {
      af::versa< double, af::c_grid<2> > result(image_mask_.accessor());
      for (std::size_t i = 0; i < reflection_pred_.size(); ++i) {
        af::const_ref< double > pred = reflection_pred_[i].const_ref();
        af::const_ref< vec3<int> > pnts = reflection_pnts_[i].const_ref();
        for (std::size_t j = 0; j < pred.size(); ++j) {
          int x = pnts[j][0];
          int y = pnts[j][1];
          result(y,x) = pred[j];
        }
      }
      return result;
    }

    double get_mosaicity() const {
      return mosaicity_;
    }

    void set_mosaicity(const double &mosaicity) {
      mosaicity_ = mosaicity;
      update_image_model_ = true;
      update_reflection_data_ = true;
    }

    double get_bandpass() const {
      return bandpass_;
    }

    void set_bandpass(const double &bandpass) {
      bandpass_ = bandpass;
      update_image_model_ = true;
      update_reflection_data_ = true;
    }

    double get_foreground_limit() const {
      return foreground_limit_;
    }

    void set_foreground_limit(const double &foreground_limit) {
      foreground_limit_ = foreground_limit;
      update_image_model_ = true;
      update_reflection_data_ = true;
    }

    double get_background_limit() const {
      return background_limit_;
    }

    void set_background_limit(const double &background_limit) {
      background_limit_ = background_limit;
      update_image_model_ = true;
      update_reflection_data_ = true;
    }

    double get_num_samples() const {
      return num_samples_;
    }

    void set_num_samples(const double &num_samples) {
      num_samples_ = num_samples;
      update_image_model_ = true;
      update_reflection_data_ = true;
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

    af::shared< vec3<double> > observed() const {
      return observed_;
    }

    af::shared< vec3<double> > predicted() const {
      return predicted_;
    }

    af::shared<double> background() const {
      return background_;
    }

    af::shared<double> intensity() const {
      return intensity_;
    }

    af::shared<double> variance() const {
      return variance_;
    }

    af::shared<std::size_t> num_foreground() const {
      return num_foreground_;
    }

    af::shared<std::size_t> num_background() const {
      return num_background_;
    }

    af::shared<double> scale() const {
      return scale_;
    }

    af::shared<bool> success() const {
      return success_;
    }

    af::shared< Shoebox<> > shoebox() const {

      af::shared< Shoebox<> > result(reflections_.size());

      for (std::size_t l = 0; l < reflection_mask_.size(); ++l) {

        // Get the reflection pixel information
        af::const_ref< vec3<int> > pnts = reflection_pnts_[l].const_ref();
        af::const_ref< double >    data = reflection_data_[l].const_ref();
        af::const_ref< int >       mask = reflection_mask_[l].const_ref();

        // Compute the bounding box
        DIALS_ASSERT(pnts.size() > 0);
        int x0 = pnts[0][0];
        int x1 = pnts[0][0]+1;
        int y0 = pnts[0][1];
        int y1 = pnts[0][1]+1;
        for (std::size_t i = 1; i < pnts.size(); ++i) {
          if (pnts[i][0] < x0)  x0 = pnts[i][0];
          if (pnts[i][0] >= x1) x1 = pnts[i][0] + 1;
          if (pnts[i][1] < y0)  y0 = pnts[i][1];
          if (pnts[i][1] >= y1) y1 = pnts[i][1] + 1;
        }
        DIALS_ASSERT(x0 < x1);
        DIALS_ASSERT(y0 < y1);
        int6 bbox(x0, x1, y0, y1, 0, 1);

        // Allocate the shoebox
        Shoebox<> sbox(0, bbox);
        sbox.allocate();

        // Fill the shobox
        for (std::size_t i = 1; i < pnts.size(); ++i) {
          int x = pnts[i][0];
          int y = pnts[i][1];
          int ii = x - x0;
          int jj = y - y0;
          DIALS_ASSERT(ii >= 0 && ii < sbox.data.accessor()[2]);
          DIALS_ASSERT(jj >= 0 && jj < sbox.data.accessor()[1]);
          sbox.data(0, jj,ii) = data[i];
          sbox.mask(0, jj,ii) = mask[i];
          sbox.background(0, jj,ii) = background_[l];
        }

        result[l] = sbox;
      }

      // Return the shoebox
      return result;
    }

    std::size_t size() const {
      return intensity_.size();
    }

    af::shared<double> wavelength_values() const {
      double wavelength = beam_.get_wavelength();
      af::shared<double> result(wavelength_values_.size());
      for (std::size_t i = 0; i < result.size(); ++i) {
        result[i] = std::sqrt(2) * bandpass_ * wavelength_values_[i] + wavelength;
      }
      return result;
    }

    af::shared<double> wavelength_weights() const {
      return wavelength_weights_;
    }

    void update(bool pixel_lookup = true) {

      // Generate the pixel index lookup
      if (update_pixel_lookup_ && (pixel_lookup || first_)) {
        generate_pixel_index_lookup();
        update_pixel_lookup_ = false;
        update_image_model_ = true;
      }

      // Generate reflection pixel model
      if (update_image_model_) {
        generate_image_model();
        update_image_model_ = false;
        update_reflection_data_ = true;
      }

      // Integrate reflections
      if (update_reflection_data_) {
        integrate_reflections();
        update_reflection_data_ = false;
      }

      // Unset the first flag
      first_ = false;
    }

    double least_squares_score() const {

      int mask_code = Valid | Foreground;

      double score = 0;

      for (std::size_t i = 0; i < intensity_.size(); ++i) {

        double B = background_[i];
        double I = intensity_[i];
        double V = variance_[i];
        double S = scale_[i];
        bool T = success_[i];

        if (T && I > 0 && S > 0 && V > 0) {

          double P = I / S;

          af::const_ref< double > c = reflection_data_[i].const_ref();
          af::const_ref< int>     m = reflection_mask_[i].const_ref();
          af::const_ref< double > p = reflection_pred_[i].const_ref();

          DIALS_ASSERT(c.size() == m.size());
          DIALS_ASSERT(c.size() == p.size());

          for (std::size_t j = 0; j < c.size(); ++j) {
            if ((m[j] & mask_code) == mask_code) {
              score += std::pow(c[j] - (B + P*p[j]), 2);
            }
          }
        }
      }

      return score;
    }

    double maximum_likelihood_score() const {

      const double TINY = 1e-20;

      int mask_code = Valid | Foreground;

      double score = 0;

      for (std::size_t i = 0; i < intensity_.size(); ++i) {

        double B = background_[i];
        double I = intensity_[i];
        double V = variance_[i];
        double S = scale_[i];
        bool T = success_[i];

        if (T && I > 0 && S > 0 && V > 0) {

          double P = I / S;

          af::const_ref< double > c = reflection_data_[i].const_ref();
          af::const_ref< int>     m = reflection_mask_[i].const_ref();
          af::const_ref< double > p = reflection_pred_[i].const_ref();

          DIALS_ASSERT(c.size() == m.size());
          DIALS_ASSERT(c.size() == p.size());

          for (std::size_t j = 0; j < c.size(); ++j) {
            if ((m[j] & mask_code) == mask_code) {
              double v = B + P*p[j];
              if (v < TINY) {
                v = TINY;
              }
              score += c[j] * std::log(v) - v;
            }
          }
        }
      }

      return score;
    }

  protected:

    void generate_pixel_index_lookup() {

      // Initialise the transform
      PixelToMillerIndex transform(beam_, detector_, crystal_);

      // Clear the map
      lookup_.clear();

      std::size_t xsize = detector_[0].get_image_size()[0];
      std::size_t ysize = detector_[0].get_image_size()[1];

      for (std::size_t j = 0; j < ysize; ++j) {
        for (std::size_t i = 0; i < xsize; ++i) {

          // Map the centre of the pixel
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

      // If predict all is set then fill the reflection table miller indices
      if (predict_all_) {

        af::shared< cctbx::miller::index<> > miller_indices;
        for (LookupMap::iterator it = lookup_.begin(); it != lookup_.end(); ++it) {
          miller_indices.push_back(it->first);
        }

        // Clear current stuff and add new miller indices
        reflections_.clear();
        reflections_.resize(miller_indices.size());
        reflections_["miller_index"] = miller_indices;
      }
    }

    void generate_image_model() {

      af::const_ref< cctbx::miller::index<> > miller_index = reflections_["miller_index"];

      std::size_t xsize = detector_[0].get_image_size()[0];

      // Initialise the transform
      PixelToMillerIndex transform(beam_, detector_, crystal_);

      // Get the A matrix
      mat3<double> UB = crystal_.get_A();

      // Compute the normal constant
      double K = std::pow(1.0 / (std::sqrt(2*pi) * mosaicity_), 3);
      double K2 = 1.0 / (std::sqrt(2*pi) * mosaicity_);

      // Resize the arrays
      reflection_pnts_ = af::shared< af::shared< vec3<int> > >(miller_index.size());
      reflection_data_ = af::shared< af::shared<double> >(miller_index.size());
      reflection_pred_ = af::shared< af::shared<double> >(miller_index.size());
      reflection_mask_ = af::shared< af::shared<int> >(miller_index.size());
      for (std::size_t l = 0; l < miller_index.size(); ++l) {

        vec3<double> h0(
            miller_index[l][0],
            miller_index[l][1],
            miller_index[l][2]);
        vec3<double> q0 = UB * h0;
        const std::vector<std::size_t> &pixels = lookup_[h0];

        // Extract pixels
        af::shared<double> reflection_data(pixels.size());
        af::shared<double> reflection_pred(pixels.size());
        af::shared<int> reflection_mask(pixels.size());
        af::shared< vec3<int> > reflection_pnts(pixels.size());
        for (std::size_t k = 0; k < pixels.size(); ++k) {
          int index = pixels[k];
          int i = (index % xsize);
          int j = ((int)(index / xsize));

          // Map the corners of the pixel
          vec3<double> A = transform.q(0, i, j);
          vec3<double> B = transform.q(0, i+1, j);
          vec3<double> C = transform.q(0, i, j+1);
          vec3<double> D = transform.q(0, i+1, j+1);
          vec3<double> E = transform.h(0, i+0.5, j+0.5);
          vec3<double> F = transform.q(0, i+0.5, j+0.5);

          // The distance from all corners
          double distance_E = (E - h0).length();
          double distance_F = (F - q0).length();

          // Compute the approximate area of the pixel in reciprocal space
          double area_abc = 0.5 * (B-A).cross(C-A).length();
          double area_bcd = 0.5 * (B-D).cross(C-D).length();
          double area = area_abc + area_bcd;

          // Integrate the 3D normal distribution over the pixel
          double sum_f = 0.0;
          for (std::size_t jj = 0; jj < num_samples_; ++jj) {
            for (std::size_t ii = 0; ii < num_samples_; ++ii) {
              vec3<double> q = transform.q(0,
                  i + (ii + 0.5) / (double)num_samples_,
                  j + (jj + 0.5) / (double)num_samples_);

              double distance = (q - q0).length();

              sum_f += std::exp(-0.5 * std::pow(distance / mosaicity_, 2));
            }
          }

          // Compute the value of the integral
          double I = area * K * sum_f / (double)(num_samples_ * num_samples_);

          // FIXME
          // The partiality is different depending on whether I do it in HKL
          // space or Q space even when the mosaicity computed in either space
          // is equivalent which is clearly wrong. However, if I do this, then
          // I get the same (more sensible) answer. The reason why is not clear
          // to me!
          //
          // Thinking about it, this is probably wrong answer - I bet it has something
          // to do with the ratio between the areas which I'm missing somewhere.
          /* I *= mosaicity_; */
          //
          // Whilst the integral of the normal density over a given range is the
          // same in whatever space it is taken, the the value of the normal
          // density function depends on the space as the variance may be
          // different under linear transform. We therefore normalize the scale
          // with respect to the (0,0,0) reflection by dividing through by a
          // constant term. This is then the reletive likelihood
          I /= K2;

          // Compute the predicted integrated intensity on the pixel
          reflection_pred[k] = I;

          // Add the point
          reflection_pnts[k] = vec3<int>(i, j, 0);

          // Set pixel data
          reflection_data[k] = image_data_[index];

          // Assign pixel as foreground or background depending on mosaicity
          if (use_mosaicity_for_mask_) {
            if (distance_F < 3 * mosaicity_) {
              reflection_mask[k] = image_mask_[index] | Foreground;
            } else if (distance_F < 5 * mosaicity_) {
              reflection_mask[k] = image_mask_[index] | Background;
            }
          } else {
            if (distance_E < foreground_limit_) {
              reflection_mask[k] = image_mask_[index] | Foreground;
            } else if (distance_E < background_limit_) {
              reflection_mask[k] = image_mask_[index] | Background;
            }
          }
        }

        // Add to arrays
        reflection_data_[l] = reflection_data;
        reflection_mask_[l] = reflection_mask;
        reflection_pred_[l] = reflection_pred;
        reflection_pnts_[l] = reflection_pnts;
      }
    }

    void integrate_reflections() {

      typedef LookupMap::iterator iterator;

      // Initialise all the arrays
      observed_ = af::shared< vec3<double> >(reflections_.size());
      predicted_ = af::shared< vec3<double> >(reflections_.size());
      background_ = af::shared< double >(reflections_.size());
      intensity_ = af::shared< double >(reflections_.size());
      variance_ = af::shared< double >(reflections_.size());
      scale_ = af::shared< double >(reflections_.size());
      num_foreground_ = af::shared< std::size_t >(reflections_.size());
      num_background_ = af::shared< std::size_t >(reflections_.size());
      success_ = af::shared< bool >(reflections_.size());

      // The mask codes
      int mask_code_bg = Valid | Background;
      int mask_code_fg = Valid | Foreground;

      // Loop through all the reflections
      for (std::size_t i = 0; i < reflections_.size(); ++i) {

        // Get the pixel data for the reflection
        af::shared< vec3<int> > reflection_pnts = reflection_pnts_[i];
        af::shared<double>      reflection_data = reflection_data_[i];
        af::shared<double>      reflection_pred = reflection_pred_[i];
        af::shared<int>         reflection_mask = reflection_mask_[i];

        // Convert to double
        af::shared< vec3<double> > reflection_pnts_double(reflection_pnts.size());
        for (std::size_t j = 0; j < reflection_pnts.size(); ++j) {
          reflection_pnts_double[j][0] = reflection_pnts[j][0] + 0.5;
          reflection_pnts_double[j][1] = reflection_pnts[j][1] + 0.5;
          reflection_pnts_double[j][2] = reflection_pnts[j][2] + 0.5;
        }

        // Compute the number of foreground and background pixels
        std::size_t count_fg = 0;
        std::size_t count_bg = 0;
        for (std::size_t j = 0; j < reflection_mask.size(); ++j) {
          if ((reflection_mask[j] & mask_code_fg) == mask_code_fg) {
            count_fg++;
          }
          if ((reflection_mask[j] & mask_code_bg) == mask_code_bg) {
            count_bg++;
          }
        }
        num_foreground_[i] = count_fg;
        num_background_[i] = count_bg;

        // Integrate the reflection
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
              reflection_pnts_double.const_ref());

          vec3<double> xobs = compute_centroid(
              reflection_data.const_ref(),
              reflection_mask.const_ref(),
              reflection_pnts_double.const_ref());

          observed_[i] = xobs;
          predicted_[i] = xcal;
          background_[i] = B;
          intensity_[i] = I;
          variance_[i] = V;
          scale_[i] = S;
          success_[i] = true;
          //std::cout << "Success" << std::endl;
        } catch (dials::error e) {
          observed_[i] = vec3<double>(-1,-1,-1);
          predicted_[i] = vec3<double>(-1,-1,-1);
          background_[i] = -1;
          intensity_[i] = -1;
          variance_[i] = -1;
          scale_[i] = -1;
          success_[i] = false;
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

    LookupMap lookup_;

    Beam beam_;
    Detector detector_;
    Crystal crystal_;
    af::reflection_table reflections_;

    af::versa< int, af::c_grid<2> > image_data_;
    af::versa< int, af::c_grid<2> > image_mask_;

    double mosaicity_;
    double bandpass_;
    double foreground_limit_;
    double background_limit_;
    std::size_t num_samples_;
    bool predict_all_;
    bool use_mosaicity_for_mask_;

    bool update_pixel_lookup_;
    bool update_image_model_;
    bool update_reflection_data_;
    bool first_;

    af::shared< double > wavelength_values_;
    af::shared< double > wavelength_weights_;

    af::shared< af::shared< double > > reflection_data_;
    af::shared< af::shared< int    > > reflection_mask_;
    af::shared< af::shared< double > > reflection_pred_;
    af::shared< af::shared< vec3<int> > > reflection_pnts_;

    af::shared< vec3<double> > observed_;
    af::shared< vec3<double> > predicted_;
    af::shared< double > intensity_;
    af::shared< double > variance_;
    af::shared< double > background_;
    af::shared< double > scale_;
    af::shared< std::size_t > num_foreground_;
    af::shared< std::size_t > num_background_;
    af::shared< bool   > success_;
  };




  BOOST_PYTHON_MODULE(dials_scratch_jmp_stills_ext)
  {
    class_<Model>("Model", no_init)
      .def(init<
        const Beam &,
        const Detector &,
        const Crystal &,
        af::reflection_table,
        const af::const_ref<int,  af::c_grid<2> >&,
        const af::const_ref<bool, af::c_grid<2> >&,
        double,
        double,
        double,
        double,
        std::size_t,
        bool,
        bool>((
            arg("beam"),
            arg("detector"),
            arg("crystal"),
            arg("reflections"),
            arg("image_data"),
            arg("image_mask"),
            arg("mosaicity")        = 0.1,
            arg("bandpass")         = 0.01,
            arg("foreground_limit") = 3.0,
            arg("background_limit") = 5.0,
            arg("num_samples")      = 5,
            arg("predict_all")      = false,
            arg("use_mosaicity_for_mask") = false)))
      .add_property("beam",
          &Model::get_beam,
          &Model::set_beam)
      .add_property("detector",
          &Model::get_detector,
          &Model::set_detector)
      .add_property("crystal",
          &Model::get_crystal,
          &Model::set_crystal)
      .add_property("reflections",
          &Model::get_reflections,
          &Model::set_reflections)
      .add_property("image_data",
          &Model::get_image_data,
          &Model::set_image_data)
      .add_property("image_mask",
          &Model::get_image_mask,
          &Model::set_image_mask)
      .add_property("image_pred",
          &Model::get_image_pred)
      .add_property("mosaicity",
          &Model::get_mosaicity,
          &Model::set_mosaicity)
      .add_property("bandpass",
          &Model::get_bandpass,
          &Model::set_bandpass)
      .add_property("foreground_limit",
          &Model::get_foreground_limit,
          &Model::set_foreground_limit)
      .add_property("background_limit",
          &Model::get_background_limit,
          &Model::set_background_limit)
      .add_property("num_samples",
          &Model::get_num_samples,
          &Model::set_num_samples)
      .def("update",
          &Model::update, (
            arg("pixel_lookup") = true))
      .def("least_squares_score",
          &Model::least_squares_score)
      .def("maximum_likelihood_score",
          &Model::maximum_likelihood_score)
      .def("wavelength_values", &Model::wavelength_values)
      .def("wavelength_weights", &Model::wavelength_weights)
      .def("observed", &Model::observed)
      .def("predicted", &Model::predicted)
      .def("background", &Model::background)
      .def("intensity", &Model::intensity)
      .def("variance", &Model::variance)
      .def("num_foreground", &Model::num_foreground)
      .def("num_background", &Model::num_background)
      .def("shoebox", &Model::shoebox)
      .def("scale", &Model::scale)
      .def("success", &Model::success)
      .def("data", &Model::data)
      .def("mask", &Model::mask)
      .def("pred", &Model::pred)
      .def("__len__", &Model::size)
      ;
  }

}}} // namespace dials_scratch::examples::boost_python
