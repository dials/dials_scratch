#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <scitbx/constants.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/array_family/reflection_table.h>
#include <dxtbx/model/detector.h>
#include <dxtbx/model/scan.h>
#include <dxtbx/model/crystal.h>
#include <dxtbx/imageset.h>

#include <boost/asio.hpp>
#include <boost/thread.hpp>
#include <boost/lockfree/queue.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include <dials/algorithms/spot_prediction/pixel_to_miller_index.h>
#include <dials/algorithms/background/glm/robust_poisson_mean.h>
#include <dials/algorithms/integration/sum/summation.h>
#include <dials/model/data/mask_code.h>
#include <dials/model/data/shoebox.h>
#include <dials/algorithms/profile_model/gaussian_rs/mask_calculator.h>
#include <dials/algorithms/background/glm/robust_poisson_mean.h>
#include <dials/algorithms/profile_model/modeller/circle_sampler.h>
#include <dials/algorithms/profile_model/gaussian_rs/transform/transform.h>
#include <dials/algorithms/integration/sum/summation.h>
#include <dials/algorithms/integration/fit/fitting.h>
#include <dials/algorithms/centroid/centroid.h>
#include <dials/error.h>
#include <map>
#include <vector>
#include <boost/unordered_map.hpp>
#include <dials/algorithms/background/glm/creator.h>



using namespace boost::python;

namespace dials { namespace algorithms { namespace boost_python {


  using dxtbx::model::Beam;
  using dxtbx::model::Detector;
  using dxtbx::model::Goniometer;
  using dxtbx::model::Scan;

  using dxtbx::ImageSweep;
  using dxtbx::format::Image;

  using dials::model::Shoebox;
  using dials::algorithms::profile_model::gaussian_rs::MaskCalculator3D;
  using dials::algorithms::profile_model::gaussian_rs::transform::TransformSpec;
  using dials::algorithms::profile_model::gaussian_rs::CoordinateSystem;
  using dials::algorithms::profile_model::gaussian_rs::transform::TransformReverseNoModel;
  using dials::algorithms::profile_model::gaussian_rs::transform::TransformForward;
  using dials::algorithms::Creator;

  using dials::model::Centroid;
  using dials::model::Intensity;


  namespace detail {

    template <typename T>
    T median(const af::const_ref<T> &x) {
      af::shared<T> temp(x.begin(), x.end());
      std::nth_element(temp.begin(), temp.begin() + temp.size() / 2, temp.end());
      return temp[temp.size() / 2];
    }

  }

  struct vec3_int_hash : std::unary_function< vec3<int>, std::size_t> {
    std::size_t operator()(vec3<int> const& p) const {
        std::size_t seed = 0;
        boost::hash_combine(seed, p[0]);
        boost::hash_combine(seed, p[1]);
        boost::hash_combine(seed, p[2]);
        return seed;
    }
  };

  /* class PixelList { */
  /* public: */

  /*   PixelList( */
  /*       af::reflection_table reflections, */
  /*       const Beam &beam, */
  /*       const Detector &detector, */
  /*       const Goniometer &goniometer, */
  /*       const Scan &scan, */
  /*       double sigma_b, */
  /*       double sigma_m, */
  /*       int zstart) : */
  /*       beam_(beam), */
  /*       detector_(detector), */
  /*       goniometer_(goniometer), */
  /*       scan_(scan), */
  /*       sigma_b_(sigma_b), */
  /*       sigma_m_(sigma_m), */
  /*       spec_( */
  /*         beam, */
  /*         detector, */
  /*         goniometer, */
  /*         scan, */
  /*         sigma_b, */
  /*         sigma_m, */
  /*         4.0, */
  /*         10), */
  /*       zstart_(zstart) { */

  /*     MaskCalculator3D masker( */
  /*         beam, */
  /*         detector, */
  /*         goniometer, */
  /*         scan, */
  /*         sigma_b * 3, */
  /*         sigma_m * 3); */

  /*     af::const_ref< vec3<double> > s1    = reflections["s1"]; */
  /*     af::const_ref< vec3<double> > xyz   = reflections["xyzcal.px"]; */
  /*     af::const_ref< std::size_t >  panel = reflections["panel"]; */
  /*     af::const_ref< int6 >         bbox  = reflections["bbox"]; */

  /*     std::size_t zsize = scan.get_num_images(); */
  /*     std::size_t xsize = detector[0].get_image_size()[0]; */
  /*     std::size_t ysize = detector[0].get_image_size()[1]; */

  /*     lookup_pixls_.resize(reflections.size()); */
  /*     std::size_t num_pixels = 0; */

  /*     for (std::size_t r = 0; r < reflections.size(); ++r) { */
  /*       Shoebox<> shoebox(panel[r], bbox[r]); */
  /*       shoebox.allocate(); */
  /*       masker.single(shoebox, s1[r], xyz[r][2], panel[r]); */

  /*       af::const_ref< int, af::c_grid<3> > mask = shoebox.mask.const_ref(); */

  /*       int x0 = bbox[r][0]; */
  /*       int y0 = bbox[r][2]; */
  /*       int z0 = bbox[r][4]-zstart_; */

  /*       std::size_t num_foreground = 0; */
  /*       for (std::size_t i = 0; i < mask.size(); ++i) { */
  /*         if (mask[i] & 4) { */
  /*           num_foreground++; */
  /*         } */
  /*       } */

  /*       lookup_pixls_[r].resize(num_foreground); */

  /*       std::size_t i = 0; */
  /*       for (std::size_t z = 0; z < mask.accessor()[0]; ++z) { */
  /*         for (std::size_t y = 0; y < mask.accessor()[1]; ++y) { */
  /*           for (std::size_t x = 0; x < mask.accessor()[2]; ++x) { */
  /*             int zz = z0 + z; */
  /*             int yy = y0 + y; */
  /*             int xx = x0 + x; */
  /*             if (zz >= 0 && zz < zsize && */
  /*                 yy >= 0 && yy < ysize && */
  /*                 xx >= 0 && xx < xsize) { */
  /*               if (mask(z,y,x) & 4) { */
  /*                 std::size_t p = xx + yy * xsize + zz * ysize * xsize; */
  /*                 //lookup_count_[p] += 1; */
  /*                 //lookup_refls_[p].push_back(r); */
  /*                 DIALS_ASSERT(i < lookup_pixls_[r].size()); */
  /*                 lookup_pixls_[r][i] = p; */
  /*                 num_pixels++; */
  /*                 i++; */
  /*               } */
  /*             } */
  /*           } */
  /*         } */
  /*       } */
  /*     } */

  /*     lookup_count_.reserve(num_pixels); */

  /*     for (std::size_t i = 0; i < lookup_pixls_.size(); ++i) { */
  /*       for (std::size_t j = 0; j < lookup_pixls_[i].size(); ++j) { */
  /*         std::size_t k = lookup_pixls_[i][j]; */
  /*         lookup_count_[k] += 1; */
  /*       } */
  /*     } */


  /*     /1* typedef lookup_refls_type::iterator iterator; *1/ */

  /*     /1* for (iterator it = lookup_refls_.begin(); it != lookup_refls_.end(); ++it) { *1/ */
  /*     /1*   const std::vector<std::size_t> &reflections = it->second; *1/ */
  /*     /1*   for (std::size_t j = 0; j < reflections.size() - 1; ++j) { *1/ */
  /*     /1*     for (std::size_t i = 1; i < reflections.size(); ++i) { *1/ */
  /*     /1*       std::size_t idx_j = reflections[j]; *1/ */
  /*     /1*       std::size_t idx_i = reflections[i]; *1/ */
  /*     /1*       lookup_overlaps_[idx_j].push_back(idx_i); *1/ */
  /*     /1*       lookup_overlaps_[idx_i].push_back(idx_j); *1/ */
  /*     /1*     } *1/ */
  /*     /1*   } *1/ */
  /*     /1* } *1/ */
  /*   } */

  /*   void update_mask(af::ref< int, af::c_grid<3> > mask) { */
  /*     typedef lookup_count_type::iterator iterator; */
  /*     for (iterator it = lookup_count_.begin(); it != lookup_count_.end(); ++it) { */
  /*       std::size_t index = it->first; */
  /*       mask[index] |= 4; */
  /*     } */
  /*   } */

  /*   af::shared<bool> compute_background( */
  /*       const af::const_ref< double, af::c_grid<3> > &data, */
  /*       const af::const_ref< int,    af::c_grid<3> > &mask, */
  /*       af::reflection_table reflections) { */

  /*     std::size_t zsize = data.accessor()[0]; */
  /*     std::size_t ysize = data.accessor()[1]; */
  /*     std::size_t xsize = data.accessor()[2]; */

  /*     af::ref< double > background = reflections["background"]; */
  /*     af::ref< std::size_t > flags = reflections["flags"]; */
  /*     af::const_ref< int6 > bbox = reflections["bbox"]; */
  /*     af::shared<bool> success(reflections.size()); */

  /*     for (std::size_t i = 0; i < reflections.size(); ++i) { */
  /*       int x0 = std::max(bbox[i][0], 0); */
  /*       int x1 = std::min(bbox[i][1], (int)xsize); */
  /*       int y0 = std::max(bbox[i][2], 0); */
  /*       int y1 = std::min(bbox[i][3], (int)ysize); */
  /*       int z0 = std::max(bbox[i][4]-zstart_, 0); */
  /*       int z1 = std::min(bbox[i][5]-zstart_, (int)zsize); */

  /*       std::size_t num_background = 0; */
  /*       for (std::size_t z = z0; z < z1; ++z) { */
  /*         for (std::size_t y = y0; y < y1; ++y) { */
  /*           for (std::size_t x = x0; x < x1; ++x) { */
  /*             if (mask(z, y, x) == 1) { */
  /*               num_background++; */
  /*             } */
  /*           } */
  /*         } */
  /*       } */

  /*       af::shared<double> pixels(num_background); */
  /*       std::size_t count = 0; */
  /*       for (std::size_t z = z0; z < z1; ++z) { */
  /*         for (std::size_t y = y0; y < y1; ++y) { */
  /*           for (std::size_t x = x0; x < x1; ++x) { */
  /*             if (mask(z, y, x) == 1) { */
  /*               pixels[count++] = data(z, y, x); */
  /*             } */
  /*           } */
  /*         } */
  /*       } */

  /*       double median = detail::median(pixels.const_ref()); */
  /*       if (median == 0) { */
  /*         median = 1.0; */
  /*       } */

  /*       try { */

  /*         // Compute the result */
  /*         RobustPoissonMean result( */
  /*             pixels.const_ref(), */
  /*             median, */
  /*             1.345, */
  /*             1e-3, */
  /*             100); */
  /*         DIALS_ASSERT(result.converged()); */

  /*         // Compute the background */
  /*         double mean_background = result.mean(); */

  /*         background[i] = mean_background; */
  /*         success[i] = true; */
  /*       } catch(dials::error) { */
  /*         flags[i] |= af::DontIntegrate; */

  /*       } catch(scitbx::error) { */
  /*         flags[i] |= af::DontIntegrate; */
  /*       } */
  /*     } */
  /*     return success; */
  /*   } */

  /*   af::versa< double, af::c_grid<3> > get_profile_on_detector( */
  /*       boost::python::list reference, */
  /*       const CircleSampler &sampler, */
  /*       vec3<double> s1, */
  /*       vec3<double> xyz, */
  /*       double phi, */
  /*       std::size_t panel, */
  /*       int6 bbox) { */

  /*     vec3<double> m2 = goniometer_.get_rotation_axis(); */
  /*     vec3<double> s0 = beam_.get_s0(); */

  /*     std::size_t index = sampler.nearest(0, xyz); */

  /*     af::const_ref< double, af::c_grid<3> > profile = boost::python::extract< */
  /*       af::const_ref< double, af::c_grid<3> > >(reference[index])(); */

  /*     CoordinateSystem cs(m2, s0, s1, phi); */

  /*     TransformReverseNoModel transform(spec_, cs, bbox, panel, profile); */
  /*     return transform.profile(); */
  /*   } */


  /*   void get_data_on_detector( */
  /*       const af::const_ref< double, af::c_grid<3> > &data, */
  /*       const af::const_ref< int, af::c_grid<3> > &mask, */
  /*       int index, */
  /*       int6 bbox, */
  /*       af::ref< double, af::c_grid<3> > r_data, */
  /*       af::ref< double, af::c_grid<3> > r_bgrd, */
  /*       af::ref< int, af::c_grid<3> > r_mask, */
  /*       double background) { */

  /*     std::size_t zsize = data.accessor()[0]; */
  /*     std::size_t ysize = data.accessor()[1]; */
  /*     std::size_t xsize = data.accessor()[2]; */

  /*     int x0 = bbox[0]; */
  /*     int x1 = bbox[1]; */
  /*     int y0 = bbox[2]; */
  /*     int y1 = bbox[3]; */
  /*     int z0 = bbox[4]-zstart_; */
  /*     int z1 = bbox[5]-zstart_; */

  /*     DIALS_ASSERT(z1 > z0); */
  /*     DIALS_ASSERT(y1 > y0); */
  /*     DIALS_ASSERT(x1 > x0); */
  /*     vec3<std::size_t> size(z1 - z0, y1 - y0, x1 - x0); */
  /*     DIALS_ASSERT(r_data.accessor().all_eq(size)); */
  /*     DIALS_ASSERT(r_bgrd.accessor().all_eq(size)); */
  /*     DIALS_ASSERT(r_mask.accessor().all_eq(size)); */

  /*     /1* DIALS_ASSERT(x0 >= 0 && x1 <= data.accessor()[2]); *1/ */
  /*     /1* DIALS_ASSERT(y0 >= 0 && y1 <= data.accessor()[1]); *1/ */
  /*     /1* DIALS_ASSERT(z0 >= 0 && z1 <= data.accessor()[0]); *1/ */

  /*     std::size_t i = 0; */
  /*     for (int z = z0; z < z1; ++z) { */
  /*       for (int y = y0; y < y1; ++y) { */
  /*         for (int x = x0; x < x1; ++x) { */
  /*           if (z >= 0 && z < zsize && */
  /*               y >= 0 && y < ysize && */
  /*               x >= 0 && x < xsize) { */
  /*             r_data[i] = data(z,y,x); */
  /*             r_mask[i] = mask(z,y,x); */
  /*             r_bgrd[i] = background; */
  /*             if (r_mask[i] & 4) { */
  /*               r_mask[i] &= ~1; */
  /*             } */
  /*           } */
  /*           i++; */
  /*         } */
  /*       } */
  /*     } */

  /*     for (std::size_t i = 0; i < lookup_pixls_[index].size(); ++i) { */
  /*       std::size_t j = lookup_pixls_[index][i]; */
  /*       if (lookup_count_[j] == 1) { */
  /*         int z = (int)std::floor(j / (xsize * ysize)); */
  /*         int k = j - z * (xsize * ysize); */
  /*         int y = (int)std::floor(k / xsize); */
  /*         int x = k - y * xsize; */
  /*         int zz = z - z0; */
  /*         int yy = y - y0; */
  /*         int xx = x - x0; */
  /*         DXTBX_ASSERT(zz >= 0 && zz < r_mask.accessor()[0]); */
  /*         DXTBX_ASSERT(yy >= 0 && yy < r_mask.accessor()[1]); */
  /*         DXTBX_ASSERT(xx >= 0 && xx < r_mask.accessor()[2]); */
  /*         r_mask(zz, yy, xx) |= 1; */
  /*       } */
  /*     } */
  /*   } */


  /*   af::shared<bool> compute_intensity( */
  /*       const af::const_ref< double, af::c_grid<3> > &data, */
  /*       const af::const_ref< int, af::c_grid<3> > &mask, */
  /*       af::reflection_table reflections, */
  /*       boost::python::list reference, */
  /*       const CircleSampler &sampler */
  /*       ) { */

  /*     af::const_ref< vec3<double> > s1 = reflections["s1"]; */
  /*     af::const_ref< vec3<double> > xyzpx = reflections["xyzcal.px"]; */
  /*     af::const_ref< vec3<double> > xyzmm = reflections["xyzcal.mm"]; */
  /*     af::const_ref< std::size_t > panel = reflections["panel"]; */
  /*     af::const_ref< int6 > bbox = reflections["bbox"]; */
  /*     af::const_ref< double > background = reflections["background"]; */

  /*     af::ref<double> Isum = reflections["intensity.sum.value"]; */
  /*     af::ref<double> Vsum = reflections["intensity.sum.variance"]; */
  /*     af::ref<double> Iprf = reflections["intensity.prf.value"]; */
  /*     af::ref<double> Vprf = reflections["intensity.prf.variance"]; */
  /*     af::ref<double> part = reflections["partiality"]; */
  /*     af::ref<std::size_t> flags = reflections["flags"]; */
  /*     af::shared<bool> success(reflections.size()); */

  /*     for (std::size_t i = 0; i < reflections.size(); ++i) { */
  /*       if (flags[i] & af::DontIntegrate) { */
  /*         continue; */
  /*       } */

  /*       af::versa<double, af::c_grid<3> > profile = get_profile_on_detector( */
  /*         reference, */
  /*         sampler, */
  /*         s1[i], */
  /*         xyzpx[i], */
  /*         xyzmm[i][2], */
  /*         panel[i], */
  /*         bbox[i]); */

  /*       int x0 = bbox[i][0]; */
  /*       int x1 = bbox[i][1]; */
  /*       int y0 = bbox[i][2]; */
  /*       int y1 = bbox[i][3]; */
  /*       int z0 = bbox[i][4]-zstart_; */
  /*       int z1 = bbox[i][5]-zstart_; */
  /*       DIALS_ASSERT(x1 > x0); */
  /*       DIALS_ASSERT(y1 > y0); */
  /*       DIALS_ASSERT(z1 > z0); */
  /*       std::size_t xsize = x1 - x0; */
  /*       std::size_t ysize = y1 - y0; */
  /*       std::size_t zsize = z1 - z0; */

  /*       af::c_grid<3> grid(zsize, ysize, xsize); */
  /*       af::versa< double, af::c_grid<3> > r_data(grid); */
  /*       af::versa< double, af::c_grid<3> > r_bgrd(grid); */
  /*       af::versa< int,    af::c_grid<3> > r_mask(grid); */

  /*       try { */
  /*         get_data_on_detector( */
  /*           data, */
  /*           mask, */
  /*           i, */
  /*           bbox[i], */
  /*           r_data.ref(), */
  /*           r_bgrd.ref(), */
  /*           r_mask.ref(), */
  /*           background[i]); */
  /*       } catch (dials::error) { */
  /*         continue; */
  /*       } */

  /*       double p = 0; */
  /*       for (std::size_t j = 0; j < profile.size(); ++j) { */
  /*         if (r_mask[j] == 5) { */
  /*           p += profile[j]; */
  /*         } */
  /*       } */
  /*       part[i] = p; */

  /*       Summation<double> summation(r_data.const_ref(), r_bgrd.const_ref(), r_mask.const_ref()); */
  /*       Isum[i] = summation.intensity(); */
  /*       Vsum[i] = summation.variance(); */
  /*       flags[i] |= af::IntegratedSum; */

  /*       /1* try { *1/ */

  /*       /1*   af::versa<bool, af::c_grid<3> > r_mask2(r_mask.accessor()); *1/ */
  /*       /1*   for (std::size_t j = 0; j < r_mask.size(); ++j) { *1/ */
  /*       /1*     r_mask2[j] = r_mask[j] == 5; *1/ */
  /*       /1*   } *1/ */

  /*       /1*   ProfileFitting<double> fitting(profile.const_ref(), r_mask2.const_ref(), r_data.const_ref(), r_bgrd.const_ref(), 1e-3, 10); *1/ */
  /*       /1*   DIALS_ASSERT(fitting.niter() < 10); *1/ */
  /*       /1*   Iprf[i] = fitting.intensity(); *1/ */
  /*       /1*   Vprf[i] = fitting.variance(); *1/ */
  /*       /1*   success[i] = true; *1/ */
  /*       /1*   flags[i] |= af::IntegratedSum; *1/ */
  /*       /1*   flags[i] |= af::IntegratedPrf; *1/ */
  /*       /1* } catch (dials::error) { *1/ */

  /*       /1* } *1/ */
  /*     } */


  /*     return success; */
  /*   } */

  /* protected: */

  /*   typedef boost::unordered_map< std::size_t, std::vector< std::size_t > > lookup_refls_type; */
  /*   typedef boost::unordered_map< std::size_t, std::size_t >                lookup_count_type; */
  /*   typedef std::vector< std::vector< std::size_t > > lookup_pixls_type; */
  /*   typedef boost::unordered_map< std::size_t, std::vector< std::size_t > > lookup_overlaps_type; */


  /*   Beam beam_; */
  /*   Detector detector_; */
  /*   Goniometer goniometer_; */
  /*   Scan scan_; */
  /*   double sigma_b_; */
  /*   double sigma_m_; */
  /*   TransformSpec spec_; */
  /*   int zstart_; */

  /*   lookup_count_type lookup_count_; */
  /*   //lookup_refls_type lookup_refls_; */
  /*   lookup_pixls_type lookup_pixls_; */
  /*   //lookup_overlaps_type lookup_overlaps_; */

  /* }; */



  class ThreadPool {
  public:

    ThreadPool(std::size_t N)
        : work_(io_service_),
          started_(0),
          finished_(0) {
      for (std::size_t i = 0; i < N; ++i) {
        threads_.create_thread(
            boost::bind(
              &boost::asio::io_service::run,
              &io_service_));
      }
    }

    ~ThreadPool() {
      io_service_.stop();
      try {
        threads_.join_all();
      } catch (const std::exception&) {
        // pass
      }
    }

    template <typename Function>
    void post(Function function) {
      started_++;
      io_service_.post(FunctionRunner<Function>(function, finished_));
    }

    void wait() {
      while (finished_ < started_);
    }

  protected:

    template <typename Function>
    class FunctionRunner {
    public:

      FunctionRunner(Function function, boost::atomic<std::size_t> &counter)
        : function_(function),
          counter_(counter) {}

      void operator()() {
        function_();
        counter_++;
      }

    protected:

      Function function_;
      boost::atomic<std::size_t> &counter_;
    };

    boost::asio::io_service io_service_;
    boost::asio::io_service::work work_;
    boost::thread_group threads_;
    std::size_t started_;
    boost::atomic<std::size_t> finished_;

  };


  class Reflection {
  public:

    typedef typename af::reflection_table_type_generator::data_type data_type;
    typedef std::map<std::string, data_type> map_type;

    typedef typename map_type::key_type key_type;
    typedef typename map_type::mapped_type mapped_type;
    typedef typename map_type::value_type map_value_type;
    typedef typename map_type::iterator iterator;
    typedef typename map_type::const_iterator const_iterator;
    typedef typename map_type::size_type size_type;

    Reflection()
      : data_(boost::make_shared<map_type>()){}

    /**
     * Access a value by key
     * @param key The column name
     * @returns The proxy object to access the value
     */
    const mapped_type& operator[](const key_type &key) const {
      return data_->operator[](key);
    }

    /**
     * Access a value by key
     * @param key The column name
     * @returns The proxy object to access the value
     */
    mapped_type& operator[](const key_type &key) {
      return data_->operator[](key);
    }

    /**
     * Access a value by key
     * @param key The column name
     * @returns The value.
     */
    template <typename T>
    T& get(const key_type &key) {
      iterator it = find(key);
      DIALS_ASSERT(it != end());
      return boost::get<T>(it->second);
    }

    /**
     * Access a value by key
     * @param key The column name
     * @returns The value.
     */
    template <typename T>
    const T& get(const key_type &key) const {
      const_iterator it = find(key);
      DIALS_ASSERT(it != end());
      return boost::get<T>(it->second);
    }

    /** @returns An iterator to the beginning of the column map */
    iterator begin() {
      return data_->begin();
    }

    /** @returns An iterator to the end of the column map */
    iterator end() {
      return data_->end();
    }

    /** @returns A const iterator to the beginning of the column map */
    const_iterator begin() const {
      return data_->begin();
    }

    /** @returns A const iterator to the end of the column map */
    const_iterator end() const {
      return data_->end();
    }

    /** @returns The number of values in the table */
    size_type size() const {
      return data_->size();
    }

    /** @returns Is the table empty */
    bool empty() const {
      return data_->empty();
    }

    /** @returns The number of columns matching the key (0 or 1) */
    size_type count(const key_type &key) const {
      return data_->count(key);
    }

    /**
     * Find a column matching the key
     * @param key The column name
     * @returns An iterator to the column
     */
    iterator find(const key_type &key) {
      return data_->find(key);
    }

    /**
     * Find a column matching the key
     * @param key The column name
     * @returns A const iterator to the column
     */
    const_iterator find(const key_type &key) const {
      return data_->find(key);
    }

    /**
     * Erase a column from the table.
     * @param key The column name
     * @returns The number of columns removed
     */
    size_type erase(const key_type &key) {
      return data_->erase(key);
    }

    /** Clear the table */
    void clear() {
      data_->clear();
    }

    /** @returns Does the table contain the key. */
    bool contains(const key_type &key) const {
      const_iterator it = find(key);
      return it != end();
    }

  protected:

    boost::shared_ptr<map_type> data_;

  };

  struct element_to_variant_visitor : public boost::static_visitor<Reflection::data_type> {
    std::size_t n_;
    element_to_variant_visitor(std::size_t n) : n_(n) {}
    template <typename T>
    Reflection::data_type operator () (T &col) {
      DIALS_ASSERT(n_ < col.size());
      return Reflection::data_type(col[n_]);
    }
  };

  struct variant_to_element_visitor : public boost::static_visitor<void> {
    af::reflection_table table_;
    std::size_t n_;
    Reflection::key_type key_;
    variant_to_element_visitor(af::reflection_table table, std::size_t n, Reflection::key_type key) :
      table_(table),
      n_(n),
      key_(key) {}
    template <typename T>
    void operator () (const T &item) {
      af::ref<T> col = table_[key_];
      DIALS_ASSERT(n_ < col.size());
      col[n_] = item;
    }
  };

  struct to_table_visitor : public boost::static_visitor<void> {
    af::reflection_table table_;
    std::string key_;
    to_table_visitor(af::reflection_table table, std::string key) :
      table_(table),
      key_(key) {
      DIALS_ASSERT(table_.size() == 1);
    }
    template <typename T>
    void operator () (const T &item) {
      table_[key_] = af::shared<T>(1, item);
    }
  };


  Reflection get_reflection(af::reflection_table table, std::size_t i) {

    typedef typename af::reflection_table::const_iterator iterator;
    DIALS_ASSERT(i < table.size());
    Reflection result;
    element_to_variant_visitor visitor(i);
    for (iterator it = table.begin(); it != table.end(); ++it) {
      result[it->first] = it->second.apply_visitor(visitor);
    }
    return result;
  }


  void set_reflection(af::reflection_table table, std::size_t i, Reflection item) {
    typedef typename Reflection::const_iterator iterator;
    DIALS_ASSERT(i < table.size());
    for (iterator it = item.begin(); it != item.end(); ++it) {
      variant_to_element_visitor visitor(table, i, it->first);
      it->second.apply_visitor(visitor);
    }
  }

  af::reflection_table Reflection_to_table(const Reflection &self) {
    typedef typename Reflection::const_iterator iterator;
    af::reflection_table result(1);
    for (iterator it = self.begin(); it != self.end(); ++it) {
      to_table_visitor visitor(result, it->first);
      it->second.apply_visitor(visitor);
    }
    return result;
  }


  struct item_to_object_visitor : public boost::static_visitor<object> {
    template <typename T>
    object operator () (T &col) {
      return object(col);
    }
  };

  boost::python::object Reflection_get(const Reflection &self, std::string name) {
    typename Reflection::mapped_type item = self[name];
    item_to_object_visitor visitor;
    return item.apply_visitor(visitor);
  }

  void Reflection_set_bool(Reflection &self, std::string name, bool item) {
    self[name] = Reflection::data_type(item);
  }

  void Reflection_set_int(Reflection &self, std::string name, int item) {
    self[name] = Reflection::data_type(item);
  }

  void Reflection_set_size_t(Reflection &self, std::string name, std::size_t item) {
    self[name] = Reflection::data_type(item);
  }

  void Reflection_set_double(Reflection &self, std::string name, double item) {
    self[name] = Reflection::data_type(item);
  }

  void Reflection_set_string(Reflection &self, std::string name, std::string item) {
    self[name] = Reflection::data_type(item);
  }

  void Reflection_set_vec2_double(Reflection &self, std::string name, vec2<double> item) {
    self[name] = Reflection::data_type(item);
  }

  void Reflection_set_vec3_double(Reflection &self, std::string name, vec3<double> item) {
    self[name] = Reflection::data_type(item);
  }

  void Reflection_set_mat3_double(Reflection &self, std::string name, mat3<double> item) {
    self[name] = Reflection::data_type(item);
  }

  void Reflection_set_int6(Reflection &self, std::string name, int6 item) {
    self[name] = Reflection::data_type(item);
  }

  void Reflection_set_miller_index(Reflection &self, std::string name, cctbx::miller::index<> item) {
    self[name] = Reflection::data_type(item);
  }

  void Reflection_set_shoebox(Reflection &self, std::string name, Shoebox<> item) {
    self[name] = Reflection::data_type(item);
  }


  class MaskCalculatorIface {
  public:

    virtual void operator()(Reflection &reflection) const = 0;

  };

  class BackgroundCalculatorIface {
  public:

    virtual void operator()(Reflection &reflection) const = 0;

  };

  class IntensityCalculatorIface {
  public:

    virtual void operator()(Reflection &reflection) const = 0;

  };


  class MaskCalculator : public MaskCalculatorIface {
  public:

    MaskCalculator(const BeamBase &beam,
                   const Detector &detector,
                   const Goniometer &gonio,
                   const Scan &scan,
                   double delta_b,
                   double delta_m)
      : func_(beam,
              detector,
              gonio,
              scan,
              delta_b,
              delta_m) {}

    virtual void operator()(Reflection &reflection) const {
      func_.single(
        reflection.get< Shoebox<> >("shoebox"),
        reflection.get< vec3<double> >("s1"),
        reflection.get< vec3<double> >("xyzcal.px")[2],
        reflection.get< std::size_t >("panel"));
    }

  protected:

    MaskCalculator3D func_;

  };

  class BackgroundCalculator : public BackgroundCalculatorIface {
  public:

    BackgroundCalculator(
          Creator::Model model,
          double tuning_constant,
          std::size_t max_iter,
          std::size_t min_pixels)
      : creator_(
          model,
          tuning_constant,
          max_iter,
          min_pixels) {}

    virtual void operator()(Reflection &reflection) const {
      creator_.single(reflection.get< Shoebox<> >("shoebox"));
    }

  protected:

    Creator creator_;
  };


  class Reference {
  public:

    typedef af::versa< double, af::c_grid<3> > data_type;
    typedef af::versa< bool, af::c_grid<3> > mask_type;

    void append(af::const_ref< double, af::c_grid<3> > data,
                af::const_ref< bool, af::c_grid<3> > mask) {
      data_type d(data.accessor());
      mask_type m(mask.accessor());
      std::copy(data.begin(), data.end(), d.begin());
      std::copy(mask.begin(), mask.end(), m.begin());
      data_.push_back(d);
      mask_.push_back(m);
    }

    af::const_ref< double, af::c_grid<3> > data(std::size_t index) const {
      DIALS_ASSERT(index < data_.size());
      return data_[index].const_ref();
    }

    af::const_ref< bool, af::c_grid<3> > mask(std::size_t index) const {
      DIALS_ASSERT(index < mask_.size());
      return mask_[index].const_ref();
    }

  protected:

    af::shared<data_type> data_;
    af::shared<mask_type> mask_;
  };

  namespace detail {

    struct check_mask_code {
      int mask_code;
      check_mask_code(int code) : mask_code(code) {}
      bool operator()(int a) const {
        return ((a & mask_code) == mask_code);
      }
    };

    struct check_either_mask_code {
      int mask_code1;
      int mask_code2;
      check_either_mask_code(int code1, int code2) : mask_code1(code1), mask_code2(code2) {}
      bool operator()(int a) const {
        return ((a & mask_code1) == mask_code1) || ((a & mask_code2) == mask_code2);
      }
    };

  }

  class IntensityCalculator : public IntensityCalculatorIface {
  public:

    IntensityCalculator(
        const Reference &reference,
        const CircleSampler &sampler,
        const TransformSpec &spec)
      : reference_(reference),
        sampler_(sampler),
        spec_(spec) {
    }

    virtual void operator()(Reflection &reflection) const {

      typedef af::const_ref< double, af::c_grid<3> > data_const_reference;
      typedef af::const_ref< bool, af::c_grid<3> > mask_const_reference;

      reflection["flags"] = reflection.get<std::size_t>("flags") & ~af::IntegratedPrf;

      vec3<double> s1 = reflection.get< vec3<double> >("s1");
      double phi = reflection.get< vec3<double> >("xyzcal.mm")[2];
      vec3<double> xyz = reflection.get< vec3<double> >("xyzcal.px");
      int6 bbox = reflection.get<int6>("bbox");
      std::size_t panel = reflection.get<std::size_t>("panel");
      Shoebox<> sbox = reflection.get< Shoebox<> >("shoebox");
      DIALS_ASSERT(sbox.is_consistent());

      if (check(sbox) == false) {
        return;
      }

      // Create the coordinate system
      vec3<double> m2 = spec_.goniometer().get_rotation_axis();
      vec3<double> s0 = spec_.beam()->get_s0();
      CoordinateSystem cs(m2, s0, s1, phi);

      // Get the reference profiles
      std::size_t index = sampler_.nearest(sbox.panel, xyz);
      data_const_reference p = reference_.data(index);
      mask_const_reference mask1 = reference_.mask(index);

      // Create the data array
      af::versa< double, af::c_grid<3> > data(sbox.data.accessor());
      std::copy(
          sbox.data.begin(),
          sbox.data.end(),
          data.begin());

      // Create the background array
      af::versa< double, af::c_grid<3> > background(sbox.background.accessor());
      std::copy(
          sbox.background.begin(),
          sbox.background.end(),
          background.begin());

      // Create the mask array
      af::versa< bool, af::c_grid<3> > mask(sbox.mask.accessor());
      std::transform(
          sbox.mask.begin(),
          sbox.mask.end(),
          mask.begin(),
          detail::check_mask_code(Valid | Foreground));

      // Compute the transform
      TransformForward<double> transform(
          spec_,
          cs,
          sbox.bbox,
          sbox.panel,
          data.const_ref(),
          background.const_ref(),
          mask.const_ref());

      // Get the transformed shoebox
      data_const_reference c = transform.profile().const_ref();
      data_const_reference b = transform.background().const_ref();
      mask_const_reference mask2 = transform.mask().const_ref();
      af::versa< bool, af::c_grid<3> > m(mask2.accessor());
      DIALS_ASSERT(mask1.size() == mask2.size());
      for (std::size_t j = 0; j < m.size(); ++j) {
        m[j] = mask1[j] && mask2[j];
      }

      ProfileFitter<double> fit(
          c,
          b,
          m.const_ref(),
          p,
          1e-3,
          100);
      DIALS_ASSERT(fit.niter() < 100);

      // Set the integrated flag
      reflection["intensity.prf.value"] = fit.intensity()[0];
      reflection["intensity.prf.variance"] = fit.variance()[0];
      reflection["intensity.prf.correlation"] = fit.correlation();
      reflection["flags"] = reflection.get<std::size_t>("flags") | af::IntegratedPrf;
    }

  protected:

    bool check(const Shoebox<> &sbox) const {

      // Check if the bounding box is in the image
      bool bbox_valid =
        sbox.bbox[0] >= 0 &&
        sbox.bbox[2] >= 0 &&
        sbox.bbox[1] <= spec_.detector()[sbox.panel].get_image_size()[0] &&
        sbox.bbox[3] <= spec_.detector()[sbox.panel].get_image_size()[1];

      // Check if all pixels are valid
      bool pixels_valid = true;
      for (std::size_t i = 0; i < sbox.mask.size(); ++i) {
        if (sbox.mask[i] & Foreground && !(sbox.mask[i] & Valid)) {
          pixels_valid = false;
          break;
        }
      }

      // Return whether to use or not
      return bbox_valid && pixels_valid;
    }

    Reference reference_;
    CircleSampler sampler_;
    TransformSpec spec_;

  };

  class ReflectionIntegrator {
  public:

    ReflectionIntegrator(
          const MaskCalculatorIface &compute_mask,
          const BackgroundCalculatorIface &compute_background,
          const IntensityCalculatorIface &compute_intensity,
          const af::const_ref< double, af::c_grid<3> > &buffer,
          int zstart,
          double underload,
          double overload)
      : compute_mask_(compute_mask),
        compute_background_(compute_background),
        compute_intensity_(compute_intensity),
        buffer_(buffer),
        zstart_(zstart),
        underload_(underload),
        overload_(overload) {}

    void operator()(Reflection &reflection) const {

      extract_shoebox(
          buffer_,
          reflection,
          zstart_,
          underload_,
          overload_);

      compute_mask_(reflection);

      try {
        compute_background_(reflection);
      } catch (dials::error) {
        return;
      }

      compute_centroid(reflection);

      compute_summed_intensity(reflection);

      try {
        compute_intensity_(reflection);
      } catch (dials::error) {
        return;
      }

      reflection.erase("shoebox");
    }


  protected:

    void extract_shoebox(
          const af::const_ref< double, af::c_grid<3> > &buffer,
          Reflection &reflection,
          int zstart,
          double underload,
          double overload) const {
      std::size_t panel = reflection.get<std::size_t>("panel");
      int6 bbox = reflection.get<int6>("bbox");
      Shoebox<> shoebox(panel, bbox);
      shoebox.allocate();
      af::ref< float, af::c_grid<3> > data = shoebox.data.ref();
      af::ref< int,   af::c_grid<3> > mask = shoebox.mask.ref();
      int x0 = bbox[0];
      int x1 = bbox[1];
      int y0 = bbox[2];
      int y1 = bbox[3];
      int z0 = bbox[4];
      int z1 = bbox[5];
      DIALS_ASSERT(x1 > x0);
      DIALS_ASSERT(y1 > y0);
      DIALS_ASSERT(z1 > z0);
      std::size_t zsize = z1 - z0;
      std::size_t ysize = y1 - y0;
      std::size_t xsize = x1 - x0;
      DIALS_ASSERT(zsize == data.accessor()[0]);
      DIALS_ASSERT(ysize == data.accessor()[1]);
      DIALS_ASSERT(xsize == data.accessor()[2]);
      for (std::size_t k = 0; k < zsize; ++k) {
        for (std::size_t j = 0; j < ysize; ++j) {
          for (std::size_t i = 0; i < xsize; ++i) {
            int kk = z0 + k - zstart;
            int jj = y0 + j;
            int ii = x0 + i;
            if (kk >= 0 &&
                jj >= 0 &&
                ii >= 0 &&
                kk < buffer.accessor()[0] &&
                jj < buffer.accessor()[1] &&
                ii < buffer.accessor()[2]) {
              double d = buffer(kk, jj, ii);
              int m = (d > underload && d < overload) ? Valid : 0;
              data(k,j,i) = d;
              mask(k,j,i) = m;
            } else {
              data(k,j,i) = 0;
              mask(k,j,i) = 0;
            }
          }
        }
      }
      reflection["shoebox"] = shoebox;
    }

    void compute_centroid(Reflection &r) const {
      Shoebox<> shoebox = r.get< Shoebox<> >("shoebox");
      Centroid centroid = shoebox.centroid_foreground_minus_background();
      r["xyzobs.px.value"] = centroid.px.position;
      r["xyzobs.px.variance"] = centroid.px.variance;
    }

    void compute_summed_intensity(Reflection &r) const {
      Shoebox<> shoebox = r.get< Shoebox<> >("shoebox");
      Intensity intensity = shoebox.summed_intensity();
      r["intensity.sum.value"] = intensity.observed.value;
      r["intensity.sum.variance"] = intensity.observed.variance;
    }

    const MaskCalculatorIface &compute_mask_;
    const BackgroundCalculatorIface &compute_background_;
    const IntensityCalculatorIface &compute_intensity_;
    af::const_ref< double, af::c_grid<3> > buffer_;
    int zstart_;
    double underload_;
    double overload_;

  };

  namespace detail {

    class sort_by_frame {
    public:
      sort_by_frame(af::const_ref<int6> bbox)
        : bbox_(bbox) {}

      bool operator()(std::size_t a, std::size_t b) const {
        return bbox_[a][5] < bbox_[b][5];
      }

    protected:

      af::const_ref<int6> bbox_;
    };
  }

  class Lookup {
  public:

    Lookup(af::const_ref<int6> bbox, int zstart, std::size_t n)
      : indices_(bbox.size()) {
      for (std::size_t i = 0; i < indices_.size(); ++i) {
        indices_[i] = i;
      }
      std::sort(indices_.begin(), indices_.end(), detail::sort_by_frame(bbox));
      DIALS_ASSERT(bbox[indices_.front()][5] - zstart >= 1);
      DIALS_ASSERT(bbox[indices_.back()][5]  - zstart <= n);
      std::size_t i = 0;
      offset_.push_back(0);
      for (std::size_t j = 0; j < n; ++j) {
        while (i < indices_.size() && bbox[indices_[i]][5] - zstart <= j+1) i++;
        offset_.push_back(i);
      }
      DIALS_ASSERT(offset_.size() == n + 1);
      DIALS_ASSERT(offset_.back() == indices_.size());
    }


    af::const_ref<std::size_t> indices(std::size_t z) const {
      DIALS_ASSERT(z < offset_.size()-1);
      DIALS_ASSERT(offset_[z+1] >= offset_[z]);
      std::size_t i = offset_[z];
      std::size_t n = offset_[z+1] - i;
      return af::const_ref<std::size_t>(&indices_[i], n);
    }

  private:

    std::vector<std::size_t> indices_;
    std::vector<std::size_t> offset_;
  };

  class Integrator {
  public:

    Integrator(
          af::reflection_table reflections,
          ImageSweep imageset,
          const MaskCalculatorIface &compute_mask,
          const BackgroundCalculatorIface &compute_background,
          const IntensityCalculatorIface &compute_intensity)
      : reflections_(reflections) {

      ThreadPool pool(8);

      DIALS_ASSERT(imageset.get_detector() != NULL);
      DIALS_ASSERT(imageset.get_scan() != NULL);

      std::size_t xsize = (*imageset.get_detector())[0].get_image_size()[0];
      std::size_t ysize = (*imageset.get_detector())[0].get_image_size()[1];
      std::size_t zsize = imageset.size();
      int zstart = imageset.get_scan()->get_array_range()[0];

      DIALS_ASSERT(xsize > 0);
      DIALS_ASSERT(ysize > 0);
      DIALS_ASSERT(zsize > 0);

      af::versa< double, af::c_grid<3> > buffer(af::c_grid<3>(zsize, ysize, xsize));

      double underload = (*imageset.get_detector())[0].get_trusted_range()[0];
      double overload  = (*imageset.get_detector())[0].get_trusted_range()[1];

      af::shared<int6> bbox = reflections.get<int6>("bbox");

      Lookup lookup(bbox.const_ref(), zstart, buffer.accessor()[0]);

      ReflectionIntegrator integrator(
          compute_mask,
          compute_background,
          compute_intensity,
          buffer.const_ref(),
          zstart,
          underload,
          overload);

      std::vector< Reflection > reflection_buffer;
      for (std::size_t i = 0; i < reflections.size(); ++i) {
        reflection_buffer.push_back(get_reflection(reflections, i));
      }

      af::ref<std::size_t> flags = reflections["flags"];

      for (std::size_t i = 0; i < zsize; ++i) {
        std::cout << zstart + i << std::endl;
        Image<double> data = imageset.get_corrected_data(i);
        DIALS_ASSERT(data.tile(0).data().size() == ysize * xsize);
        for (std::size_t j = 0; j < ysize * xsize; ++j) {
          buffer[i * (xsize*ysize) + j] = data.tile(0).data()[j];
        }
        af::const_ref<std::size_t> indices = lookup.indices(i);
        for (std::size_t j = 0; j < indices.size(); ++j) {

          std::size_t k = indices[j];

          if (flags[k] & af::DontIntegrate) {
            flags[k] &= ~af::IntegratedSum;
            flags[k] &= ~af::IntegratedPrf;
            continue;
          }

          pool.post(
                boost::bind(
                  &ReflectionIntegrator::operator(),
                  integrator,
                  boost::ref(reflection_buffer[k])));
        }
      }

      pool.wait();

      for (std::size_t i = 0; i < reflections.size(); ++i) {
        set_reflection(reflections, i, reflection_buffer[i]);
      }
    }

    af::reflection_table reflections() const {
      return reflections_;
    }

  protected:

    void read_image_data(
          ImageSweep imageset,
          af::ref< double, af::c_grid<3> > buffer) const {
      int zstart = imageset.get_scan()->get_array_range()[0];
      std::size_t zsize = buffer.accessor()[0];
      std::size_t ysize = buffer.accessor()[1];
      std::size_t xsize = buffer.accessor()[2];
      for (std::size_t i = 0; i < zsize; ++i) {
        std::cout << zstart + i << std::endl;
        Image<double> data = imageset.get_corrected_data(i);
        DIALS_ASSERT(data.tile(0).data().size() == ysize * xsize);
        for (std::size_t j = 0; j < ysize * xsize; ++j) {
          buffer[i * (xsize*ysize) + j] = data.tile(0).data()[j];
        }
      }
    }


    af::reflection_table reflections_;
  };



  BOOST_PYTHON_MODULE(dials_scratch_jmp_fitting_ext)
  {

    /* class_<PixelList>("PixelList", no_init) */
    /*   .def(init< */
    /*     af::reflection_table, */
    /*     const Beam &, */
    /*     const Detector &, */
    /*     const Goniometer &, */
    /*     const Scan &, */
    /*     double, */
    /*     double, */
    /*     int>()) */
    /*   .def("update_mask", &PixelList::update_mask) */
    /*   .def("compute_background", &PixelList::compute_background) */
    /*   .def("compute_intensity", &PixelList::compute_intensity) */
    /*   ; */

    class_<Reflection>("Reflection")
      .def("get",
          &Reflection_get)
      .def("set_bool",
          &Reflection_set_bool)
      .def("set_int",
          &Reflection_set_int)
      .def("set_size_t",
          &Reflection_set_size_t)
      .def("set_double",
          &Reflection_set_double)
      .def("set_string",
          &Reflection_set_string)
      .def("set_vec2_double",
          &Reflection_set_vec2_double)
      .def("set_vec3_double",
          &Reflection_set_vec3_double)
      .def("set_mat3_double",
          &Reflection_set_mat3_double)
      .def("set_int6",
          &Reflection_set_int6)
      .def("set_miller_index",
          &Reflection_set_miller_index)
      .def("set_shoebox",
          &Reflection_set_shoebox)
      .def("to_table",
          &Reflection_to_table)
      ;

    def("set_reflection", &set_reflection);
    def("get_reflection", &get_reflection);

    class_<MaskCalculatorIface, boost::noncopyable>("MaskCalculatorIface", no_init)
      ;

    class_<BackgroundCalculatorIface, boost::noncopyable>("BackgroundCalculatorIface", no_init)
      ;

    class_<IntensityCalculatorIface, boost::noncopyable>("IntensityCalculatorIface", no_init)
      ;

    class_<MaskCalculator, bases<MaskCalculatorIface> >("MaskCalculator", no_init)
      .def(init<const BeamBase&,
                const Detector&,
                const Goniometer&,
                const Scan&,
                double,
                double>())
      ;

    class_<BackgroundCalculator, bases<BackgroundCalculatorIface> >("BackgroundCalculator", no_init)
      .def(init<Creator::Model,
                double,
                std::size_t,
                std::size_t>())
      ;


    class_<Reference>("Reference")
      .def("append", &Reference::append)
      ;

    class_<IntensityCalculator, bases<IntensityCalculatorIface> >("IntensityCalculator", no_init)
      .def(init<
          const Reference&,
          const CircleSampler&,
          const TransformSpec&>())
      ;


    class_<Integrator>("Integrator", no_init)
      .def(init<
          const af::reflection_table&,
          ImageSweep,
          const MaskCalculatorIface&,
          const BackgroundCalculatorIface&,
          const IntensityCalculatorIface&>((
              arg("reflections"),
              arg("imageset"),
              arg("compute_mask"),
              arg("compute_background"),
              arg("compute_intensity"))))
      .def("reflections",
          &Integrator::reflections)
      ;

  }

}}} // namespace dials_scratch::examples::boost_python
