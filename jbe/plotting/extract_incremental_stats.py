import gc
import resource
import json
from dxtbx.serialize import load
from dials.array_family import flex
from dials.algorithms.scaling.scaling_library import scaled_data_as_miller_array
from dials.report.plots import ResolutionPlotsAndStats
import iotbx


def merging_stats_from_scaled_array(scaled_miller_array):
    result = iotbx.merging_statistics.dataset_statistics(
        i_obs=scaled_miller_array,
        n_bins=20,
        anomalous=False,
        sigma_filtering=None,
        eliminate_sys_absent=False,
        use_internal_variance=False,
        cc_one_half_significance_level=0.01,
    )

    intensities_anom = scaled_miller_array.as_anomalous_array()
    intensities_anom = intensities_anom.map_to_asu().customized_copy(
        info=scaled_miller_array.info()
    )
    anom_result = iotbx.merging_statistics.dataset_statistics(
        i_obs=intensities_anom,
        n_bins=20,
        anomalous=True,
        sigma_filtering=None,
        cc_one_half_significance_level=0.01,
        eliminate_sys_absent=False,
        use_internal_variance=False,
    )
    return result, anom_result


filenames = (
    "/Users/whi10850/Documents/test_data/multi_example/integrated_files/inc_test/1/",
    "/Users/whi10850/Documents/test_data/multi_example/integrated_files/inc_test/2/",
)

cc_one_half_data = []
i_over_sig_data = []

for filepath in filenames:
    experiments = load.experiment_list(
        filepath + "scaled_experiments.json", check_format=False
    )
    reflections = flex.reflection_table.from_pickle(filepath + "scaled.pickle")
    print("memory: %s" % int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
    arr = scaled_data_as_miller_array([reflections], experiments)
    norm_stats, anom_stats = merging_stats_from_scaled_array(arr)
    plotter = ResolutionPlotsAndStats(norm_stats, anom_stats)
    resolution_plots = plotter.make_all_plots()
    cc_one_half_data.append(resolution_plots["cc_one_half"])
    i_over_sig_data.append(resolution_plots["i_over_sig_i"])
    del experiments
    del reflections
    gc.collect()
    print("added data for %s" % filepath)


data = {"cc_one_half": cc_one_half_data, "i_over_sigma": i_over_sig_data}
with open("incremental_data.json", "w") as f:
    json.dump(data, f, indent=True)
