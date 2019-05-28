import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
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
    experiments = load.experiment_list(filepath+"scaled_experiments.json", check_format=False)
    reflections = flex.reflection_table.from_pickle(filepath+"scaled.pickle")

    arr = scaled_data_as_miller_array([reflections], experiments)
    norm_stats, anom_stats = merging_stats_from_scaled_array(arr)
    plotter = ResolutionPlotsAndStats(norm_stats, anom_stats)
    resolution_plots = plotter.make_all_plots()
    cc_one_half_data.append(resolution_plots["cc_one_half"]["data"])
    i_over_sig_data.append(resolution_plots["i_over_sig_i"]["data"])

fontsize = 10

fig = plt.figure(figsize=(10, 4.5))
gs = GridSpec(1, 2)
ax1 = fig.add_subplot(gs[0, 0])
for i, d in enumerate(cc_one_half_data):
    ax1.plot(d[0]["x"], d[0]["y"], label=r"CC$_{1/2}$"+" (%s)" % i)
    ax1.plot(d[2]["x"], d[2]["y"], label=r"CC$_{1/2}$-anom."+" (%s)" % i)
ax1.set_xticks(resolution_plots["cc_one_half"]["layout"]["xaxis"]["tickvals"])
ax1.set_xticklabels(resolution_plots["cc_one_half"]["layout"]["xaxis"]["ticktext"])
ax1.legend(fontsize=fontsize)
ax1.set_xlabel(r"Resolution ($\AA$)", fontsize=fontsize)
ax1.set_ylabel(r"CC$_{1/2}$", fontsize=fontsize)

ax2 = fig.add_subplot(gs[0, 1])
for i, d in enumerate(i_over_sig_data):
    ax2.plot(d[0]["x"], d[0]["y"], label="I / sigma (%s)" % i)
ax2.set_xticks(resolution_plots["i_over_sig_i"]["layout"]["xaxis"]["tickvals"])
ax2.set_xticklabels(resolution_plots["i_over_sig_i"]["layout"]["xaxis"]["ticktext"])
ax2.set_xlabel(r"Resolution ($\AA$)", fontsize=fontsize)
ax2.set_ylabel(r"< I / $\sigma$ >", fontsize=fontsize)
ax2.legend(fontsize=fontsize)

plt.savefig("incremental_stats.pdf")
