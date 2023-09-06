import json
import matplotlib.pyplot as plt

fig,ax = plt.subplots(nrows=1, ncols=2)
dirs = ["xia2-ellipsoids6", "xia2-stills"]
labels = ["ellipsoid", "stills"]
for d, plotlabel in zip(dirs, labels):

    with open(f"{d}/dials.merge.json") as f:
        data = json.load(f)
    keys = data.keys()
    assert len(keys) == 1
    data = data[list(keys)[0]]
    normal_cchalf = data["resolution_plots"]["cc_one_half"]["data"][0]
    xticks=data["resolution_plots"]["cc_one_half"]["layout"]["xaxis"]["tickvals"]
    xtickvals = data["resolution_plots"]["cc_one_half"]["layout"]["xaxis"]["ticktext"]


    ax[0].scatter(normal_cchalf["x"], normal_cchalf["y"], label=plotlabel)
    ax[0].set_xticks(xticks)
    ax[0].set_xticklabels(xtickvals)
    ax[0].set_ylabel("CC1/2")
    ax[0].set_xlabel("Resolution (A)")
    ax[0].legend()

    rsplit = data["resolution_plots"]["r_split"]["data"][0]
    xticks=data["resolution_plots"]["r_split"]["layout"]["xaxis"]["tickvals"]
    xtickvals = data["resolution_plots"]["r_split"]["layout"]["xaxis"]["ticktext"]

    ax[1].scatter(rsplit["x"], rsplit["y"], label=plotlabel)
    ax[1].set_xticks(xticks)
    ax[1].set_xticklabels(xtickvals)
    ax[1].set_ylabel("R-split")
    ax[1].set_xlabel("Resolution (A)")
    ax[1].legend()
plt.savefig("compare_plots.png")