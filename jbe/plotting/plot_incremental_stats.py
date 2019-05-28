import json
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
matplotlib.use("Agg")

with open("incremental_data.json") as f:
    data = json.load(f)
    cc_one_half_data = data["cc_one_half"]
    i_over_sig_data = data["i_over_sigma"]

    fontsize = 10

    fig = plt.figure(figsize=(10, 4.5))
    gs = GridSpec(1, 2)
    ax1 = fig.add_subplot(gs[0, 0])
    for i, d in enumerate(cc_one_half_data):
        ax1.plot(d["data"][0]["x"], d["data"][0]["y"], label=r"CC$_{1/2}$"+" (%s)" % i)
        ax1.plot(d["data"][2]["x"], d["data"][2]["y"], label=r"CC$_{1/2}$-anom."+" (%s)" % i)
    ax1.set_xticks(cc_one_half_data[-1]["layout"]["xaxis"]["tickvals"])
    ax1.set_xticklabels(cc_one_half_data[-1]["layout"]["xaxis"]["ticktext"])
    ax1.legend(fontsize=fontsize)
    ax1.set_xlabel(r"Resolution ($\AA$)", fontsize=fontsize)
    ax1.set_ylabel(r"CC$_{1/2}$", fontsize=fontsize)

    ax2 = fig.add_subplot(gs[0, 1])
    for i, d in enumerate(i_over_sig_data):
        ax2.plot(d["data"][0]["x"], d["data"][0]["y"], label="I / sigma (%s)" % i)
    ax2.set_xticks(i_over_sig_data[-1]["layout"]["xaxis"]["tickvals"])
    ax2.set_xticklabels(i_over_sig_data[-1]["layout"]["xaxis"]["ticktext"])
    ax2.set_xlabel(r"Resolution ($\AA$)", fontsize=fontsize)
    ax2.set_ylabel(r"< I / $\sigma$ >", fontsize=fontsize)
    ax2.legend(fontsize=fontsize)

    plt.savefig("incremental.pdf")
