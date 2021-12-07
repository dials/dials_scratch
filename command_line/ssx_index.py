"""
Simple script to run the stills indexer on the spotfinding results from a
still sequence i.e. SSX data.

Saves the indexed result into a single experiment list/reflection table
with a joint detector and beam model.
"""

import logging
import time
import sys, os, json
import math
import concurrent.futures
from collections import defaultdict

from libtbx import phil, Auto

from dials.util import log, show_mail_handle_errors, tabulate
from dials.util.options import OptionParser, reflections_and_experiments_from_files
from dials.util.version import dials_version

from dxtbx.model import ExperimentList, Experiment
from dials.array_family import flex
from dials.algorithms.indexing.indexer import Indexer
from dials.algorithms.indexing.max_cell import find_max_cell
from dials.algorithms.clustering import plots as cluster_plotter

from xfel.clustering.cluster import Cluster
from xfel.clustering.cluster_groups import unit_cell_info
from cctbx import crystal
from dials.command_line.combine_experiments import CombineWithReference
from jinja2 import ChoiceLoader, Environment, PackageLoader
import numpy as np


try:
    from typing import List
except ImportError:
    pass

logger = logging.getLogger("dials")

from dials.command_line.index import phil_scope
from dials.algorithms.indexing import DialsIndexError
from dials.util.log import LoggingContext

RAD2DEG = 180 / math.pi

program_defaults_phil_str = """
indexing {
  method = fft1d
}
output.log = dials.ssx_index.log
refinement {
  parameterisation {
    auto_reduction {
      min_nref_per_parameter = 1
      action = fix
    }
    beam.fix = all
    detector.fix = all
  }
  reflections {
    weighting_strategy.override = stills
    outlier.algorithm = null
  }
}
"""

phil_scope = phil.parse(
    """
method = *fft1d *real_space_grid_search
    .type = choice(multi=True)
output.html = dials.ssx_index.html
    .type = str
output.json = None
    .type = str
include scope dials.command_line.index.phil_scope
""",
    process_includes=True,
).fetch(phil.parse(program_defaults_phil_str))

loggers_to_disable = [
    "dials.algorithms.refinement.reflection_processor",
    "dials.algorithms.refinement.refiner",
    "dials.algorithms.refinement.reflection_manager",
    "dials.algorithms.indexing.stills_indexer",
    "dials.algorithms.indexing.nave_parameters",
    "dials.algorithms.indexing.basis_vector_search.real_space_grid_search",
    "dials.algorithms.indexing.basis_vector_search.combinations",
    "dials.algorithms.indexing.indexer",
]

def _index_one(experiment, refl, params, method_list, expt_no):
    for name in loggers_to_disable:
        logging.getLogger(name).disabled = True
    elist = ExperimentList([experiment])
    for method in method_list:
        params.indexing.method = method
        idxr = Indexer.from_parameters(refl, elist, params=params)
        try:
            idxr.index()
        except (DialsIndexError, AssertionError) as e:
            logger.info(
                f"Image {expt_no+1}: Failed to index with {method} method, error: {e}"
            )
            if method == method_list[-1]:
                return None, None
        else:
            logger.info(
                f"Image {expt_no+1}: Indexed {idxr.refined_reflections.size()}/{refl.size()} spots with {method} method."
            )
            return idxr.refined_experiments, idxr.refined_reflections


def generate_plots(summary_data):
    # n_indexed_arrays are cumulative n_indexed for nth lattice
    n_indexed_arrays = [np.zeros(len(summary_data))]
    rmsd_x_arrays = [np.zeros(len(summary_data))]
    rmsd_y_arrays = [np.zeros(len(summary_data))]
    rmsd_z_arrays = [np.zeros(len(summary_data))]
    n_total_indexed = np.zeros(len(summary_data))
    n_strong_array = np.zeros(len(summary_data))
    images = np.arange(1, len(summary_data) + 1)
    n_lattices = 1

    for k in sorted(summary_data.keys()):
        n_lattices_this = len(summary_data[k])
        n_strong_array[k] = summary_data[k][0]["n_strong"]
        for j, cryst in enumerate(summary_data[k]):
            if not cryst["n_indexed"]:
                continue
            if n_lattices_this > n_lattices:
                for _ in range(n_lattices_this - n_lattices):
                    n_indexed_arrays.append(np.zeros(len(summary_data)))
                    rmsd_x_arrays.append(np.zeros(len(summary_data)))
                    rmsd_y_arrays.append(np.zeros(len(summary_data)))
                    rmsd_z_arrays.append(np.zeros(len(summary_data)))
                n_lattices = n_lattices_this
            n_indexed_arrays[j][k] = cryst["n_indexed"]
            rmsd_x_arrays[j][k] = cryst["RMSD_X"]
            rmsd_y_arrays[j][k] = cryst["RMSD_Y"]
            rmsd_z_arrays[j][k] = cryst["RMSD_dPsi"]
            n_total_indexed[k] += cryst["n_indexed"]

    n_indexed_data = [
        {
            "x": images.tolist(),
            "y": n_indexed_arrays[0].tolist(),
            "type": "scatter",
            "mode": "markers",
            "name": "N indexed",
        },
    ]
    rmsd_data = [
        {
            "x": images[rmsd_x_arrays[0] > 0].tolist(),
            "y": rmsd_x_arrays[0][rmsd_x_arrays[0] > 0].tolist(),
            "type": "scatter",
            "mode": "markers",
            "name": "RMSD X",
        },
        {
            "x": images[rmsd_y_arrays[0] > 0].tolist(),
            "y": rmsd_y_arrays[0][rmsd_y_arrays[0] > 0].tolist(),
            "type": "scatter",
            "mode": "markers",
            "name": "RMSD Y",
        },
    ]
    rmsdz_data = [
        {
            "x": images[rmsd_z_arrays[0] > 0].tolist(),
            "y": rmsd_z_arrays[0][rmsd_z_arrays[0] > 0].tolist(),
            "type": "scatter",
            "mode": "markers",
            "name": "RMSD dPsi",
        },
    ]
    if n_lattices > 1:
        n_indexed_data[0]["name"] += " (lattice 1)"
        rmsd_data[0]["name"] += " (lattice 1)"
        rmsd_data[1]["name"] += " (lattice 1)"
        rmsdz_data[0]["name"] += " (lattice 1)"
        for i, arr in enumerate(n_indexed_arrays[1:]):
            sub_images = images[arr > 0]
            sub_data = arr[arr > 0]
            n_indexed_data.append(
                {
                    "x": sub_images.tolist(),
                    "y": sub_data.tolist(),
                    "type": "scatter",
                    "mode": "markers",
                    "name": f"N indexed (lattice {i+2})",
                }
            )
        for i, arr in enumerate(rmsd_x_arrays[1:]):
            sub_images = images[arr > 0]
            sub_data_x = arr[arr > 0]
            sub_data_y = rmsd_y_arrays[i+1][arr > 0]
            rmsd_data.append(
                {
                    "x": sub_images.tolist(),
                    "y": sub_data_x.tolist(),
                    "type": "scatter",
                    "mode": "markers",
                    "name": f"RMSD X (lattice {i+2})",
                },
            )
            rmsd_data.append(
                {
                    "x": sub_images.tolist(),
                    "y": sub_data_y.tolist(),
                    "type": "scatter",
                    "mode": "markers",
                    "name": f"RMSD Y (lattice {i+2})",
                },
            )
        for i, arr in enumerate(rmsd_z_arrays[1:]):
            sub_images = images[arr > 0]
            sub_data = arr[arr > 0]
            rmsdz_data.append(
                {
                    "x": sub_images.tolist(),
                    "y": sub_data.tolist(),
                    "type": "scatter",
                    "mode": "markers",
                    "name": f"RMSD dPsi (lattice {i+2})",
                },
            )
    percent_indexed = 100 * n_total_indexed / n_strong_array
    images = images.tolist()
    n_indexed_data.append(
        {
            "x": images,
            "y": n_strong_array.tolist(),
            "type": "scatter",
            "mode": "markers",
            "name": "N strong",
        },
    )

    percent_bins = np.linspace(0, 100, 51)
    percent_hist = np.histogram(percent_indexed, percent_bins)[0]

    def _generate_hist_data(rmsd_arrays, step=0.01):
        all_rmsd = np.concatenate(rmsd_arrays)
        all_rmsd = all_rmsd[all_rmsd > 0]
        mult = int(1/0.01)
        start = math.floor(np.min(all_rmsd) * mult) / mult
        stop = math.ceil(np.max(all_rmsd) * mult) / mult
        nbins = int((stop - start) / step)
        hist, bin_edges = np.histogram(
            all_rmsd,
            bins=nbins,
            range=(start, stop),
        )
        bin_centers = bin_edges[:-1] + np.diff(bin_edges) / 2
        return hist, bin_centers

    hist_x, bin_centers_x = _generate_hist_data(rmsd_x_arrays)
    hist_y, bin_centers_y = _generate_hist_data(rmsd_y_arrays)
    hist_z, bin_centers_z = _generate_hist_data(rmsd_z_arrays, 0.001)

    plots = {
        "n_indexed": {
            "data": n_indexed_data,
            "layout": {
                "title": "Number of indexed reflections per image",
                "xaxis": {"title": "image number"},
                "yaxis": {"title": "N reflections"},
            },
        },
        "percent_indexed": {
            "data": [{
                "x": images,
                "y": percent_indexed.tolist(),
                "type": "scatter",
                "mode": "markers",
                "name": "Percentage of strong spots indexed",
            }],
            "layout": {
                "title": "Percentage of strong spots indexed per image",
                "xaxis": {"title": "image number"},
                "yaxis": {"title": "Percentage"},
            },
        },
        "percent_indexed_hist":{
            "data": [{
                "x" : percent_bins.tolist(),
                "y" : percent_hist.tolist(),
                "type" : "bar",
            }],
            "layout": {
                "title": "Distribution of percentage indexed",
                "xaxis": {"title": "Percentage indexed"},
                "yaxis": {"title": "Number of images"},
                "bargap": 0,
            },
        },
        "rmsds": {
            "data": rmsd_data,
            "layout": {
                "title": "RMSDs (x, y) per image",
                "xaxis": {"title": "image number"},
                "yaxis": {"title": "RMSD (px)"},
            },
        },
        "rmsdz": {
            "data": rmsdz_data,
            "layout": {
                "title": "RMSD (dPsi) per image",
                "xaxis": {"title": "image number"},
                "yaxis": {"title": "RMSD dPsi (deg)"},
            },
        },
        "rmsdxy_hist":{
            "data": [
                {
                    "x": bin_centers_x.tolist(),
                    "y": hist_x.tolist(),
                    "type" : "bar",
                    "name" : "RMSD X",
                    "opacity": 0.6,
                },
                {
                    "x": bin_centers_y.tolist(),
                    "y": hist_y.tolist(),
                    "type" : "bar",
                    "name" : "RMSD Y",
                    "opacity": 0.6,
                },
            ],
            "layout": {
                "title": "Distribution of RMSDs (x, y)",
                "xaxis": {"title": "RMSD (px)"},
                "yaxis": {"title": "Number of images"},
                "bargap": 0,
                "barmode" : "overlay",
            },
        },
        "rmsdz_hist":{
            "data": [
                {
                    "x": bin_centers_z.tolist(),
                    "y": hist_z.tolist(),
                    "type" : "bar",
                    "name" : "RMSD dPsi",
                },
            ],
            "layout": {
                "title": "Distribution of RMSDs (dPsi)",
                "xaxis": {"title": "RMSD dPsi (deg)"},
                "yaxis": {"title": "Number of images"},
                "bargap": 0,
            },
        },
    }
    return plots


def generate_html_report(plots, filename):
    loader = ChoiceLoader(
        [
            PackageLoader("dials", "templates"),
            PackageLoader("dials", "static", encoding="utf-8"),
        ]
    )
    env = Environment(loader=loader)
    template = env.get_template("simple_report.html")
    html = template.render(
        page_title="DIALS SSX indexing report",
        panel_title="Indexing plots",
        graphs=plots,
    )
    with open(filename, "wb") as f:
        f.write(html.encode("utf-8", "xmlcharrefreplace"))


def index_all_concurrent(experiments, reflections, params, method_list):

    # first determine n_strong per image:
    n_strong_per_image = {}
    for i, (expt, table) in enumerate(zip(experiments, reflections)):
        img = expt.imageset.get_path(i).split("/")[-1]
        n_strong = table.get_flags(table.flags.strong).count(True)
        n_strong_per_image[img] = n_strong

    with concurrent.futures.ProcessPoolExecutor(
        max_workers=params.indexing.nproc
    ) as pool:
        sys.stdout = open(os.devnull, "w")  # block printing from rstbx
        futures = {
            pool.submit(_index_one, expt, table, params, method_list, i): i
            for i, (table, expt) in enumerate(zip(reflections, experiments))
        }
        tables_list = [None] * len(reflections)
        expts_list = [None] * len(reflections)

        for future in concurrent.futures.as_completed(futures):
            try:
                expts, refls = future.result()
                j = futures[future]
            except Exception as e:
                logger.info(e)
            else:
                if refls and expts:
                    tables_list[j] = refls
                    elist = ExperimentList()
                    for jexpt in expts:
                        elist.append(
                            Experiment(
                                identifier=jexpt.identifier,
                                beam=jexpt.beam,
                                detector=jexpt.detector,
                                scan=jexpt.scan,
                                goniometer=jexpt.goniometer,
                                crystal=jexpt.crystal,
                                imageset=jexpt.imageset[j : j + 1],
                            )
                        )
                    expts_list[j] = elist

    sys.stdout = sys.__stdout__  # restore printing

    # now postprocess - record generic information
    results_summary = defaultdict(list)
    indexed_experiments = ExperimentList()
    indexed_reflections = flex.reflection_table()
    n_tot = 0
    for idx, (elist, table) in enumerate(zip(expts_list, tables_list)):
        if not (elist and table):
            img = experiments[idx].imageset.get_path(idx).split("/")[-1]
            n_strong = n_strong_per_image[img]
            results_summary[idx].append({
                   "Image" : img,
                    "n_indexed" : 0,
                    "n_strong" : n_strong,
            })
            continue
        path = elist[0].imageset.get_path(0)
        img = path.split("/")[-1]
        n_strong = n_strong_per_image[img]
        indexed_experiments.extend(elist)
        ids_map = dict(table.experiment_identifiers())
        for k in table.experiment_identifiers().keys():
            del table.experiment_identifiers()[k]
        table["id"] += n_tot
        for k, v in ids_map.items():
            table.experiment_identifiers()[k + n_tot] = v
        n_tot += len(ids_map.keys())
        indexed_reflections.extend(table)

        # record some things for printing to output log/html
        for id_ in table.experiment_identifiers().keys():
            selr = table.select(table["id"] == id_)
            calx, caly, _ = selr["xyzcal.px"].parts()
            obsx, obsy, _ = selr["xyzobs.px.value"].parts()
            delpsi = selr["delpsical.rad"]
            rmsd_x = flex.mean((calx - obsx) ** 2) ** 0.5
            rmsd_y = flex.mean((caly - obsy) ** 2) ** 0.5
            rmsd_z = flex.mean(((delpsi) * RAD2DEG) ** 2) ** 0.5
            n_id_ = calx.size()
            results_summary[idx].append(
                {
                   "Image" : img,
                    "n_indexed" : n_id_,
                    "n_strong" : n_strong,
                    "RMSD_X" : rmsd_x,
                    "RMSD_Y" : rmsd_y,
                    "RMSD_dPsi" : rmsd_z,
                }
            )

    indexed_reflections.assert_experiment_identifiers_are_consistent(
        indexed_experiments
    )

    return indexed_experiments, indexed_reflections, results_summary

def make_summary_table(results_summary):
    # make a summary table
    overall_summary_header = [
        "Image",
        "expt_id",
        "n_indexed",
        "RMSD X",
        "RMSD Y",
        "RMSD dPsi",
    ]

    rows = []
    total = 0
    if any(len(v) > 1 for v in results_summary.values()):
        show_lattices = True
        overall_summary_header.insert(1, "lattice")
    else:
        show_lattices = False
    for k in sorted(results_summary.keys()):
        for j, cryst in enumerate(results_summary[k]):
            if not cryst["n_indexed"]:
                continue
            n_idx, n_strong = (cryst["n_indexed"], cryst["n_strong"])
            frac_idx =  f"{n_idx}/{n_strong} ({100*n_idx/n_strong:2.1f}%)"
            row = [cryst["Image"], str(total), frac_idx, cryst["RMSD_X"], cryst["RMSD_Y"], cryst["RMSD_dPsi"]]
            if show_lattices:
                row.insert(1, j + 1)
            rows.append(row)
            total += 1

    summary_table = tabulate(rows, overall_summary_header)
    return summary_table

def make_cluster_plots(large_clusters):
    cluster_plots = {}
    for n, cluster in enumerate(large_clusters):
        uc_params = [flex.double() for i in range(6)]
        for c in cluster.members:
            ucp = c.crystal_symmetry.unit_cell().parameters()
            for i in range(6):
                uc_params[i].append(ucp[i])
        d_this = cluster_plotter.plot_uc_histograms(uc_params)
        d_this["uc_scatter"]["layout"]["title"] += f" cluster {n+1}"
        d_this["uc_hist"]["layout"]["title"] += f" cluster {n+1}"
        d_this[f"uc_scatter_{n}"] = d_this.pop("uc_scatter")
        d_this[f"uc_hist_{n}"] = d_this.pop("uc_hist")
        cluster_plots.update(d_this)
    return cluster_plots


def index(experiments, observed, params):
    params.refinement.parameterisation.scan_varying = False
    params.indexing.stills.indexer = "stills"

    reflections = observed.split_by_experiment_id()
    # Calculate necessary quantities
    for refl, experiment in zip(reflections, experiments):
        elist = ExperimentList([experiment])
        refl["imageset_id"] = flex.int(refl.size(), 0)  # needed for centroid_px_to_mm
        refl.centroid_px_to_mm(elist)
        refl.map_centroids_to_reciprocal_space(elist)

    if (params.indexing.max_cell is Auto) and (
        not params.indexing.known_symmetry.unit_cell
    ):
        max_cells = []
        for refl in reflections:
            try:
                max_cells.append(find_max_cell(refl).max_cell)
            except (DialsIndexError, AssertionError):
                pass
        logger.info(f"Setting max cell to {max(max_cells):.1f} " + "\u212B")
        params.indexing.max_cell = max(max_cells)

    method_list = params.method
    if "real_space_grid_search" in method_list:
        if not params.indexing.known_symmetry.unit_cell:
            logger.info("No unit cell given, real_space_grid_search will not be used")
            method_list.remove("real_space_grid_search")
    methods = ", ".join(method_list)
    pl = "s" if (len(method_list) > 1) else ""
    logger.info(f"Attempting indexing with {methods} method{pl}")

    indexed_experiments, indexed_reflections, results_summary = index_all_concurrent(
        experiments,
        reflections,
        params,
        method_list,
    )

    summary_table = make_summary_table(results_summary)
    logger.info("\nSummary of images sucessfully indexed\n" + summary_table)

    n_images = len(set(e.imageset.get_path(0) for e in indexed_experiments))
    logger.info(f"{indexed_reflections.size()} spots indexed on {n_images} images\n")

    # combine beam and detector models if not already
    if (len(indexed_experiments.detectors())) > 1 or (
        len(indexed_experiments.beams())
    ) > 1:
        combine = CombineWithReference(
            detector=indexed_experiments[0].detector, beam=indexed_experiments[0].beam
        )
        elist = ExperimentList()
        for expt in indexed_experiments:
            elist.append(combine(expt))
        indexed_experiments = elist

    return indexed_experiments, indexed_reflections, results_summary


@show_mail_handle_errors()
def run(args: List[str] = None, phil: phil.scope = phil_scope) -> None:
    parser = OptionParser(
        usage="",
        read_experiments=True,
        read_reflections=True,
        phil=phil_scope,
        check_format=False,
        epilog=__doc__,
    )
    params, _ = parser.parse_args(args=args, show_diff_phil=False)

    if not params.input.experiments or not params.input.reflections:
        parser.print_help()
        sys.exit()

    reflections, experiments = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )
    log.config(verbosity=1, logfile=params.output.log)
    logger.info(dials_version())

    diff_phil = parser.diff_phil.as_str()
    if diff_phil:
        logger.info("The following parameters have been modified:\n%s", diff_phil)

    st = time.time()
    indexed_experiments, indexed_reflections, summary_data = index(experiments, reflections[0], params)

    # print some clustering information
    ucs = Cluster.from_crystal_symmetries([
        crystal.symmetry(
            unit_cell=expt.crystal.get_unit_cell(),
            space_group=expt.crystal.get_space_group(),
        )
        for expt in indexed_experiments
    ])
    clusters, _ = ucs.ab_cluster(
        5000, log=None, write_file_lists=False, doplot=False
    )
    large_clusters = []
    cluster_plots = {}
    threshold = math.floor(0.05 * len(indexed_experiments))
    for cluster in clusters:
        if len(cluster.members) > threshold:
            large_clusters.append(cluster)
    large_clusters.sort(key=lambda x:len(x.members), reverse=True)

    if large_clusters:
        logger.info(f"""
Unit cell clustering analysis, clusters with >5% of the number of crystals indexed
""" + unit_cell_info(large_clusters))
        if params.output.html or params.output.json:
            cluster_plots = make_cluster_plots(large_clusters)
    else:
        logger.info(f"No clusters found with >5% of the number of crystals.")

    logger.info(f"Saving indexed experiments to {params.output.experiments}")
    indexed_experiments.as_file(params.output.experiments)
    logger.info(f"Saving indexed reflections to {params.output.reflections}")
    indexed_reflections.as_file(params.output.reflections)

    if params.output.html or params.output.json:
        summary_plots = generate_plots(summary_data)
        if cluster_plots:
            summary_plots.update(cluster_plots)
        if params.output.html:
            generate_html_report(summary_plots, params.output.html)
        if params.output.json:
            with open(params.output.json, "w") as outfile:
                json.dump(summary_plots, outfile)

    logger.info(f"Total time: {time.time() - st:.2f}s")


if __name__ == "__main__":
    run()
