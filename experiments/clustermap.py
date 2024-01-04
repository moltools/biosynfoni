import sys, os, argparse

import numpy as np
import seaborn as sns
import scipy.cluster.hierarchy as sch
import scipy.cluster.vq as scv
import matplotlib as mpl

# matplotlib.use('Agg')       #if in background
import matplotlib.pylab as plt
from matplotlib.patches import Patch
import pandas as pd

from biosynfoni.def_biosynfoni import FP_VERSIONS


def cli():
    parser = argparse.ArgumentParser()

    parser.add_argument("fingerprint", help="Fingerprint file")
    parser.add_argument("labels", help="Labels file")
    parser.add_argument(
        "-s",
        "--subsample",
        required=False,
        type=int,
        help="subsample size",
        default=None,
    )
    parser.add_argument(
        "-r",
        "--seed",
        "--randomseed",
        required=False,
        type=int,
        help="seed for subsampling",
        default=None,
    )
    parser.add_argument(
        "-S",
        "--synthetic compounds",
        required=False,
        type=str,
        help="format of output file",
        default=None,
    )
    return parser.parse_args()


class recursion_depth:
    def __init__(self, limit):
        self.limit = limit
        self.default_limit = sys.getrecursionlimit()

    def __enter__(self):
        sys.setrecursionlimit(self.limit)

    def __exit__(self, type, value, traceback):
        sys.setrecursionlimit(self.default_limit)


class ClusterMap:
    def __init__(self, df, labels, metric, method) -> None:
        self.df = df
        self.indexes = df.index.values
        self.labels = labels
        self.colordict = self.get_colordict()
        self.metric = metric
        self.method = method
        # calculations:
        self.distances = self.get_distances()
        self.clustering = self.get_clustering()
        # self.distances = None
        # self.clustering = None
        # self.tree = self.get_tree()
        self.colors, self.handles = self._get_category_colors_handles(self.labels)
        self.clustermap = self.seacluster()
        self.clusterfig = self.get_clusterplot()
        plt.close()
        pass

    def get_distances(self):
        """calculates distance with scipy cluster hierarchy"""
        return sch.distance.pdist(self.df, metric=self.metric)

    def set_distances(self, distances):
        self.distances = distances
        return None

    def get_clustering(self):
        """calculates dendogram from data frame and distances"""
        # plt.title(out_file)
        with recursion_depth(10000):
            clustering = sch.linkage(self.distances, method=self.method)
        # plt.close()
        return clustering

    def get_tree(self):
        return sch.dendrogram(
            self.clustering, leaf_font_size=2, color_threshold=4, labels=self.indexes
        )

    def seacluster(self):
        # comp_color, comp_handles = get_gnsr_diff_color() #compounds

        cmap = _cmap_makezerowhite("mako")

        # self.distances = self.set_distances(self.get_distances())
        # self.clustering = self.get_clustering()
        # fig, ax = plt.subplots(figsize=(15,6))
        # sns.set(font_scale=0.5)
        with recursion_depth(20000):
            seacluster = sns.clustermap(
                self.df,
                method=self.method,
                metric=self.metric,
                xticklabels=1,  # every 1: label
                # robust = True,
                # square = True,
                row_colors=self.colors,
                # col_colors = phyl_color,
                # z_score = 1,
                # figsize= (len(list(data_frame)), min(200, len(data_frame.index.values))),
                cmap=cmap,
                cbar_kws={"shrink": 0.3, "label": "counts"},
            )
        return seacluster

    def get_clusterplot(self, legend_title="NP-Classifier prediction"):
        """returns plt of clustermap with legend of categories"""

        # plt.figure(figsize=(15,6))

        # plt.setp(seacluster.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
        # plt.setp(seacluster.ax_heatmap.xaxis.get_majorticklabels(), rotation=30, fontsize=8)
        # plt.setp(seacluster.ax_heatmap.xaxis.get_majorticklabels(), rotation=-30)
        # seacluster.set_xlim([0,2])

        # seacluster.set_title(TITLE, fontsize=16, fontdict={})
        plt.title(
            f"Hierarchical clustering of compounds and substructures {self.method}, {self.metric}"
        )  # + 'z scored')

        # seacluster = self.seacluster()

        # improve layout
        # plt.tight_layout()

        # add legend
        handles = self.handles
        legend_colors = [Patch(facecolor=handles[name]) for name in handles.keys()]
        plt.legend(
            legend_colors,
            handles.keys(),
            title=legend_title,
            bbox_to_anchor=(1, 1),
            bbox_transform=plt.gcf().transFigure,
            loc="upper right",
        )

        # plt.show()
        # return dendrogram_linkage
        return

    def save_clustermap(self, fmt="pdf"):
        out_file = f"clustermap_{self.method}_{self.metric}.pdf"
        # self.clusterfig.savefig(out_file, format=fmt)
        self.clustermap.savefig(out_file, format=fmt)
        # plt.savefig(out_file, format=fmt)
        return None

    def get_dendogram_tree(self):
        # for cutting
        return self.clustermap.dendrogram_row.dendrogram

    def get_dendogram_linkage(self):
        return self.clustermap.dendrogram_row.linkage

    def get_colordict(self):
        colordict = {
            "Terpenoids": sns.color_palette("Set3")[6],  # green
            "Alkaloids": sns.color_palette("Set3")[9],  # purple
            "Shikimates and Phenylpropanoids": sns.color_palette("Set3")[4],  # blue
            "Fatty acids": sns.color_palette("Set3")[5],  # orange
            "Carbohydrates": sns.color_palette("Set3")[7],  # pink
            "Polyketides": sns.color_palette("Set3")[3],  # light red
            "Amino acids and Peptides": "bisque",
            # "No NP-Classifier prediction": "grey",
            "None": "grey",
            "Synthetic": "black",
        }
        return colordict

    def _get_category_colors_handles(self, categories: pd.Series):
        """uses colour dictionary to assign colors to the categories"""
        network_dict = {}
        categories.fillna("None", inplace=True)

        for ind, cat in categories.items():
            # network_dict[ind] = [self.colordict[str(cat).split(',')[0]]] #in case of multiple categories, take the first one
            network_dict[ind] = [self.colordict[str(cat)]]
        network_colors = pd.DataFrame.from_dict(network_dict, orient="index")
        network_colors.columns = [""]
        handles = self.colordict
        return network_colors, handles


def _cmap_makezerowhite(default_cmap: str = "mako"):
    # define color map:-------------------------------------------------
    cmap = sns.color_palette(default_cmap, as_cmap=True)  # define the colormap
    # extract all colors from the .jet map
    cmaplist = [cmap(i) for i in range(cmap.N)]
    # force the first color entry to be white (to distinguish 0 from low values)
    cmaplist[0] = (1.0, 1.0, 1.0, 0)
    # create the new map
    cmap = mpl.colors.LinearSegmentedColormap.from_list("Custom cmap", cmaplist, cmap.N)
    return cmap


def main():
    blocked = "../../input/coco/coconut_full1103.bsf"
    overlapped = "../../input/coco/coconut_full1103noblock.bsf"

    input_file = blocked
    filetype = "pdf"
    version = input_file.split("/")[-1].split("_")[-1].split(".")[0]
    substructure_names = FP_VERSIONS["full_1103"]

    df = pd.read_csv(blocked, sep=",", header=None, dtype=int)
    df.columns = substructure_names

    npcs = pd.read_csv(
        "../../input/coco/coconut_classes.tsv",
        sep="\t",
        header=None,
        dtype=str,
        usecols=[0],
    )

    # filter out multiple-prediction compounds
    npcs.fillna(",", inplace=True)
    df = df[~npcs[0].str.contains(",")]
    npcs = npcs[~npcs[0].str.contains(",")]
    # oneprediction_df_overlap = df_overlap[~npcs[0].str.contains(',')]

    # filter out only-zero columns in df ~~~~~~~~~~~~~~~ CHECK IF THIS APPLIES TO YOUR PURPOSES ~~~~~~~~~~~~~~~
    df = df.loc[:, (df != 0).any(axis=0)]
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # subsample indexes
    if args.subsample:
        np.random.seed(args.seed)
        idx = np.random.choice(df.index, args.subsample, replace=False)
        df = df.loc[idx]
        npcs = npcs.loc[idx]

    npcs_series = npcs.iloc[:, 0]

    # check if indexes are the same
    assert df[df.index != npcs.index].empty
    assert npcs[npcs.index != npcs_series.index].empty
    assert type(npcs_series) == pd.Series

    # if args.synthetic: #under construction
    #     synthetic_fp = pd.read_csv(args.synthetic, dtype=int)
    #     synthetic_fp = np.random.choice(
    #         synthetic_fp.shape[0], fp.shape[0], replace=False
    #     )
    #     fp = np.concatenate((fp, synthetic_fp))
    #     synthetic_labels = np.array(["synthetic" for _ in range(synthetic_fp.shape[0])])
    #     labels = np.concatenate((labels, synthetic_labels))

    cwd = os.getcwd()
    # make a directory in grandparent directory called clustermaps
    # os.chdir("../../")
    os.makedirs("clustermaps", exist_ok=True)
    os.chdir("clustermaps")

    for method in ["average", "complete", "single", "weighted"]:
        for metric in [
            "euclidean",
            "cityblock",
            "cosine",
            "correlation",
            "hamming",
            "jaccard",
            "mahalanobis",
            "chebyshev",
            "canberra",
            "braycurtis",
            "dice",
            "kulsinski",
            "matching",
            "rogerstanimoto",
            "russellrao",
            "sokalmichener",
            "sokalsneath",
            "yule",
        ]:
            # clustermap = ClusterMap(df, npcs_series, metric, method)
            # clustermap.save_clustermap(fmt=filetype)
            try:
                clustermap = ClusterMap(df, npcs_series, metric, method)
                # clustermap.save_clustermap(fmt=filetype)
            except:
                # errors can occur for some metrics if they have too small sample sets, or with certain combinations
                print(f"failed for {method} and {metric}")
            pass
        pass

    os.chdir(cwd)
    return None


if __name__ == "__main__":
    main()
