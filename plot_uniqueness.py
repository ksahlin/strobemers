import os,sys
import argparse

import random

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except (ImportError, RuntimeError):
    print("COULD not import matplotlib")

import seaborn as sns
import pandas as pd
from matplotlib import pyplot


def plot(input_csv, outfolder, acc):
    sns.set(font_scale=1.5)
    indata = pd.read_csv(input_csv)

    dashes = { "kmers" : "", 
                "spaced_sparse": (5,5),
                "spaced_dense": (5,5),
                 "minstrobes2" : (1,1),
                 "minstrobes3" : (1,1),
                 "randstrobes2" : (1,1),
                 "randstrobes3" : (1,1),
                 "hybridstrobes2" : (1,1),
                 "hybridstrobes3" : (1,1)}
    print(indata)
    g = sns.relplot(
        data=indata, x="k", y="unique",
        col="chr", hue="datastructure", style="datastructure", kind="line",  dashes = dashes,
        col_wrap=3, col_order=["chr1", "chr2", "chr3"])
    # ax = sns.lineplot(data=indata, x="k", y="unique", hue="datastructure", style="chr", palette = sns.color_palette()[:7])
    axes = g.axes
    g.set_axis_labels("k", "% unique")
    # g.set_xticklabels([18,24,30,36])
    # ax.set_ylabel("% unique")
    # ax.set_xlabel("k")
    # axes.set_xticks([18,24,30,36] )
    # ax.set_ylim((75, 100))
    g.set(ylim=(80, 100), xticks=[18,24,30,36])
    # ax.set_xticks([18,24,30,36])

    axes = g.axes.flatten()
    axes[0].set_title("Chr1")
    axes[1].set_title("Chr2")
    axes[2].set_title("Chr3")
    # plt.tight_layout()
    # g.set_xticklabels()
    # plt.xlabel(fontsize=14)
    # plt.ylabel(fontsize=14)
    plt.savefig(os.path.join(outfolder, "uniqueness_{0}.eps".format(acc)))
    plt.savefig(os.path.join(outfolder, "uniqueness_{0}.pdf".format(acc)))
    plt.close()




def main(args):

    plot(args.csv, args.outfolder, args.acc)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calc identity", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('csv', type=str, help='Path to consensus fastq file(s)')
    parser.add_argument('outfolder', type=str,  help='A fasta file with transcripts that are shared between samples and have perfect illumina support.')
    parser.add_argument('acc', type=str,  help='A fasta file with transcripts that are shared between samples and have perfect illumina support.')
    # parser.add_argument('--pickled_subreads', type=str, help='Path to an already parsed subreads file in pickle format')
    # parser.set_defaults(which='main')
    args = parser.parse_args()



    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()

    main(args)