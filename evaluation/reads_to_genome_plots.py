## Various plots from large table

import sys
import argparse
import os
import random
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except (ImportError, RuntimeError):
    print("COULD not import matplotlib")
# import matplotlib.pyplot as plt
# import matplotlib

import numpy as np
import seaborn as sns
import pandas as pd

def plot_coverage(input_csv, outfolder):
    sns.set(font_scale=1.0)
    indata = pd.read_csv(input_csv)

    ax = sns.lmplot(data=indata, x="r_len", y="r_frac_cov", hue="method", legend_out = True)

    plt.xlabel('Read length',fontsize=10)
    plt.ylabel('Fraction covered',fontsize=12)
    # plt.tick_params(rotation=90)
    # ax.set_xticklabels(size = 10)
    # ax.set_xticklabels(ax.get_xmajorticklabels(), fontsize = 8)
    plt.ylim(0, 1)
    # plt.tight_layout()
    # plt.legend(loc='lower left', prop={'size': 6})
    plt.savefig(os.path.join(outfolder, "plot_coverage.eps"))
    plt.savefig(os.path.join(outfolder, "plot_coverage.pdf"))
    plt.clf()
    plt.close()

def plot_nr_hits(input_csv, outfolder):
    sns.set(font_scale=1.0)
    indata = pd.read_csv(input_csv)

    ax = sns.lmplot(data=indata, x="r_len", y="r_hits", hue="method", legend_out = True)

    plt.xlabel('Read length',fontsize=10)
    plt.ylabel('Number of NAMs',fontsize=12)
    # plt.tight_layout()
    plt.savefig(os.path.join(outfolder, "plot_nr_hits.eps"))
    plt.savefig(os.path.join(outfolder, "plot_nr_hits.pdf"))
    plt.clf()
    plt.close()



def main(args):
    
    sns.set(style="whitegrid")
    flatui = ["#2ecc71", "#e74c3c"] # https://chrisalbon.com/python/data_visualization/seaborn_color_palettes/
    sns.set_palette(flatui)    # total_error_rate(args.input_csv, args.outfolder)
    plot_coverage(args.results, args.outfolder)
    plot_nr_hits(args.results, args.outfolder)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evalueate randstrobes.")
    parser.add_argument('results', type=str, help='Path to coverage stats file')
    parser.add_argument('outfolder', type=str, help='Path to all stats file')

    args = parser.parse_args()

    outfolder = args.outfolder
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
    main(args)