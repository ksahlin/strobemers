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
    indata = pd.read_csv(input_csv,sep='\t')
    # ax = sns.violinplot(x="day", y="total_bill", hue="smoker",
    #                 data=tips, palette="muted")
    ax = sns.lineplot(x="ref_id", y="coverage", hue="method", ci = "sd", hue_order= ["randstrobes-(3,10,20,70)", "randstrobes-(2,15,20,70)", "hybridstrobes-(3,10,20,70)", "hybridstrobes-(2,15,20,70)", "minstrobes-(3,10,20,70)", "minstrobes-(2,15,20,70)", "kmers"], 
                         data=indata, markers=True)
    # ax = sns.violinplot(x="ref_id", y="coverage", hue="method",
    #                       hue_order= ["randstrobes", "kmers"], data=indata)
    plt.xlabel('SIRV', fontsize=14)
    plt.ylabel('Fraction covered',fontsize=16)
    plt.tick_params(rotation=90)
    # ax.set_xticklabels(size = 10)
    ax.set_xticklabels(ax.get_xmajorticklabels(), fontsize = 8)
    plt.ylim(0, 1)
    plt.tight_layout()
    plt.legend(loc='lower left', prop={'size': 6})
    # ax.set_ylabel("Fraction correct")
    # ax.set_xlabel("Exon size")
    plt.savefig(os.path.join(outfolder, "plot_coverage.eps"))
    plt.savefig(os.path.join(outfolder, "plot_coverage.pdf"))
    plt.clf()
    plt.close()

def plot_nr_hits(input_csv, outfolder):
    sns.set(font_scale=1.0)
    indata = pd.read_csv(input_csv,sep='\t')
    # ax = sns.barplot(x="ref_id", y="nr_hits", hue="method", data=indata,
    #                       hue_order= ["randstrobes", "kmers"])
    ax = sns.lineplot(x="ref_id", y="nr_hits", hue="method", ci = "sd", hue_order= ["randstrobes-(3,10,20,70)", "randstrobes-(2,15,20,70)", "hybridstrobes-(3,10,20,70)", "hybridstrobes-(2,15,20,70)", "minstrobes-(3,10,20,70)", "minstrobes-(2,15,20,70)", "kmers"], 
                          data=indata, markers=True)
    plt.xlabel('SIRV', fontsize=14)
    plt.ylabel('Number NAMs',fontsize=16)
    plt.tick_params(rotation=90)
    # ax.set_xticklabels(size = 10)
    ax.set_xticklabels(ax.get_xmajorticklabels(), fontsize = 8)
    # plt.ylim(0, 1)
    plt.tight_layout()
    # plt.legend(loc='lower right')
    # ax.set_ylabel("Fraction correct")
    # ax.set_xlabel("Exon size")
    plt.savefig(os.path.join(outfolder, "plot_nr_hits.eps"))
    plt.savefig(os.path.join(outfolder, "plot_nr_hits.pdf"))
    plt.clf()
    plt.close()


def plot_normalized_match_length(input_csv, outfolder):
    sns.set(font_scale=1.0)
    indata = pd.read_csv(input_csv,sep='\t')
    # ax = sns.barplot(x="ref_id", y="nr_hits", hue="method", data=indata,
    #                       hue_order= ["randstrobes", "kmers"])
    plt.tick_params(axis='x', which='minor', labelsize=7)
    ax = sns.lineplot(x="ref_id", y="normalized_match_length", hue="method", ci = "sd", hue_order= ["randstrobes-(3,10,20,70)", "randstrobes-(2,15,20,70)", "hybridstrobes-(3,10,20,70)", "hybridstrobes-(2,15,20,70)", "minstrobes-(3,10,20,70)", "minstrobes-(2,15,20,70)", "kmers"], 
                        data=indata, markers=True)
    plt.xlabel('SIRV', fontsize=14)
    plt.ylabel('Normalized NAM length',fontsize=16)
    plt.tick_params(rotation=90)
    # ax.set_xticklabels(size = 10)
    ax.set_xticklabels(ax.get_xmajorticklabels(), fontsize = 8)
    plt.ylim(0, 1)
    plt.tight_layout()
    # plt.legend(loc='lower right')
    # ax.set_ylabel("Fraction correct")
    # ax.set_xlabel("Exon size")
    plt.savefig(os.path.join(outfolder, "plot_normalized_match_length.eps"))
    plt.savefig(os.path.join(outfolder, "plot_normalized_match_length.pdf"))
    plt.clf()
    plt.close()



def main(args):
    
    sns.set(style="whitegrid")
    flatui = ["#2ecc71", "#e74c3c"] # https://chrisalbon.com/python/data_visualization/seaborn_color_palettes/
    sns.set_palette(flatui)    # total_error_rate(args.input_csv, args.outfolder)
    plot_coverage(args.coverage, args.outfolder)
    plot_nr_hits(args.nr_hits, args.outfolder)
    plot_normalized_match_length(args.normalized_match_length, args.outfolder)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evalueate randstrobes.")
    parser.add_argument('coverage', type=str, help='Path to coverage stats file')
    parser.add_argument('nr_hits', type=str, help='Path to nr_hits file')
    parser.add_argument('normalized_match_length', type=str, help='Path to normalized_match_length file')
    parser.add_argument('outfolder', type=str, help='Path to all stats file')

    args = parser.parse_args()

    outfolder = args.outfolder
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
    main(args)