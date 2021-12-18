import os,sys
import argparse
import errno

import random
from collections import defaultdict, Counter
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except (ImportError, RuntimeError):
    print("COULD not import matplotlib")


def mkdir_p(path):
    try:
        os.makedirs(path)
        print("creating", path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def plot_histogram(x, outfolder):
    # plt.bar(x.keys(), x.values())
    # plt.hist(x, density=True, bins=30)  # density=False would make counts
    plt.hist(x, density=False, bins=50) 
    # plt.ylabel('Probability')
    plt.ylabel('Count')
    plt.xlabel('Data')
    plt.xlim(0, 50)
    plt.yscale('log')
    outfile = os.path.join(outfolder, "sampled_position_distribution.pdf")
    plt.savefig(outfile)


# def plot_histogram(input_csv, outfolder):
#     pd.set_option("display.precision", 8)
#     df = pd.read_csv(input_csv)

#     # df_corr = df.loc[df['read_type'] == 'corrected']
#     # df_orig = df.loc[df['read_type'] == 'original']
#     # error_rate_orig = df_orig['error_rate'].tolist()
#     # error_rate_corr = df_corr['error_rate'].tolist()

#     sampled_pos = df_orig['error_rate'].tolist()
#     # bins = [0.1*i for i in range(300)]
#     pyplot.hist(error_rate_corr, 100, range=[0, 20], alpha=0.5, label='Corrected')
#     pyplot.hist(error_rate_orig, 100, range=[0, 20], alpha=0.5, label='Original')
#     pyplot.legend(loc='upper right')
#     # pyplot.xlabel("Difference to genome (%)")
#     pyplot.xlabel("Error rate (%)")
#     pyplot.ylabel("Read count")
#     plt.savefig(os.path.join(outfolder, dataset+ "_full.eps"))
#     plt.savefig(os.path.join(outfolder, dataset+ "_full.pdf"))
#     plt.close()

def main(args):

    sampling_distribution = {}
    p2_sampled = []
    for line in open(args.positions, "r"):
        # nohash,Sahlin1,contig-120_0,75,152
        h, l, ref, p1, p2 = line.strip().split(",")
        p2_sampled.append(p2)

    C = Counter(p2_sampled)
    unique_p2_sampled = len(C)
    most_repetitive = C.most_common(1)
    print(h, l, unique_p2_sampled, most_repetitive)
    plot_histogram(C.values(), args.outfolder)
    # plot_histogram(args.positions, args.outfolder)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evaluate sampling dispersity", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('positions', type=str,  default=False, help='Input positions file')
    parser.add_argument('outfolder', type=str,  default=False, help='Input positions file')

    args = parser.parse_args()
    mkdir_p(args.outfolder)

    main(args)


