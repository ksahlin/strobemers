import os,sys
import argparse

import random
from collections import defaultdict

from modules import indexing, help_functions

def print_stats(acc, datastructure, all_mers, k_size, total_mers):
    abundances = list(all_mers.values())
    del all_mers
    unique = abundances.count(1)
    percent_unique = round(100*unique/total_mers, 1)
    mean = sum(abundances)/len(abundances)
    ab_sorted = sorted(abundances)
    lower_75 = ab_sorted[1*len(abundances)//4]
    median = ab_sorted[len(abundances)//2]
    upper_75 = ab_sorted[3*len(abundances)//4]
    data = ",".join([str(d) for d in [k_size, datastructure, acc, mean, median, lower_75, upper_75, percent_unique]])
    print(data)


def compute_uniqueness(args, acc, seq, k_size, total_mers):
    # w_1 = 25
    # w_2 = 25 
    w = 1
    w_low = 25
    w_high = 50
    all_mers = defaultdict(int)
    if args.kmers:
        datastructure = "kmers"
        for i in range(len(seq) - k_size +1):
            all_mers[hash(seq[i:i+k_size])] += 1

    if args.spaced_dense:
        datastructure = "spaced_dense"
        span_size = k_size+k_size//2
        positions = set(random.sample(range(1, span_size - 1 ), k_size-2)) 
        positions.add(0)
        positions.add(span_size - 1) # asserts first and last position is sampled so that we have a spaced kmer of length span size
        for s in indexing.spaced_kmers_iter(seq, k_size, span_size, positions):
            all_mers[s] += 1

    if args.spaced_sparse:
        datastructure = "spaced_sparse"
        span_size = 3*k_size
        positions = set(random.sample(range(1, span_size - 1 ), k_size-2)) 
        positions.add(0)
        positions.add(span_size - 1) # asserts first and last position is sampled so that we have a spaced kmer of length span size
        for s in indexing.spaced_kmers_iter(seq, k_size, span_size, positions):
            all_mers[s] += 1

    elif args.minstrobes2:
        datastructure = "minstrobes2"
        for s in indexing.minstrobes_iter(seq, k_size, w_low, w_high, w, order = 2, buffer_size = 10000000):
            all_mers[s] += 1

    elif args.minstrobes3:
        datastructure = "minstrobes3"
        for s in indexing.minstrobes_iter(seq, k_size, w_low, w_high, w, order = 3, buffer_size = 10000000):
            all_mers[s] += 1

    elif args.randstrobes2:
        datastructure = "randstrobes2"
        for s in indexing.randstrobes_iter(seq, k_size, w_low, w_high, w, order = 2, buffer_size = 10000000): # (seq, k_size, order = 2, w_1 = 50 ):
            all_mers[s] += 1

    elif args.randstrobes3:
        datastructure = "randstrobes3"
        for s in indexing.randstrobes_iter(seq, k_size, w_low, w_high, w, order = 3, buffer_size = 10000000):
            all_mers[s] += 1


    elif args.hybridstrobes2:
        datastructure = "hybridstrobes2"
        for s in indexing.hybridstrobes_iter(seq, k_size, w_low, w_high, w, order = 2, buffer_size = 10000000): # (seq, k_size, order = 2, w_1 = 50 ):
            all_mers[s] += 1

    elif args.hybridstrobes3:
        datastructure = "hybridstrobes3"
        for s in indexing.hybridstrobes_iter(seq, k_size, w_low, w_high, w, order = 3, buffer_size = 10000000):
            all_mers[s] += 1


    print_stats(acc, datastructure, all_mers, k_size, total_mers)



def main(args):

    genome = {acc: seq for (acc, (seq, _)) in help_functions.readfq(open(args.fasta, 'r'))}

    for acc,seq in genome.items():
        acc = acc.split()[0]
        # print(acc)
        genome[acc] = seq.replace("N", "") # remove Ns

    print(len(genome), sum([len(v) for k,v in genome.items()]))

    print("datastructure,k,acc,mean,median,lower_75,upper_75,\%-unique")

    total_mers = {}

    for acc, seq in genome.items():
        # print(acc)
        for k_size in [18,24,30,36]: #[18,24,30]:
            total_mers[acc] = len(seq) - k_size + 1
            if acc == "chr1" or acc == "chr2" or acc == "chr3": # or acc == "chr4" or acc == "chr5":
                compute_uniqueness(args, acc, seq, k_size, total_mers[acc])







if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calc identity", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--fasta', type=str,  default=False, help='Path to consensus fastq file(s)')
    parser.add_argument('--kmers',  action="store_true", help='Kmer size')
    parser.add_argument('--minstrobes2', action="store_true", help='Kmer size')
    parser.add_argument('--minstrobes3',  action="store_true", help='Kmer size')
    parser.add_argument('--randstrobes2',  action="store_true", help='Kmer size')
    parser.add_argument('--randstrobes3',  action="store_true", help='Kmer size')
    parser.add_argument('--hybridstrobes2',  action="store_true", help='Kmer size')
    parser.add_argument('--hybridstrobes3',  action="store_true", help='Kmer size')
    parser.add_argument('--spaced_dense',  action="store_true", help='Kmer size')
    parser.add_argument('--spaced_sparse',  action="store_true", help='Kmer size')
    # parser.add_argument('--k', type=int, default=13, help='Kmer size')
    # parser.add_argument('--w', type=int, default=20, help='Window size')
    # parser.add_argument('--outfolder', type=str,  default=None, help='A fasta file with transcripts that are shared between samples and have perfect illumina support.')
    # parser.add_argument('--pickled_subreads', type=str, help='Path to an already parsed subreads file in pickle format')
    # parser.set_defaults(which='main')
    args = parser.parse_args()



    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()

    main(args)