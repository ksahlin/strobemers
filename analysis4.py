import os,sys
import argparse

import random
from collections import defaultdict

from modules import indexing, help_functions

def print_stats(acc, datastructure, all_mers):
    abundances = list(all_mers.values())
    del all_mers
    unique = abundances.count(1)
    mean = sum(abundances)/len(abundances)
    ab_sorted = sorted(abundances)
    lower_75 = ab_sorted[1*len(abundances)//4]
    median = ab_sorted[len(abundances)//2]
    upper_75 = ab_sorted[3*len(abundances)//4]
    data = " & ".join([str(d) for d in [datastructure, acc, mean, median, lower_75, upper_75, unique]])
    print(data)



def main(args):
    # cluster ONT sirv sequences?
    k_size = 12
    N_1 = 25
    N_2 = 25 
    genome = {acc: seq for (acc, (seq, _)) in help_functions.readfq(open(args.fasta, 'r'))}
    for acc,seq in genome.items():
        genome[acc] = seq.replace("N", "") # remove Ns
    print(len(genome), sum([len(v) for k,v in genome.items()]))

    print("datastructure & acc & mean & median & lower_75 & upper_75 & number unique")
    

    if args.kmers:
        datastructure = "kmers"
        # s = {p : genome[i:i+k_size] for p, i in enumerate(range(len(genome) - k_size +1))}
        for acc, seq in genome.items():
            all_mers = defaultdict(int)
            for i in range(len(seq) - k_size +1):
                all_mers[hash(seq[i:i+k_size])] += 1
            print_stats(acc, datastructure, all_mers)

    elif args.minstrobes2:
        datastructure = "minstrobes2"
        for acc, seq in genome.items():
            all_mers = defaultdict(int)
            for i,s in enumerate(indexing.minstrobes_iter(seq, k_size, order = 2, N_1 = 50 )):
                all_mers[hash(s)] += 1
            print_stats(acc, datastructure, all_mers)

    elif args.minstrobes3:
        datastructure = "minstrobes3"
        for acc, seq in genome.items():
            all_mers = defaultdict(int)
            for i, s in enumerate(indexing.minstrobes_iter(seq, k_size, order = 3, N_1 = N_1, N_2 = N_2 )):
                all_mers[hash(s)] += 1
            print_stats(acc, datastructure, all_mers)

    elif args.randstrobes2:
        datastructure = "randstrobes2"
        for acc, seq in genome.items():
            all_mers = defaultdict(int)
            for s in indexing.randstrobes_iter(seq, k_size, order = 2, N_1 = 50 ):
                all_mers[s] += 1
            print_stats(acc, datastructure, all_mers)

    elif args.randstrobes3:
        datastructure = "randstrobes3"
        for acc, seq in genome.items():
            all_mers = defaultdict(int)
            for s in indexing.randstrobes_iter(seq, k_size, order = 3, N_1 = N_1, N_2 = N_2 ):
                all_mers[s] += 1
            print_stats(acc, datastructure, all_mers)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calc identity", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--fasta', type=str,  default=False, help='Path to consensus fastq file(s)')
    parser.add_argument('--kmers',  action="store_true", help='Kmer size')
    parser.add_argument('--minstrobes2', action="store_true", help='Kmer size')
    parser.add_argument('--minstrobes3',  action="store_true", help='Kmer size')
    parser.add_argument('--randstrobes2',  action="store_true", help='Kmer size')
    parser.add_argument('--randstrobes3',  action="store_true", help='Kmer size')
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