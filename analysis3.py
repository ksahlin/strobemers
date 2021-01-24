import os,sys
import argparse

import random
from collections import defaultdict

from modules import indexing, help_functions


def main(args):
    k_size = 30

    genome = {acc: seq for (acc, (seq, _)) in help_functions.readfq(open(args.fasta, 'r'))}
    for acc,seq in genome.items():
        genome[acc] = seq.replace("N", "") # remove Ns
    print(len(genome), sum([len(v) for k,v in genome.items()]))
    

    kmers = defaultdict(int)
    # s = {p : genome[i:i+k_size] for p, i in enumerate(range(len(genome) - k_size +1))}
    for acc, seq in genome.items():
        for i in range(len(seq) - k_size +1):
            kmers[hash(seq[i:i+k_size])] += 1
    # for k,v in s.items():
    #     kmers[v] += 1
    abundances = list(kmers.values())
    del kmers
    mean = sum(abundances)/len(abundances)
    ab_sorted = sorted(abundances)
    lower_75 = ab_sorted[1*len(abundances)//4]
    median = ab_sorted[len(abundances)//2]
    upper_75 = ab_sorted[3*len(abundances)//4]
    print("kmers Mean:", mean, "median", median, "lower_75", lower_75, "upper_75", upper_75)

    # s = indexing.minstrobes_iter(genome, k_size, order = 2, N_1 = 50 )
    minstrobes2 = defaultdict(int)
    for acc, seq in genome.items():
        for i,s in enumerate(indexing.minstrobes_iter(seq, k_size, order = 2, N_1 = 50 )):
            minstrobes2[hash(s)] += 1
        # print(acc)

    abundances = list(minstrobes2.values())
    del minstrobes2
    mean = sum(abundances)/len(abundances)
    ab_sorted = sorted(abundances)
    lower_75 = ab_sorted[1*len(abundances)//4]
    median = ab_sorted[len(abundances)//2]
    upper_75 = ab_sorted[3*len(abundances)//4]
    print("minstrobes2 Mean:", mean, "median", median, "lower_75", lower_75, "upper_75", upper_75)


    # s = indexing.minstrobes(genome, k_size, order = 3, N_1 = 25, N_2 = 25)
    minstrobes3 = defaultdict(int)
    for acc, seq in genome.items():
        for i, s in enumerate(indexing.minstrobes_iter(seq, k_size, order = 3, N_1 = 25, N_2 = 25 )):
            minstrobes3[hash(s)] += 1
    abundances = list(minstrobes3.values())
    del minstrobes3
    mean = sum(abundances)/len(abundances)
    ab_sorted = sorted(abundances)
    lower_75 = ab_sorted[1*len(abundances)//4]
    median = ab_sorted[len(abundances)//2]
    upper_75 = ab_sorted[3*len(abundances)//4]
    print("minstrobes3 Mean:", mean, "median", median, "lower_75", lower_75, "upper_75", upper_75)




    # s = indexing.randstrobes(genome, k_size, order = 2, N_1 = 50 )
    randstrobes2 = defaultdict(int)
    for acc, seq in genome.items():
        for s in indexing.randstrobes_iter(seq, k_size, order = 2, N_1 = 50 ):
            randstrobes2[s] += 1

    abundances = list(randstrobes2.values())
    del randstrobes2
    mean = sum(abundances)/len(abundances)
    ab_sorted = sorted(abundances)
    lower_75 = ab_sorted[1*len(abundances)//4]
    median = ab_sorted[len(abundances)//2]
    upper_75 = ab_sorted[3*len(abundances)//4]
    print("randstrobes2 Mean:", mean, "median", median, "lower_75", lower_75, "upper_75", upper_75)


    # s = indexing.randstrobes(genome, k_size, order = 3, N_1 = 25, N_2 = 25)
    randstrobes3 = defaultdict(int)
    for acc, seq in genome.items():
        for s in indexing.minstrobes_iter(seq, k_size, order = 3, N_1 = 25, N_2 = 25 ):
            randstrobes3[s] += 1

    abundances = list(randstrobes3.values())
    del randstrobes3
    mean = sum(abundances)/len(abundances)
    ab_sorted = sorted(abundances)
    lower_75 = ab_sorted[1*len(abundances)//4]
    median = ab_sorted[len(abundances)//2]
    upper_75 = ab_sorted[3*len(abundances)//4]
    print("randstrobes3 Mean:", mean, "median", median, "lower_75", lower_75, "upper_75", upper_75)





if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calc identity", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--fasta', type=str,  default=False, help='Path to consensus fastq file(s)')
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