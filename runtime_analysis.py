import os,sys
import argparse

import random
from collections import defaultdict

from modules import indexing2, help_functions

from time import time

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


def time_datastructure(seq, k_size, w_size, data_structure):
    w = 1
    w_low = 1
    w_high = w_size+1
    prime = 997
    all_mers = defaultdict(int)


    if data_structure == "kmer":
        datastructure = "kmers"
        for p, hash_value in indexing2.kmer_iter(seq, k_size, w):
            all_mers[hash_value] += 1


    elif data_structure == "minstrobes2":
        for p, p2, hash_value in indexing2.seq_to_minstrobes2_iter(seq, k_size, w_low, w_high, prime, w):
            all_mers[hash_value] += 1

    elif  data_structure == "minstrobes3":
        for p, p2, p3, hash_value in indexing2.seq_to_minstrobes3_iter(seq, k_size, w_low, w_high, prime, w):
            all_mers[hash_value] += 1

    elif data_structure == "randstrobes2":
        for p, p2, hash_value in indexing2.seq_to_randstrobes2_iter(seq, k_size, w_low, w_high, prime, w):
            all_mers[hash_value] += 1

    elif data_structure == "randstrobes3":
        for p, p2, p3, hash_value  in indexing2.seq_to_randstrobes3_iter(seq, k_size, w_low, w_high, prime, w):
            all_mers[hash_value] += 1

    elif data_structure == "hybridstrobes2":
        if w_size > 4: # not defined below number of active queues which are 4 in this implementation
            for p, p2, hash_value in indexing2.seq_to_hybridstrobes2_iter(seq, k_size, w_low, w_high):
                all_mers[hash_value] += 1

def main(args):

    all_strings = [ "".join([random.choice("ACGT") for i in range(args.l)]) for j in range(args.n) ]

    timings = {"kmer" : 0,
                "minstrobes2": 0,
                "minstrobes3": 0,
                "randstrobes2": 0,
                "randstrobes3": 0,
                "hybridstrobes2": 0}


    for k_size in [18,36,54,60,72]: #[18,24,30]:
        for w_size in [1,10,20,30,40,50, 100]: #[18,24,30]:
            # print(w_size)
            start = time()
            for s in all_strings:
                time_datastructure(s, k_size, w_size, "kmer")
            elapsed = time() - start
            timings["kmer"] = elapsed
            # print("kmers", k_size, w_size,  elapsed)

            start = time()
            for s in all_strings:
                time_datastructure(s, k_size, w_size, "minstrobes2")
            elapsed = time() - start
            timings["minstrobes2"] = elapsed
            # print("minstrobes2", k_size, w_size,  elapsed)

            start = time()
            for s in all_strings:
                time_datastructure(s, k_size, w_size, "minstrobes3")
            elapsed = time() - start
            timings["minstrobes3"] = elapsed
            # print("minstrobes3", k_size, w_size,  elapsed)

            start = time()
            for s in all_strings:
                time_datastructure(s, k_size, w_size, "randstrobes2")
            elapsed = time() - start
            timings["randstrobes2"] = elapsed
            # print("randstrobes2", k_size, w_size,  elapsed)

            start = time()
            for s in all_strings:
                time_datastructure(s, k_size, w_size, "randstrobes3")
            elapsed = time() - start
            timings["randstrobes3"] = elapsed
            # print("randstrobes3", k_size, w_size,  elapsed)

            start = time()
            for s in all_strings:
                time_datastructure(s, k_size, w_size, "hybridstrobes2")
            elapsed = time() - start
            timings["hybridstrobes2"] = elapsed
            # print("minstrobes2", k_size, w_size,  elapsed)

            ms2 = round(timings["minstrobes2"]  / timings["kmer"], 1)
            ms3 = round(timings["minstrobes3"]  / timings["kmer"], 1)
            km = round(timings["kmer"]  / timings["kmer"], 1)
            rs2 = round(timings["randstrobes2"]  / timings["kmer"], 1)
            rs3 = round(timings["randstrobes3"]  / timings["kmer"], 1)
            hs2 = round(timings["hybridstrobes2"]  / timings["kmer"], 1)

            # for ds, v in timings.items():
            #     # print(k,v, timings["kmer"])
            print("{0} & {1} & {2} & {3} & {4} & {5} & {6} & {7}".format(k_size, w_size,  km, ms2, ms3, rs2, rs3, hs2) )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calc identity", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # parser.add_argument('--fasta', type=str,  default=False, help='Path to consensus fastq file(s)')
    # parser.add_argument('--kmers',  action="store_true", help='Kmer size')
    # parser.add_argument('--minstrobes2', action="store_true", help='Kmer size')
    # parser.add_argument('--minstrobes3',  action="store_true", help='Kmer size')
    # parser.add_argument('--randstrobes2',  action="store_true", help='Kmer size')
    # parser.add_argument('--randstrobes3',  action="store_true", help='Kmer size')
    # parser.add_argument('--spaced_dense',  action="store_true", help='Kmer size')
    # parser.add_argument('--spaced_sparse',  action="store_true", help='Kmer size')
    parser.add_argument('--n', type=int, default=1000, help='Number of strings (experiments) to run')
    parser.add_argument('--l', type=int, default=10000, help='Length of strings')
    # parser.add_argument('--w', type=int, default=20, help='Window size')
    # parser.add_argument('--outfolder', type=str,  default=None, help='A fasta file with transcripts that are shared between samples and have perfect illumina support.')
    # parser.add_argument('--pickled_subreads', type=str, help='Path to an already parsed subreads file in pickle format')
    # parser.set_defaults(which='main')
    args = parser.parse_args()



    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()

    main(args)