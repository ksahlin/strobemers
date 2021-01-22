import os,sys
import argparse

import random

from modules import indexing, help_functions

def get_intervals(mers, matches):
    if not matches:
        return []
    tmp = [p for (p, k) in mers.items() if k in matches ]
    ivls = []
    iv_start = tmp[0]
    p_prev = tmp[0]
    for p in tmp[1:]:
        if p == p_prev + 1:
            p_prev += 1
            continue
        else:
            ivls.append((iv_start, p_prev))
            p_prev = p
            iv_start = p
    ivls.append((iv_start, p_prev))
    return ivls

def positions_matching_kmers(seq1, seq2, k_size):

    print("length seq1:", len(seq1))
    print("length seq2:", len(seq2))

    #kmers
    kmers_pos1 = {p : seq1[i:i+k_size] for p, i in enumerate(range(len(seq1) - k_size +1))}
    kmers_seq1 = set([seq1[i:i+k_size] for i in range(len(seq1) - k_size +1)])
    kmers_seq2 = set([seq2[i:i+k_size] for i in range(len(seq2) - k_size +1)])
    matches  = kmers_seq1 & kmers_seq2
    ivls = get_intervals(kmers_pos1, matches)
    print("kmers nr_hits:", len(matches))    
    print("kmer intervals:", ivls)


    assert k_size % 2 == 0, "Not even kmer length, results will be different"
    randomers1 = indexing.randomers(seq1, k_size, order = 2, N_1 = 50 )
    randomers2 = indexing.randomers(seq2, k_size, order = 2, N_1 = 50 )    
    matches2 = set(randomers1.values()) & set(randomers2.values())
    ivls = get_intervals(randomers1, matches2)
    print("2-spaced kmers nr_hits:", len(matches2))    
    print("2-spaced kmers intervals:", ivls)


    assert k_size % 3 == 0, "Not div by 3 kmer length, results will be different"
    randomers1 = indexing.randomers(seq1, k_size, order = 3, N_1 = 25, N_2 = 25)
    randomers2 = indexing.randomers(seq2, k_size, order = 3, N_1 = 25, N_2 = 25 )
    matches3 = set(randomers1.values()) & set(randomers2.values())
    ivls = get_intervals(randomers1, matches3)
    print("3-spaced kmers nr_hits:", len(matches3))    
    print("3-spaced kmers intervals:", ivls)

def main(args):
    # reads = { acc: seq for i, (acc, (seq, _)) in enumerate(readfq(open(args.fasta, 'r')))}
    # seq1, seq2 = list(reads.values())
    L = 1000
    k = 30
    mut_freq = 0.07
    seq1 = "".join([random.choice("ACGT") for i in range(L)])
    muts = set(random.sample(range(len(seq1)),int(L*mut_freq)))
    muts = set(range(20,1000,20)) #set([20,40,60,80])
    # print(muts)
    seq2 = "".join([seq1[i] if i not in muts else random.choice(help_functions.reverse_complement(seq1[i])) for i in range(len(seq1))])
    positions_matching_kmers(seq1, seq2, k)
    # strobe_kmers = get_strobe_kmers(reads)
    # kmers = get_kmers(reads)

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