#! /usr/bin/env python

from __future__ import print_function
import os,sys
import argparse

import errno
from time import time
import itertools
import shutil
import subprocess
import math
import re

import random
import parasail
import pysam
# For eventual De Bruijn graph approach
import itertools
from collections import defaultdict, deque
from sys import stdout

from scipy.stats import poisson


'''
    Below code taken from https://github.com/lh3/readfq/blob/master/readfq.py
'''
def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].replace(" ", "_"), [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, (''.join(seqs), None) # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, (seq, ''.join(seqs)); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, (seq, None) # yield a fasta record instead
                break



def reverse_complement(string):
    rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'t', 'c':'g', 'g':'c', 't':'a', 'N':'N', 'X':'X', 'n':'n', 'Y':'R', 'R':'Y', 'K':'M', 'M':'K', 'S':'S', 'W':'W', 'B':'V', 'V':'B', 'H':'D', 'D':'H', 'y':'r', 'r':'y', 'k':'m', 'm':'k', 's':'s', 'w':'w', 'b':'v', 'v':'b', 'h':'d', 'd':'h'}
    rev_comp = ''.join([rev_nuc[nucl] for nucl in reversed(string)])
    return(rev_comp)

def minimizers_comb_iterator(minimizers, k, x_low, x_high):
    # print("read")
    for i, (m1, p1) in enumerate(minimizers[:-1]):
        m1_curr_spans = []
        for j, (m2, p2) in enumerate(minimizers[i+1:]):
            if x_low < p2 - p1 and p2 - p1 <= x_high:
                m1_curr_spans.append( (m2, p2) )
                # yield (m1,p1), (m2, p2) 
            elif p2 - p1 > x_high:
                break
        yield (m1, p1), m1_curr_spans[::-1]


def get_minimizers_and_positions(reads, w, k, hash_fcn = 'lex'):
    # 1. homopolymenr compress read and obtain minimizers
    M = {}
    for r_id in reads:
        (acc, seq, qual) = reads[r_id]
        if hash_fcn == "lex":
            minimizers = get_kmer_minimizers(seq, k, w)
        # elif hash_fcn == "rev_lex":
        #     minimizers = get_kmer_maximizers(seq, k, w)

        M[r_id] = minimizers

    return M

def get_kmer_minimizers(seq, k_size, w_size):
    # kmers = [seq[i:i+k_size] for i in range(len(seq)-k_size) ]
    w = w_size - k_size
    window_kmers = deque([seq[i:i+k_size] for i in range(w +1)])
    curr_min = min(window_kmers)
    minimizers = [ (curr_min, list(window_kmers).index(curr_min)) ]

    for i in range(w+1,len(seq) - k_size):
        new_kmer = seq[i:i+k_size]
        # updateing window
        discarded_kmer = window_kmers.popleft()
        window_kmers.append(new_kmer)

        # we have discarded previous windows minimizer, look for new minimizer brute force
        if curr_min == discarded_kmer: 
            curr_min = min(window_kmers)
            minimizers.append( (curr_min, list(window_kmers).index(curr_min) + i - w ) )

        # Previous minimizer still in window, we only need to compare with the recently added kmer 
        elif new_kmer < curr_min:
            curr_min = new_kmer
            minimizers.append( (curr_min, i) )

    return minimizers

def kmer_counter(reads):
    k_size = 7
    count = defaultdict(int)
    count_pos = defaultdict(list)
    reads_kmers = {}
    for r_i in reads:
        seq = reads[r_i]
        # seq_hpol_comp = ''.join(ch for ch, _ in itertools.groupby(seq))
        read_kmers = deque([seq[i:i+k_size] for i in range(len(seq) - k_size +1)])
        kmer_prev = ""
        for p, kmer in enumerate(read_kmers):

            count[kmer] += 1
            count_pos[kmer].append(p)
            # if kmer_prev != kmer:
            #     count_pos[kmer].append(p)
            #     kmer_prev = kmer               
            # else:
            #     count_pos[kmer][-1] = p

        reads_kmers[r_i] = read_kmers

    cnt_sorted = sorted(count.items(), key = lambda x: x[1], reverse=True)
    # print(len(cnt_sorted),cnt_sorted)
    # for acc in reads_kmers:
    #     print([count[k] for k in reads_kmers[acc]])
        # sys.exit()
    # if "ATCAAGG" in count:
    #     print("ATCAAGG",count["ATCAAGG"])
    #     sys.exit()

    # for k in list(count.keys()):
    #     if count[k] > 3:
    #         print(sum(sorted(count_pos[k]))/float(len(count_pos[k])), count_pos[k], k)
    #     else:
    #         del count[k]
    #         del count_pos[k]

    # sys.exit()
    return count, count_pos


def h(n, f , k):
    return f(n) % k 

import operator

def argmin(values):
    min_index, min_value = min(enumerate(values), key=operator.itemgetter(1))
    return min_index, min_value


def get_spaced_kmer_order2(subseq, m_size, h_table, h_table_inv):
    k1 = subseq[0:m_size]
    f = lambda x: x
    mod = 2**26
    # h_val1 = h(h_table_inv[k1], f , mod) 
    min_index, min_value = argmin([ h(h_table_inv[k1] - h_table_inv[subseq[i:i+m_size]], f, mod) for i in range(m_size, len(subseq) - m_size + 1)])
    # print(min_index)
    min_k2 = subseq[m_size+ min_index:m_size+ min_index+m_size]
    # print(len(k1 + min_k2))
    return k1 + min_k2


def get_spaced_kmer_order3(subseq, m_size, h_table, h_table_inv, N_1, N_2):
    # print(len(subseq),m_size, N_1, N_2, m_size, m_size+ N_1 - m_size + 1, [i for i in range(m_size + N_1, m_size + N_1 + N_2 - m_size + 1)])
    k1 = subseq[0:m_size]
    f = lambda x: x
    mod = 2**26
    min_index, min_value = argmin([ h(h_table_inv[k1] - h_table_inv[subseq[i:i+m_size]], f, mod) for i in range(m_size, m_size+ N_1 - m_size + 1)])
    min_k2 = subseq[m_size + min_index: m_size+ min_index+m_size]

    min_index, min_value = argmin([ h(h_table_inv[k1] - h_table_inv[min_k2] + h_table_inv[subseq[i:i+m_size]], f, mod) for i in range(m_size + N_1, m_size + N_1 + N_2 - m_size + 1)])
    min_k3 = subseq[m_size + N_1 + min_index: m_size + N_1+ min_index+m_size]

    return k1 + min_k2 + min_k3


def number_matching_kmers(seq1, seq2, k_size):
    kmers_seq1 = set([seq1[i:i+k_size] for i in range(len(seq1) - k_size +1)])
    kmers_seq2 = set([seq2[i:i+k_size] for i in range(len(seq2) - k_size +1)])
    n = len(kmers_seq1 & kmers_seq2)
    print("kmers length set1:", len(kmers_seq1))
    print("kmers length set2:", len(kmers_seq2))
    print("common kmers:", n)

    assert k_size % 2 == 0, "Not even kmer length, results will be different"
    N_1 = 50 
    m_size = k_size//2
    h_table_inv = { seq1[i:i+m_size] : random.getrandbits(32) for i in range(len(seq1) - m_size +1)}
    h_table_inv2 = { seq2[i:i+m_size] : random.getrandbits(32)for i in range(len(seq2) - m_size +1)}
    # print(h_table_inv)
    # h_table_inv = { seq1[i:i+m_size] : seq1[i:i+m_size] for i in range(len(seq1) - m_size +1)}
    # h_table_inv2 = { seq2[i:i+m_size] : seq2[i:i+m_size] for i in range(len(seq2) - m_size +1)}

    h_table_inv.update(h_table_inv2)
    h_table = {v: k for k,v in h_table_inv.items() }
    sp2_kmers_seq1 = set([get_spaced_kmer_order2(seq1[i:min(i+N_1, len(seq1))], m_size, h_table, h_table_inv) for i in range(len(seq1) - k_size +1)])
    sp2_kmers_seq2 = set([get_spaced_kmer_order2(seq2[i:min(i+N_1, len(seq2))], m_size, h_table, h_table_inv) for i in range(len(seq2) - k_size +1)])
    n2 = len(sp2_kmers_seq1 & sp2_kmers_seq2)
    print("length set1:", len(sp2_kmers_seq1))
    print("length set2:", len(sp2_kmers_seq2))
    print("Common 2-spaced kmers:", n2)
    print()
    print()

    N_1 = 25 
    N_2 = 25     
    m_size = k_size//3
    h_table_inv = { seq1[i:i+m_size] : random.getrandbits(32) for i in range(len(seq1) - m_size +1)}
    h_table_inv2 = { seq2[i:i+m_size] : random.getrandbits(32)for i in range(len(seq2) - m_size +1)}
    h_table_inv.update(h_table_inv2)
    h_table = {v: k for k,v in h_table_inv.items() }

    sp3_kmers_seq1 = set([get_spaced_kmer_order3(seq1[i:min(i+m_size+N_1+N_2, len(seq1))], m_size, h_table, h_table_inv, min(N_1 + N_2, len(seq1)-i - m_size)//2, min(N_1 + N_2, len(seq1)-i - m_size)//2) for i in range(len(seq1) - k_size +1)])
    sp3_kmers_seq2 = set([get_spaced_kmer_order3(seq2[i:min(i+m_size+N_1+N_2, len(seq2))], m_size, h_table, h_table_inv, min(N_1 + N_2, len(seq1)-i - m_size)//2, min(N_1 + N_2, len(seq1)-i - m_size)//2) for i in range(len(seq2) - k_size +1)])
    n2 = len(sp3_kmers_seq1 & sp3_kmers_seq2)
    print("length set1:", len(sp3_kmers_seq1))
    print("length set2:", len(sp3_kmers_seq2))
    print("Common 3-spaced kmers:", n2)
    print()
    print()

    return n, n2

def get_intervals(kmers_pos1, matches):
    if not matches:
        return []
    tmp = [p for (p, k) in kmers_pos1 if k in matches ]
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
    kmers_pos1 = [(p, seq1[i:i+k_size]) for p, i in enumerate(range(len(seq1) - k_size +1))]
    kmers_seq1 = set([seq1[i:i+k_size] for i in range(len(seq1) - k_size +1)])
    kmers_seq2 = set([seq2[i:i+k_size] for i in range(len(seq2) - k_size +1)])

    matches  = kmers_seq1 & kmers_seq2
    print("length seq1:", len(seq1))
    print("length seq2:", len(seq2))
    ivls = get_intervals(kmers_pos1, matches)

    print("kmer intervals:", ivls)

    assert k_size % 2 == 0, "Not even kmer length, results will be different"
    N_1 = 50 
    m_size = k_size//2
    h_table_inv = { seq1[i:i+m_size] : random.getrandbits(32) for i in range(len(seq1) - m_size +1)}
    h_table_inv2 = { seq2[i:i+m_size] : random.getrandbits(32)for i in range(len(seq2) - m_size +1)}
    # h_table_inv = { seq1[i:i+m_size] : seq1[i:i+m_size] for i in range(len(seq1) - m_size +1)}
    # h_table_inv2 = { seq2[i:i+m_size] : seq2[i:i+m_size] for i in range(len(seq2) - m_size +1)}

    h_table_inv.update(h_table_inv2)
    h_table = {v: k for k,v in h_table_inv.items() }
    sp2_kmers_pos1 = [ (p, get_spaced_kmer_order2(seq1[i:min(i+N_1, len(seq1))], m_size, h_table, h_table_inv)) for p, i in enumerate(range(len(seq1) - k_size +1))]

    sp2_kmers_seq1 = set([get_spaced_kmer_order2(seq1[i:min(i+N_1, len(seq1))], m_size, h_table, h_table_inv) for i in range(len(seq1) - k_size +1)])
    sp2_kmers_seq2 = set([get_spaced_kmer_order2(seq2[i:min(i+N_1, len(seq2))], m_size, h_table, h_table_inv) for i in range(len(seq2) - k_size +1)])
    matches2 = sp2_kmers_seq1 & sp2_kmers_seq2
    print("length seq1:", len(seq1))
    print("length seq2:", len(seq2))
    ivls = get_intervals(sp2_kmers_pos1, matches2)
    print("2-spaced kmers intervals:", ivls)

    N_1 = 25 
    N_2 = 25     
    m_size = k_size//3
    h_table_inv = { seq1[i:i+m_size] : random.getrandbits(32) for i in range(len(seq1) - m_size +1)}
    h_table_inv2 = { seq2[i:i+m_size] : random.getrandbits(32)for i in range(len(seq2) - m_size +1)}
    h_table_inv.update(h_table_inv2)
    h_table = {v: k for k,v in h_table_inv.items() }
    sp3_kmers_pos1 = [ (p, get_spaced_kmer_order3(seq1[i:min(i+m_size+N_1+N_2, len(seq1))], m_size, h_table, h_table_inv,min(N_1 + N_2, len(seq1)-i - m_size)//2, min(N_1 + N_2, len(seq1)-i - m_size)//2)) for p, i in enumerate(range(len(seq1) - k_size +1))]

    sp3_kmers_seq1 = set([get_spaced_kmer_order3(seq1[i:min(i+m_size+N_1+N_2, len(seq1))], m_size, h_table, h_table_inv,min(N_1 + N_2, len(seq1)-i - m_size)//2, min(N_1 + N_2, len(seq1)-i - m_size)//2) for i in range(len(seq1) - k_size +1)])
    sp3_kmers_seq2 = set([get_spaced_kmer_order3(seq2[i:min(i+m_size+N_1+N_2, len(seq2))], m_size, h_table, h_table_inv,min(N_1 + N_2, len(seq1)-i - m_size)//2, min(N_1 + N_2, len(seq1)-i - m_size)//2) for i in range(len(seq2) - k_size +1)])
    matches3 = sp3_kmers_seq1 & sp3_kmers_seq2
    print("length seq1:", len(seq1))
    print("length seq2:", len(seq2))
    ivls = get_intervals(sp3_kmers_pos1, matches3)
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
    seq2 = "".join([seq1[i] if i not in muts else random.choice(reverse_complement(seq1[i])) for i in range(len(seq1))])
    n,n2 = number_matching_kmers(seq1, seq2, k)
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

