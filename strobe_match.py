#! /usr/bin/env python

from __future__ import print_function
import os,sys
import argparse

# import errno
# from time import time
# import re

# import random
# import parasail
# import pysam

from collections import defaultdict, deque
from sys import stdout
from array import array
from itertools import zip_longest

from modules import help_functions




import operator
def argmin(values):
    min_index, min_value = min(enumerate(values), key=operator.itemgetter(1))
    return min_index, min_value

def rc(string):
    rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'t', 'c':'g', 'g':'c', 't':'a', 'N':'N', 'X':'X', 'n':'n', 'Y':'R', 'R':'Y', 'K':'M', 'M':'K', 'S':'S', 'W':'W', 'B':'V', 'V':'B', 'H':'D', 'D':'H', 'y':'r', 'r':'y', 'k':'m', 'm':'k', 's':'s', 'w':'w', 'b':'v', 'v':'b', 'h':'d', 'd':'h'}
    rev_comp = ''.join([rev_nuc[nucl] for nucl in reversed(string)])
    return(rev_comp)


def randstrobe_order2(hash_seq_list, start, stop, hash_m1, k_size, prime):
    min_index, min_value = argmin([ (hash_m1+ hash_seq_list[i]) % prime for i in range(start, stop)])
    min_hash_val = hash_m1 + hash_seq_list[start + min_index]
    return min_index, min_hash_val

def seq_to_strobes_iter(seq, k_size, w_min, w_max, prime):
    hash_seq_list = [hash(seq[i:i+k_size]) for i in range(len(seq) - k_size +1)]
    for p in range(len(seq) - 2*k_size + 1):
        hash_m1 = hash_seq_list[p]
        window_p_start = p + k_size + w_min if p + w_max <= len(hash_seq_list) else max( (p + k_size + w_min) -  (p+k_size+w_max - len(hash_seq_list)), p+ k_size )
        window_p_end = min(p + w_max, len(hash_seq_list))
        min_index, hash_value = randstrobe_order2(hash_seq_list, window_p_start, window_p_end, hash_m1, k_size, prime)
        p2 = window_p_start + min_index
        yield p, p2, hash_value


def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


def build_strobemer_index(refs, k_size, w, prime):
    idx = defaultdict(lambda :array("L"))
    ref_id_to_accession = {}
    cntr = 0
    for r_id, (ref_acc, (seq, _)) in enumerate(help_functions.readfq(refs)):
        ref_id_to_accession[r_id] = ref_acc
        for p1, p2, hash_val in seq_to_strobes_iter(seq, k_size, 0, w, prime):
            idx[hash_val].append(r_id)
            idx[hash_val].append(p1)
            idx[hash_val].append(p2)
            cntr += 1
            if cntr % 1000000 == 0:
                print("{0} strobemers created from references".format(cntr))
        # print(hash_val, r_id, pos)
    return idx, ref_id_to_accession, cntr


def sort_merge(sorted_list):
    sort_merged_list = []
    curr_merge = sorted_list[0]
    for i, t1 in enumerate( sorted_list[:-1] ):
        r_id, r_pos, q_pos, length = t1
        r2_id, r2_pos, q2_pos, length2 = sorted_list[i+1]
        # print(i, r_id, r_pos, r2_id, r2_pos)
        print(r2_pos, q2_pos)
        if r_id == r2_id:  
            # print("OK", q2_pos <= q_pos + length <= q2_pos+ length, r2_pos <= r_pos + length <= r2_pos + length)
            print("2", q2_pos, q_pos + length, q2_pos+ length, r2_pos, r_pos + length, r2_pos + length)
            # overlapping on both query and ref
            if q2_pos <= q_pos + length <= q2_pos+ length  and r2_pos <= r_pos + length <= r2_pos + length:
                curr_merge = (r_id, curr_merge[1], curr_merge[2], curr_merge[3] + (q2_pos + length - (q_pos + length) )  )
                # print("HERER")

            else:
                # time to add old element
                sort_merged_list.append(curr_merge)
                curr_merge = sorted_list[i+1]

        else:
            # time to add old element
            sort_merged_list.append(curr_merge)
            curr_merge = sorted_list[i+1]
        # print(curr_merge)
        # print(sort_merged_list)
    # print(curr_merge)
    sort_merged_list.append(curr_merge)
    return sort_merged_list

def get_matches(strobes, idx, k, dont_merge_matches,  ref_id_to_accession, acc):
    """
        The merging of matches is a simple linear merging. If there are repetitive matches across e.g. a chromosome
        the merging will be broken up at the repetitive kmer. To solve the merging exactly, we would need
        to solve the collinear chaining problem after we have out matches. There is no such functionality here.

        Another way to solve this is to do a post merging after sorting the merged matches.
        If two merged matches also overlaps, they can be merged again.
    """
    if dont_merge_matches:
        matches = []
        for q_p1, q_p2, h in strobes:
            if h in idx:
                for r_id, r_p1, r_p2 in grouper(idx[h], 3):
                    matches.append( (r_id, r_p1, q_p1, k) )
        return matches
    else:
        cpm = {} # current potential merges
        merged_matches = []
        for q_p1, q_p2, h in strobes:
            if h in idx:
                for r_id, r_p1, r_p2 in grouper(idx[h], 3):
                    # if ref_id_to_accession[r_id] == acc:
                    #     continue
                    if r_id in cpm:
                        # # there is overlap in both reference and query to previous hit and of exact same distance
                        # if q_p1 < cpm[r_id][1] and r_p1 < cpm[r_id][3] and (q_p1 - cpm[r_id][0]) == (r_p1 - cpm[r_id][2]):

                        # there is overlap in both reference and query to previous hit
                        if q_p1 < cpm[r_id][1] and cpm[r_id][0] < q_p2  and cpm[r_id][2] <= r_p1 <= cpm[r_id][3]:
                            cpm[r_id][1] = max(cpm[r_id][1], q_p2 + k)
                            cpm[r_id][3] = max(cpm[r_id][3], r_p2 + k)
                            # if cpm[r_id][1] > q_p2 + k:
                            #     print("LOOL")
                            # if cpm[r_id][3] > r_p2 + k:
                            #     print("LOOL222")                                
                        else: # no overlap in at least one sequence. Output previous match region and add beginning of new match
                            prev_q_p1, prev_q_p2, prev_ref_p1, prev_ref_p2 = cpm[r_id]
                            # assert  prev_q_p2 - prev_q_p1 == prev_ref_p2 - prev_ref_p1
                            merged_matches.append( (r_id, prev_ref_p1, prev_q_p1, prev_q_p2 - prev_q_p1) )
                            cpm[r_id] = [q_p1, q_p2 + k, r_p1, r_p2 + k ]
                    else:
                        cpm[r_id] = [q_p1, q_p2 + k, r_p1, r_p2 + k ]

        # close all open merge intervals
        for r_id, (q_p1, q_pos_stop, r_pos, r_pos_stop) in cpm.items():
            merged_matches.append( (r_id, r_pos, q_p1, q_pos_stop - q_p1) )

        # print(acc, merged_matches)
        # return sorted(merged_matches, key = lambda x: x[2]) 

        # If there are repetitive matches across e.g. a chromosome
        # the merging will be broken up at the repetitive kmer.
        # here we post merge such spuriously broken up overlapping matches

        # sort first by reference id then by sum of reference and query position to resolve perfect repeats!
        new_sort = sorted(merged_matches, key = lambda x: (x[0], x[1]+x[2], x[1] ) )
        merged_matches = sort_merge(new_sort)

        # sort first by reference id then by reference position
        new_sort = sorted(merged_matches, key = lambda x: (x[0], x[1] ) )
        merged_matches = sort_merge(new_sort)

        # sort first by reference id then by query position
        new_sort = sorted(merged_matches, key = lambda x: (x[0], x[2] ) )
        merged_matches = sort_merge(new_sort)

        return sorted(merged_matches, key = lambda x: (x[0], x[2], x[1]) )


        

def seq_to_kmer_iter(seq, k_size):
    hash_seq_list = [hash(seq[i:i+k_size]) for i in range(len(seq) - k_size +1)]
    for p in range(len(seq) - k_size + 1):
        yield p, p, hash_seq_list[p]

def build_kmer_index(refs, k_size):
    idx = defaultdict(lambda :array("L"))
    ref_id_to_accession = {}
    cntr = 0
    for r_id, (ref_acc, (seq, _)) in enumerate(help_functions.readfq(refs)):
        ref_id_to_accession[r_id] = ref_acc
        for p1, p2, hash_val in seq_to_kmer_iter(seq, k_size):
            idx[hash_val].append(r_id)
            idx[hash_val].append(p1)
            idx[hash_val].append(p2)
            cntr += 1
            if cntr % 1000000 == 0:
                print("{0} kmers created from references".format(cntr))
        # print(hash_val, r_id, pos)
    return idx, ref_id_to_accession, cntr


def print_matches_to_file(query_matches, ref_id_to_accession, outfile):
    for q_acc, read_matches in query_matches:
        outfile.write(">{0}\n".format(q_acc))
        for (r_id, ref_p, q_pos, k) in read_matches:
                ref_acc = ref_id_to_accession[r_id]
                outfile.write("{0}\t{1}\t{2}\t{3}\n".format(ref_acc, ref_p, q_pos, k))

def main(args):
    PRIME = 97

    if args.kmer_index:
        idx, ref_id_to_accession, cntr = build_kmer_index(open(args.references,'r'), args.k)
        print("{0} kmers created from references".format(cntr))
    else:
        idx, ref_id_to_accession, cntr = build_strobemer_index(open(args.references,'r'),args.k, args.strobe_w_size, PRIME)
        print("{0} strobemers created from references".format(cntr))


    outfile = open(args.outfile, 'w')
    query_matches = []

    if args.rev_comp:
        outfile_rc = open(args.outfile + "_revcomp", 'w')
        matches_rc = []

    for i, (acc, (seq, _)) in enumerate(help_functions.readfq(open(args.queries, 'r'))):
        if args.kmer_index:
            strobes = [(p1, p2, h) for p1, p2, h in seq_to_kmer_iter(seq, args.k)] 
        else:    
            strobes = [(p1, p2, h) for p1, p2, h in seq_to_strobes_iter(seq, args.k, 0, args.strobe_w_size, PRIME)]

        read_matches = get_matches(strobes, idx, args.k, args.dont_merge_matches, ref_id_to_accession, acc)
        query_matches.append( (acc, read_matches) )

        if i % 1000 == 0:
            print("Finished processing {0} query sequences.".format(i))
            print_matches_to_file(query_matches, ref_id_to_accession, outfile)
            query_matches = []

        if args.rev_comp:
            if args.kmer_index:
                strobes_rc = [(p1, p2, h) for p1, p2, h in seq_to_kmer_iter(rc(seq), args.k)] 
            else:    
                strobes_rc = [(p1, p2, h) for p1, p2, h in seq_to_strobes_iter(rc(seq), args.k, 0, args.strobe_w_size, PRIME)]
            read_matches_rc = get_matches(strobes_rc, idx, args.k, args.dont_merge_matches,ref_id_to_accession, acc)
            matches_rc.append((acc, read_matches_rc))
            if i % 1000 == 0:
                print_matches_to_file(matches_rc, ref_id_to_accession, outfile_rc)
                matches_rc = []

        # sys.exit()
    print_matches_to_file(query_matches, ref_id_to_accession, outfile)
    
    if args.rev_comp:
        print_matches_to_file(matches_rc, ref_id_to_accession, outfile_rc)

    outfile.close()



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calc identity", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--queries', type=str,  default=False, help='Path to query fasta or fastq file')
    parser.add_argument('--references', type=str,  default=False, help='Path to reference fasta or fastq file')
    parser.add_argument('--k', type=int, default=15, help='Strobe size')
    parser.add_argument('--strobe_w_size', type=int, default=50, help='Strobemer window sizes')
    parser.add_argument('--w', type=int, default=1, help='Thinning window size (default = 1, i.e., no thinning)')
    parser.add_argument('--n', type=int, default=2, help='Order on strobes')
    parser.add_argument('--dont_merge_matches', action="store_true",  help='Do not merge matches with this option. It is seriously advised to\
                                                                     merge matches as the files can become huge otherwise and fill up all diskspace.\
                                                                     Do not specify this option unless you know what you are doing! Moslty here for\
                                                                     development/bugchecking purposas. The default option is to merge matches if they\
                                                                     are consectutive on both query and reference to create MAM-like matches \
                                                                     (maximal approximate matches) of various lengths, much like the output of MUMmer. This is\
                                                                     disk space frendilier, although these files can get large too.')
    parser.add_argument('--outfile', type=str,  default=None, help='TSV match file.')
    parser.add_argument('--compress', type=str,  default=None, help='Compress output')
    parser.add_argument('--rev_comp', action="store_true",  help='Match reverse complement of reads (output to separate file)')
    parser.add_argument('--kmer_index', action="store_true",  help='Kmers can be used instead of strobemers, used for performance comparison')
    # parser.add_argument('--pickled_subreads', type=str, help='Path to an already parsed subreads file in pickle format')
    # parser.set_defaults(which='main')
    args = parser.parse_args()



    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()

    main(args)

