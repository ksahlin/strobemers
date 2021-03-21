#! /usr/bin/env python

from __future__ import print_function
import os,sys
import argparse

import copy
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

MAX = sys.maxsize

def argmin(values):
    min_index, min_value = min(enumerate(values), key=operator.itemgetter(1))
    return min_index, min_value

def rc(string):
    rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'t', 'c':'g', 'g':'c', 't':'a', 'N':'N', 'X':'X', 'n':'n', 'Y':'R', 'R':'Y', 'K':'M', 'M':'K', 'S':'S', 'W':'W', 'B':'V', 'V':'B', 'H':'D', 'D':'H', 'y':'r', 'r':'y', 'k':'m', 'm':'k', 's':'s', 'w':'w', 'b':'v', 'v':'b', 'h':'d', 'd':'h'}
    rev_comp = ''.join([rev_nuc[nucl] for nucl in reversed(string)])
    return(rev_comp)


def thinner(hash_list, w):
    """
        Input: a list with hash values 
        Output: A list with tuples: (pos in original list, minimim hash value) for each window of w hashes
    """
    window_hashes = deque(hash_list[:w])
    min_index, curr_min_hash = argmin(window_hashes)
    thinned_hash_list = [ (min_index, curr_min_hash) ]

    for i in range(w, len(hash_list) + w-1):
        if i >= len(hash_list):
            new_hash = MAX
        else:
            new_hash = hash_list[i]
        # updating window
        discarded_hash = window_hashes.popleft()
        window_hashes.append(new_hash)

        # we have discarded previous windows minimizer, look for new minimizer brute force
        if curr_min_hash == discarded_hash: 
            min_index, curr_min_hash = argmin(window_hashes)
            thinned_hash_list.append( (min_index + i + 1 - w, curr_min_hash) )

        # Previous minimizer still in window, we only need to compare with the recently added kmer 
        elif new_hash < curr_min_hash:
            curr_min_hash = new_hash
            thinned_hash_list.append( (i, curr_min_hash) )


    return thinned_hash_list


def randstrobe_order3(hash_seq_list, start1, stop1, start2, stop2, hash_m1, k_size, prime):
    min_index1, min_value = argmin([ (hash_m1 - hash_seq_list[i][1]) % prime for i in range(start1, stop1)])
    min_hash_val = hash_m1 - hash_seq_list[start1 + min_index1][1]

    min_index2, min_value = argmin([ (min_hash_val - hash_seq_list[i][1]) % prime for i in range(start2, stop2)])
    min_hash_val = min_hash_val + 2*hash_seq_list[start2 + min_index2][1]

    return min_index1, min_index2, min_hash_val

def seq_to_randstrobes3_iter(seq, k_size, strobe_w_min_offset, strobe_w_max_offset, prime, w):
    hash_seq_list = [(i, hash(seq[i:i+k_size])) for i in range(len(seq) - k_size +1)]
    if w > 1:
        hash_seq_list_thinned = thinner([h for i,h in hash_seq_list], w) # produce a subset of positions, still with samme index as in full sequence
    else:
        hash_seq_list_thinned = hash_seq_list
    
    # assert len(hash_seq_list[:-k_size]) == len(hash_seq_list) - k_size

    for (p, hash_m1) in hash_seq_list_thinned: #[:-k_size]:
        if p >= len(hash_seq_list) - 2*k_size:
            break
        # hash_m1 = hash_seq_list[p]
        window_p2_start = p + k_size + strobe_w_max_offset + strobe_w_min_offset if p + 2*strobe_w_max_offset <= len(hash_seq_list) else max( (p + k_size + strobe_w_max_offset + strobe_w_min_offset) -  (p+k_size+2*strobe_w_max_offset - len(hash_seq_list)), p + 2*k_size )
        window_p2_end = min(p + 2*strobe_w_max_offset, len(hash_seq_list))

        window_p1_start = p + k_size + strobe_w_min_offset if p + 2*strobe_w_max_offset <= len(hash_seq_list) else max(p+ k_size,  len(hash_seq_list)  + 2*(strobe_w_min_offset - strobe_w_max_offset))
        window_p1_end = min(p + strobe_w_max_offset, len(hash_seq_list)- k_size)
        # print(window_p1_start, window_p1_end,  window_p2_start, window_p2_end, len(seq))
        # assert window_p1_start < window_p1_end
        # print(window_p1_start, window_p1_end)
        min_index_s1, min_index_s2, hash_value = randstrobe_order3(hash_seq_list, window_p1_start, window_p1_end, window_p2_start, window_p2_end, hash_m1, k_size, prime)
        p2 = window_p1_start + min_index_s1
        p3 = window_p2_start + min_index_s2
        yield p, p2, p3, hash_value


def randstrobe_order2(hash_seq_list, start, stop, hash_m1, k_size, prime):
    min_index, min_value = argmin([ (hash_m1 - hash_seq_list[i][1]) % prime for i in range(start, stop)])
    min_hash_val = hash_m1 - hash_seq_list[start + min_index][1]
    return min_index, min_hash_val


def seq_to_randstrobes2_iter(seq, k_size, strobe_w_min_offset, strobe_w_max_offset, prime, w):
    hash_seq_list = [(i, hash(seq[i:i+k_size])) for i in range(len(seq) - k_size +1)]
    if w > 1:
        hash_seq_list_thinned = thinner([h for i,h in hash_seq_list], w) # produce a subset of positions, still with samme index as in full sequence
    else:
        hash_seq_list_thinned = hash_seq_list
    
    # assert len(hash_seq_list[:-k_size]) == len(hash_seq_list) - k_size

    for (p, hash_m1) in hash_seq_list_thinned: #[:-k_size]:
        if p >= len(hash_seq_list) - k_size:
            break
        # hash_m1 = hash_seq_list[p]
        window_p_start = p + k_size + strobe_w_min_offset if p + strobe_w_max_offset <= len(hash_seq_list) else max( (p + k_size + strobe_w_min_offset) -  (p+k_size+strobe_w_max_offset - len(hash_seq_list)), p+ k_size )
        window_p_end = min(p + strobe_w_max_offset, len(hash_seq_list))
        # print(window_p_start, window_p_end)
        min_index, hash_value = randstrobe_order2(hash_seq_list, window_p_start, window_p_end, hash_m1, k_size, prime)
        p2 = window_p_start + min_index
        yield p, p2, hash_value


def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


def build_randstrobe3_index(refs, k_size, strobe_w_min_offset, strobe_w_max_offset, prime, w):
    idx = defaultdict(lambda :array("L"))
    ref_id_to_accession = {}
    cntr = 0
    for r_id, (ref_acc, (seq, _)) in enumerate(help_functions.readfq(refs)):
        ref_id_to_accession[r_id] = ref_acc
        for p1, p2, p3, hash_val in seq_to_randstrobes3_iter(seq, k_size, strobe_w_min_offset, strobe_w_max_offset, prime, w):
            idx[hash_val].append(r_id)
            idx[hash_val].append(p1)
            idx[hash_val].append(p2)
            idx[hash_val].append(p3)
            cntr += 1
            if cntr % 1000000 == 0:
                print("{0} randstrobes created from references".format(cntr))
        # print(hash_val, r_id, pos)
    return idx, ref_id_to_accession, cntr


def build_randstrobe2_index(refs, k_size, strobe_w_min_offset, strobe_w_max_offset, prime, w):
    idx = defaultdict(lambda :array("L"))
    ref_id_to_accession = {}
    cntr = 0
    for r_id, (ref_acc, (seq, _)) in enumerate(help_functions.readfq(refs)):
        ref_id_to_accession[r_id] = ref_acc
        for p1, p2, hash_val in seq_to_randstrobes2_iter(seq, k_size, strobe_w_min_offset, strobe_w_max_offset, prime, w):
            idx[hash_val].append(r_id)
            idx[hash_val].append(p1)
            idx[hash_val].append(p2)
            cntr += 1
            if cntr % 1000000 == 0:
                print("{0} randstrobes created from references".format(cntr))
        # print(hash_val, r_id, pos)
    return idx, ref_id_to_accession, cntr



def seq_to_minstrobes2_iter(seq, k_size, strobe_w_min_offset, strobe_w_max_offset, prime, w):
    hash_seq_list = [(i, hash(seq[i:i+k_size])) for i in range(len(seq) - k_size +1)]
    strobes = deque(thinner([h for i,h in hash_seq_list], strobe_w_max_offset - strobe_w_min_offset)) # produce a subset of positions, still with samme index as in full sequence

    if w > 1:
        hash_seq_list_thinned = thinner([h for i,h in hash_seq_list], w) # produce a subset of positions, still with same index as in full sequence
    else:
        hash_seq_list_thinned = hash_seq_list
    
    # assert len(hash_seq_list[:-k_size]) == len(hash_seq_list) - k_size

    for (p, hash_m1) in hash_seq_list_thinned: #[:-k_size]:
        if p >= len(hash_seq_list) - k_size:
            break
        # print(p,len(hash_seq_list) - k_size, len(hash_seq_list), len(seq), len(strobes), strobes)
        if p + k_size + strobe_w_min_offset < len(seq):
            while strobes[0][0] < min(p + k_size + strobe_w_min_offset, len(hash_seq_list)-1):
                l = strobes.popleft()

        # print(p, len(hash_seq_list) - k_size, len(hash_seq_list), len(seq), len(strobes), strobes[0])
        p2, hash_val = strobes[0]
        hash_value = hash_m1 - hash_val
        yield p, p2, hash_value


def build_minstrobe2_index(refs, k_size, strobe_w_min_offset, strobe_w_max_offset, prime, w):
    idx = defaultdict(lambda :array("L"))
    ref_id_to_accession = {}
    cntr = 0
    for r_id, (ref_acc, (seq, _)) in enumerate(help_functions.readfq(refs)):
        ref_id_to_accession[r_id] = ref_acc
        for p1, p2, hash_val in seq_to_minstrobes2_iter(seq, k_size, strobe_w_min_offset, strobe_w_max_offset, prime, w):
            idx[hash_val].append(r_id)
            idx[hash_val].append(p1)
            idx[hash_val].append(p2)
            cntr += 1
            if cntr % 1000000 == 0:
                print("{0} minstrobes created from references".format(cntr))
        # print(hash_val, r_id, pos)
    return idx, ref_id_to_accession, cntr


def seq_to_minstrobes3_iter(seq, k_size, strobe_w_min_offset, strobe_w_max_offset, prime, w):
    hash_seq_list = [(i, hash(seq[i:i+k_size])) for i in range(len(seq) - k_size +1)]
    strobes = deque(thinner([h for i,h in hash_seq_list], strobe_w_max_offset - strobe_w_min_offset)) # produce a subset of positions, still with samme index as in full sequence
    strobes2 = copy.deepcopy(strobes)
    if w > 1:
        hash_seq_list_thinned = thinner([h for i,h in hash_seq_list], w) # produce a subset of positions, still with same index as in full sequence
    else:
        hash_seq_list_thinned = hash_seq_list
    
    # assert len(hash_seq_list[:-k_size]) == len(hash_seq_list) - k_size

    for (p, hash_m1) in hash_seq_list_thinned: #[:-k_size]:
        if p >= len(hash_seq_list) - 2*k_size:
            break
        # print(p,len(hash_seq_list) - k_size, len(hash_seq_list), len(seq), len(strobes), strobes)

        if p + strobe_w_max_offset + strobe_w_min_offset < len(seq):
            while strobes2[0][0] <  min(p + strobe_w_min_offset + strobe_w_max_offset, len(hash_seq_list)-1):
                l = strobes2.popleft()

        if p + k_size + strobe_w_min_offset < len(seq):
            while strobes[0][0] <  min(p + k_size + strobe_w_min_offset, len(hash_seq_list)-1):
                l = strobes.popleft()
        # print(p, len(hash_seq_list) - k_size, len(hash_seq_list), len(seq), len(strobes), strobes[0], len(strobes2))

        p2, hash_val2 = strobes[0]
        p3, hash_val3 = strobes2[0]
        hash_value = hash_m1 - hash_val2 + 2*hash_val3
        yield p, p2, p3, hash_value


def build_minstrobe3_index(refs, k_size, strobe_w_min_offset, strobe_w_max_offset, prime, w):
    idx = defaultdict(lambda :array("L"))
    ref_id_to_accession = {}
    cntr = 0
    for r_id, (ref_acc, (seq, _)) in enumerate(help_functions.readfq(refs)):
        ref_id_to_accession[r_id] = ref_acc
        for p1, p2, p3, hash_val in seq_to_minstrobes3_iter(seq, k_size, strobe_w_min_offset, strobe_w_max_offset, prime, w):
            idx[hash_val].append(r_id)
            idx[hash_val].append(p1)
            idx[hash_val].append(p2)
            idx[hash_val].append(p3)
            cntr += 1
            if cntr % 1000000 == 0:
                print("{0} minstrobes created from references".format(cntr))
        # print(hash_val, r_id, pos)
    return idx, ref_id_to_accession, cntr

# def sort_merge(sorted_list):
#     sort_merged_list = []
#     curr_merge = sorted_list[0]
#     for i, t1 in enumerate( sorted_list[:-1] ):
#         r_id, r_pos, q_pos, length = t1
#         r2_id, r2_pos, q2_pos, length2 = sorted_list[i+1]
#         # print(i, r_id, r_pos, r2_id, r2_pos)
#         # print(r2_pos, q2_pos)
#         if r_id == r2_id:  
#             # print("OK", q2_pos <= q_pos + length <= q2_pos+ length, r2_pos <= r_pos + length <= r2_pos + length)
#             # print("2", q2_pos, q_pos + length, q2_pos+ length, r2_pos, r_pos + length, r2_pos + length)
#             # overlapping on both query and ref
#             # print(q2_pos + length2, q_pos + length, curr_merge[3])
#             if q2_pos <= q_pos + length <= q2_pos+ length  and r2_pos <= r_pos + length <= r2_pos + length:
#                 # curr_merge = (r_id, curr_merge[1], curr_merge[2], max(q2_pos + length2, q_pos + length ) -  q_pos ) # hit length on query sequence
#                 curr_merge = (r_id, curr_merge[1], curr_merge[2], max(r2_pos + length2, r_pos + length ) -  r_pos ) # hit length on reference sequence
#                 # print("HERER")

#             else:
#                 # time to add old element
#                 sort_merged_list.append(curr_merge)
#                 curr_merge = sorted_list[i+1]

#         else:
#             # time to add old element
#             sort_merged_list.append(curr_merge)
#             curr_merge = sorted_list[i+1]
#         # print(curr_merge)
#         # print(sort_merged_list)
#     # print(curr_merge)
#     sort_merged_list.append(curr_merge)
#     return sort_merged_list


def get_matches3(strobes, idx, k, dont_merge_matches,  ref_id_to_accession, acc, selfalign):
    """
        The merging of matches is a simple linear merging. If there are repetitive matches across e.g. a chromosome
        the merging will be broken up at the repetitive kmer. To solve the merging exactly, we would need
        to solve the collinear chaining problem after we have out matches. There is no such functionality here.

        Another way to solve this is to do a post merging after sorting the merged matches.
        If two merged matches also overlaps, they can be merged again.
    """
    if dont_merge_matches:
        matches = []
        for q_p1, q_p2, q_p3, h in strobes:
            # print()
            # print("Q", q_p1)
            if h in idx:
                for r_id, r_p1, r_p2, r_p3 in grouper(idx[h], 4):
                    # print("R", r_id, r_p1)
                    matches.append( (r_id, r_p1, q_p1, r_p3 - r_p1 + k) )
        return sorted(matches, key = lambda x: (x[0], x[2], x[1]) )
    else:
        cpm = {} # current potential merges
        merged_matches = []
        for q_p1, q_p2, q_p3, h in strobes: # iterate over query in ascending order
            if h in idx:
                for r_id, r_p1, r_p2, r_p3 in grouper(idx[h], 4): # iterate over references, all in ascending order
                    # remove self matches with below if statement, for now commented out to find eventual bugs
                    if not selfalign and ref_id_to_accession[r_id] == acc:
                        continue
                    if r_id in cpm:
                        is_added_to_an_interval_query = False
                        # print(q_p1, list(cpm[r_id].keys()))
                        for end_q in list(cpm[r_id].keys()):
                            # print()
                            # print("r_id",r_id, "end_q", end_q)
                            if q_p1 <= end_q: # overlap on query
                                is_added_to_an_interval_query = True  
                                is_added_to_an_interval_reference = False  
                                # print(list(cpm[r_id][end_q].keys()))  
                                for end_r in list(cpm[r_id][end_q].keys()):
                                    # print("Case1 end_r", end_r)
                                    # print(q_p1, )
                                    prev_q_p1, prev_q_p2, prev_ref_p1, prev_ref_p2 = cpm[r_id][end_q][end_r]
                                    # print(r_id,q_p1, "CRUCIAL:",prev_q_p1, prev_q_p2, prev_ref_p1, prev_ref_p2)
                                    # print(r_id, q_p1, cpm[r_id][end_q][end_r])
                                    # check all refs
                                    new_q_p2 = max(prev_q_p2, q_p3 + k)
                                    if prev_ref_p1 <= r_p1 <= end_r: # Overlap on reference
                                        is_added_to_an_interval_reference = True  
                                        # print("OK", prev_ref_p1, r_p1, end_r)
                                        # print("lol", prev_ref_p1, r_p1, end_r)
                                        new_r_p2 = max(end_r, r_p3 + k)
                                        del cpm[r_id][end_q][end_r]
                                        if not cpm[r_id][end_q]:
                                            del cpm[r_id][end_q]
                                        if new_q_p2 not in cpm[r_id]:
                                            cpm[r_id][ new_q_p2 ] = {}
                                            cpm[r_id][ new_q_p2 ][new_r_p2] = ( prev_q_p1, new_q_p2, prev_ref_p1, new_r_p2)
                                            # print("new:", ( prev_q_p1, new_q_p2, prev_ref_p1, new_r_p2) )
                                        elif new_r_p2 not in cpm[r_id][ new_q_p2 ]:
                                            cpm[r_id][ new_q_p2 ][new_r_p2] = ( prev_q_p1, new_q_p2, prev_ref_p1, new_r_p2)
                                            # print("appended:", ( prev_q_p1, new_q_p2, prev_ref_p1, new_r_p2) )
                                        else:
                                            # print("Was already present:", cpm[r_id][ new_q_p2 ][new_r_p2], "attempted new:", ( prev_q_p1, new_q_p2, prev_ref_p1, new_r_p2) )
                                            ( old_q_p1, new_q_p2, old_ref_p1, new_r_p2) = cpm[r_id][ new_q_p2 ][new_r_p2]
                                            cpm[r_id][ new_q_p2 ][new_r_p2] = ( min(old_q_p1, prev_q_p1), new_q_p2, min(old_ref_p1, prev_ref_p1), new_r_p2)

                                        # cpm[r_id][ new_q_p2 ][new_r_p2] = [ prev_q_p1, new_q_p2, prev_ref_p1, new_r_p2]
                                    
                                if not is_added_to_an_interval_reference:
                                    if new_q_p2 not in cpm[r_id]:
                                        cpm[r_id][ new_q_p2 ] = {} 
                                        cpm[r_id][ new_q_p2 ][r_p3 + k] = (q_p1, new_q_p2, r_p1, r_p3 + k)
                                        # print("new added1:", (q_p1, new_q_p2, r_p1, r_p3 + k) )

                                    elif r_p3 + k not in cpm[r_id][new_q_p2]:
                                        cpm[r_id][ new_q_p2 ][r_p3 + k] = (q_p1, new_q_p2, r_p1, r_p3 + k )
                                        # print("new added2:", (q_p1, new_q_p2, r_p1, r_p3 + k ) )
                                    else:
                                        # print("Was already present:", cpm[r_id][ new_q_p2 ][r_p3 + k], "attempted new:", (q_p1, new_q_p2, r_p1, r_p3 + k ) )
                                        ( old_q_p1, new_q_p2, old_ref_p1, new_r_p2) = cpm[r_id][ new_q_p2 ][r_p3 + k]
                                        cpm[r_id][ new_q_p2 ][new_r_p2] = ( min(old_q_p1, q_p1), new_q_p2, min(old_ref_p1, r_p1), new_r_p2)

                            else:
                                # print("Case2 end_r", end_r)
                            # revove the intervals that we have passed on the query here to not make the cpm dict too large...
                            # add to merged_matches dict
                                for r_end in cpm[r_id][end_q]:
                                    (q_pos_start, q_pos_stop, r_pos, r_pos_stop) = cpm[r_id][end_q][r_end]
                                    merged_matches.append( (r_id, r_pos, q_pos_start, r_pos_stop - r_pos) )
                                del cpm[r_id][end_q]
                              
                        if not is_added_to_an_interval_query: # no overlap with prev query sequences
                            cpm[r_id][q_p3 + k] = {}
                            cpm[r_id][q_p3 + k][r_p3 + k] = (q_p1, q_p3 + k, r_p1, r_p3 + k)
                    else:
                        cpm[r_id] = { q_p3 + k : {r_p3 + k : (q_p1, q_p3 + k, r_p1, r_p3 + k) }}
                        # cpm[r_id] = [q_p1, q_p2 + k, r_p1, r_p3 + k ]

        # close all open merge intervals
        for r_id  in cpm.keys():
            for q_stop in cpm[r_id]:
                for r_stop in cpm[r_id][q_stop]:
                    (q_p1, q_pos_stop, r_pos, r_pos_stop) = cpm[r_id][q_stop][r_stop]
                    merged_matches.append( (r_id, r_pos, q_p1, r_pos_stop - r_pos) )
        # print(merged_matches)
        if not merged_matches:
            return []


        return sorted(set(merged_matches), key = lambda x: (x[0], x[2], x[1]) )




def get_matches(strobes, idx, k, dont_merge_matches,  ref_id_to_accession, acc, selfalign):
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
            # print()
            # print("Q", q_p1)
            if h in idx:
                for r_id, r_p1, r_p2 in grouper(idx[h], 3):
                    # print("R", r_id, r_p1)
                    matches.append( (r_id, r_p1, q_p1, r_p2 - r_p1 + k) )
        return sorted(matches, key = lambda x: (x[0], x[2], x[1]) )
    else:
        cpm = {} # current potential merges
        merged_matches = []
        for q_p1, q_p2, h in strobes: # iterate over query in ascending order
            if h in idx:
                # print()
                # print("----------------", q_p1)
                # print(cpm)
                # print("All pos:", idx[h])
                # prev_r_id, prev_hit_r_p1,prev_hit_r_p2 = 0,0,0 # these only keep track of identical consecutive kmers/strobes
                for r_id, r_p1, r_p2 in grouper(idx[h], 3): # iterate over references, all in ascending order
                    # if prev_r_id == r_id and r_p1 == prev_hit_r_p1 + 1 and r_p2 == prev_hit_r_p2+1:
                    #     prev_r_id = r_id
                    #     prev_hit_r_p1 = r_p1
                    #     prev_hit_r_p2 = r_p2
                    #     update_relevant_pos()

                    #     continue

                    # print(r_p1, r_p2)
                    # remove self matches with below if statement, for now commented out to find eventual bugs
                    if not selfalign and ref_id_to_accession[r_id] == acc:
                        continue
                    if r_id in cpm:
                        is_added_to_an_interval_query = False
                        # print(q_p1, list(cpm[r_id].keys()))
                        for end_q in list(cpm[r_id].keys()):
                            # print()
                            # print("r_id",r_id, "end_q", end_q)
                            if q_p1 <= end_q: # overlap on query
                                is_added_to_an_interval_query = True  
                                is_added_to_an_interval_reference = False  
                                # print(list(cpm[r_id][end_q].keys()))  
                                for end_r in list(cpm[r_id][end_q].keys()):
                                    # print("Case1 end_r", end_r)
                                    # print(q_p1, )
                                    prev_q_p1, prev_q_p2, prev_ref_p1, prev_ref_p2 = cpm[r_id][end_q][end_r]
                                    # print(r_id,q_p1, "CRUCIAL:",prev_q_p1, prev_q_p2, prev_ref_p1, prev_ref_p2)
                                    # print(r_id, q_p1, cpm[r_id][end_q][end_r])
                                    # check all refs
                                    new_q_p2 = max(prev_q_p2, q_p2 + k)
                                    if prev_ref_p1 <= r_p1 <= end_r: # Overlap on reference
                                        is_added_to_an_interval_reference = True  
                                        # print("OK", prev_ref_p1, r_p1, end_r)
                                        # print("lol", prev_ref_p1, r_p1, end_r)
                                        new_r_p2 = max(end_r, r_p2 + k)
                                        del cpm[r_id][end_q][end_r]
                                        if not cpm[r_id][end_q]:
                                            del cpm[r_id][end_q]
                                        if new_q_p2 not in cpm[r_id]:
                                            cpm[r_id][ new_q_p2 ] = {}
                                            cpm[r_id][ new_q_p2 ][new_r_p2] = ( prev_q_p1, new_q_p2, prev_ref_p1, new_r_p2)
                                            # print("new:", ( prev_q_p1, new_q_p2, prev_ref_p1, new_r_p2) )
                                        elif new_r_p2 not in cpm[r_id][ new_q_p2 ]:
                                            cpm[r_id][ new_q_p2 ][new_r_p2] = ( prev_q_p1, new_q_p2, prev_ref_p1, new_r_p2)
                                            # print("appended:", ( prev_q_p1, new_q_p2, prev_ref_p1, new_r_p2) )
                                        else:
                                            # print("Was already present:", cpm[r_id][ new_q_p2 ][new_r_p2], "attempted new:", ( prev_q_p1, new_q_p2, prev_ref_p1, new_r_p2) )
                                            ( old_q_p1, new_q_p2, old_ref_p1, new_r_p2) = cpm[r_id][ new_q_p2 ][new_r_p2]
                                            cpm[r_id][ new_q_p2 ][new_r_p2] = ( min(old_q_p1, prev_q_p1), new_q_p2, min(old_ref_p1, prev_ref_p1), new_r_p2)

                                        # cpm[r_id][ new_q_p2 ][new_r_p2] = [ prev_q_p1, new_q_p2, prev_ref_p1, new_r_p2]
                                    
                                if not is_added_to_an_interval_reference:
                                    if new_q_p2 not in cpm[r_id]:
                                        cpm[r_id][ new_q_p2 ] = {} 
                                        cpm[r_id][ new_q_p2 ][r_p2 + k] = (q_p1, new_q_p2, r_p1, r_p2 + k)
                                        # print("new added1:", (q_p1, new_q_p2, r_p1, r_p2 + k) )

                                    elif r_p2 + k not in cpm[r_id][new_q_p2]:
                                        cpm[r_id][ new_q_p2 ][r_p2 + k] = (q_p1, new_q_p2, r_p1, r_p2 + k )
                                        # print("new added2:", (q_p1, new_q_p2, r_p1, r_p2 + k ) )
                                    else:
                                        # print("Was already present:", cpm[r_id][ new_q_p2 ][r_p2 + k], "attempted new:", (q_p1, new_q_p2, r_p1, r_p2 + k ) )
                                        ( old_q_p1, new_q_p2, old_ref_p1, new_r_p2) = cpm[r_id][ new_q_p2 ][r_p2 + k]
                                        cpm[r_id][ new_q_p2 ][new_r_p2] = ( min(old_q_p1, q_p1), new_q_p2, min(old_ref_p1, r_p1), new_r_p2)

                            else:
                                # print("Case2 end_r", end_r)
                            # revove the intervals that we have passed on the query here to not make the cpm dict too large...
                            # add to merged_matches dict
                                for r_end in cpm[r_id][end_q]:
                                    (q_pos_start, q_pos_stop, r_pos, r_pos_stop) = cpm[r_id][end_q][r_end]
                                    merged_matches.append( (r_id, r_pos, q_pos_start, r_pos_stop - r_pos) )
                                del cpm[r_id][end_q]


                                # # print(end_q, cpm[r_id][end_q][1])
                                # # assert end_q == cpm[r_id][end_q][1]
                                # # there is overlap in both reference and query to previous hit
                                # # `q_1 <= q_2 <= q'_1 +k` and `r_1 <= r_2 <= r'_2+k`
                                # if cpm[r_id][end_q][0] < q_p1 and q_p1 < cpm[r_id][end_q][1] and cpm[r_id][end_q][2] <= r_p1 <= cpm[r_id][end_q][3]:
                                #     prev_q_p1, prev_q_p2, prev_ref_p1, prev_ref_p2 = cpm[r_id][end_q]
                                #     new_q_p2 = max(cpm[r_id][end_q][1], q_p2 + k)
                                #     new_r_p2 = max(cpm[r_id][end_q][3], r_p2 + k)
                                #     del cpm[r_id][end_q]
                                #     cpm[r_id][ new_q_p2 ] = [ prev_q_p1, new_q_p2, prev_ref_p1, new_r_p2]
                                #     is_added_to_an_interval_query = True
                                #     # break
                                # else:
                                #     prev_q_p1, prev_q_p2, prev_ref_p1, prev_ref_p2 = cpm[r_id][end_q]
                                #     new_q_p2 = max(cpm[r_id][end_q][1], q_p2 + k)
                                #     new_r_p2 = max(cpm[r_id][end_q][3], r_p2 + k)
                                #     del cpm[r_id][end_q]
                                #     cpm[r_id][ new_q_p2 ] = [ prev_q_p1, new_q_p2, prev_ref_p1, new_r_p2]
                                #     is_added_to_an_interval_query = True

                                #     # cpm[r_id][end_q][1] = max(cpm[r_id][end_q][1], q_p2 + k)
                                #     # cpm[r_id][end_q][3] = max(cpm[r_id][end_q][3], r_p2 + k)
                                #     # print(cpm[r_id][0], q_p2 + k)
                                #     # print(cpm[r_id][2], r_p2 + k)
                                #     # if cpm[r_id][1] > q_p2 + k:
                                #     #     print("LOOL")
                                #     # if cpm[r_id][3] > r_p2 + k:
                                #     #     print("LOOL222")                                
                        if not is_added_to_an_interval_query: # no overlap with prev query sequences
                            # prev_q_p1, prev_q_p2, prev_ref_p1, prev_ref_p2 = cpm[r_id][end_q]
                            # assert  prev_q_p2 - prev_q_p1 == prev_ref_p2 - prev_ref_p1
                            # print(prev_q_p1,prev_q_p2, prev_q_p2 - prev_q_p1)
                            # print(prev_ref_p1,prev_ref_p2, prev_ref_p2 - prev_ref_p1)
                            # merged_matches.append( (r_id, prev_ref_p1, prev_q_p1, prev_ref_p2 - prev_ref_p1) )
                            # cpm[r_id] = [q_p1, q_p2 + k, r_p1, r_p2 + k ]
                            # print("HERE")
                            cpm[r_id][q_p2 + k] = {}
                            cpm[r_id][q_p2 + k][r_p2 + k] = (q_p1, q_p2 + k, r_p1, r_p2 + k)
                    else:
                        cpm[r_id] = { q_p2 + k : {r_p2 + k : (q_p1, q_p2 + k, r_p1, r_p2 + k) }}
                        # cpm[r_id] = [q_p1, q_p2 + k, r_p1, r_p2 + k ]

        # close all open merge intervals
        for r_id  in cpm.keys():
            for q_stop in cpm[r_id]:
                for r_stop in cpm[r_id][q_stop]:
                    (q_p1, q_pos_stop, r_pos, r_pos_stop) = cpm[r_id][q_stop][r_stop]
                    merged_matches.append( (r_id, r_pos, q_p1, r_pos_stop - r_pos) )
        # print(merged_matches)
        if not merged_matches:
            return []
        # print(acc, merged_matches)
        # return sorted(merged_matches, key = lambda x: x[2]) 

        # # If there are repetitive matches across e.g. a chromosome
        # # the merging will be broken up at the repetitive kmer.
        # # here we post merge such spuriously broken up overlapping matches

        # # sort first by reference id then by sum of reference and query position to resolve perfect repeats!
        # new_sort = sorted(merged_matches, key = lambda x: (x[0], x[1]+x[2], x[1] ) )
        # merged_matches = sort_merge(new_sort)
        # # print(merged_matches)

        # # sort first by reference id then by reference position
        # new_sort = sorted(merged_matches, key = lambda x: (x[0], x[1] ) )
        # merged_matches = sort_merge(new_sort)
        # # print(merged_matches)

        # # sort first by reference id then by query position
        # new_sort = sorted(merged_matches, key = lambda x: (x[0], x[2] ) )
        # merged_matches = sort_merge(new_sort)
        # # print(merged_matches)

        return sorted(set(merged_matches), key = lambda x: (x[0], x[2], x[1]) )


        

def seq_to_kmer_iter(seq, k_size, w):
    hash_seq_list = [(i, hash(seq[i:i+k_size])) for i in range(len(seq) - k_size +1)]
    if w > 1:
        hash_seq_list = thinner([h for i,h in hash_seq_list], w)
    # assert range(len(seq) - k_size + 1) == range(len(hash_seq_list))
    for (p, hash_val) in hash_seq_list:
        yield p, p, hash_val

def build_kmer_index(refs, k_size, w):
    idx = defaultdict(lambda :array("L"))
    ref_id_to_accession = {}
    cntr = 0
    for r_id, (ref_acc, (seq, _)) in enumerate(help_functions.readfq(refs)):
        ref_id_to_accession[r_id] = ref_acc
        for p1, p2, hash_val in seq_to_kmer_iter(seq, k_size, w):
            idx[hash_val].append(r_id)
            idx[hash_val].append(p1)
            idx[hash_val].append(p2)
            cntr += 1
            if cntr % 1000000 == 0:
                print("{0} kmers created from references".format(cntr))
        # print(hash_val, r_id, pos)
    return idx, ref_id_to_accession, cntr


def print_matches_to_file(query_matches, ref_id_to_accession, outfile, reverse):
    for q_acc, read_matches in query_matches:
        if reverse:
            outfile.write("> {0} Reverse\n".format(q_acc))
        else:
            outfile.write("> {0}\n".format(q_acc))
        for (r_id, ref_p, q_pos, k) in read_matches:
                ref_acc = ref_id_to_accession[r_id]
                outfile.write("  {0} {1} {2} {3}\n".format(ref_acc, ref_p, q_pos, k))

def main(args):
    PRIME = 997
    w = args.w

    if args.kmer_index:
        idx, ref_id_to_accession, cntr = build_kmer_index(open(args.references,'r'), args.k, w)
        print("{0} kmers created from references".format(cntr))
        # print(idx)
    elif args.minstrobe_index:
        if args.n == 2:
            idx, ref_id_to_accession, cntr = build_minstrobe2_index(open(args.references,'r'),args.k, args.strobe_w_min_offset, args.strobe_w_max_offset, PRIME, w)
            print("{0} minstrobes created from references".format(cntr))
        elif args.n == 3:
            idx, ref_id_to_accession, cntr = build_minstrobe3_index(open(args.references,'r'),args.k, args.strobe_w_min_offset, args.strobe_w_max_offset, PRIME, w)
            print("{0} minstrobes created from references".format(cntr))            
    else:
        if args.n == 2:
            idx, ref_id_to_accession, cntr = build_randstrobe2_index(open(args.references,'r'),args.k, args.strobe_w_min_offset, args.strobe_w_max_offset, PRIME, w)
            print("{0} randstrobes created from references".format(cntr))
        elif args.n == 3:
            idx, ref_id_to_accession, cntr = build_randstrobe3_index(open(args.references,'r'),args.k, args.strobe_w_min_offset, args.strobe_w_max_offset, PRIME, w)
            print("{0} randstrobes created from references".format(cntr))            


    outfile = open(os.path.join(args.outfolder, args.prefix + ".tsv"), 'w')
    query_matches = []

    if args.rev_comp:
        # outfile_rc = open(os.path.join(args.outfolder, args.prefix + "_revcomp.tsv"), 'w')
        matches_rc = []

    for i, (acc, (seq, _)) in enumerate(help_functions.readfq(open(args.queries, 'r'))):
        if args.kmer_index:
            strobes = [(p1, p2, h) for p1, p2, h in seq_to_kmer_iter(seq, args.k, w)] 
            read_matches = get_matches(strobes, idx, args.k, args.dont_merge_matches, ref_id_to_accession, acc, args.selfalign)

        elif args.minstrobe_index:
            if args.n == 2: 
                strobes = [(p1, p2, h) for p1, p2, h in seq_to_minstrobes2_iter(seq, args.k, args.strobe_w_min_offset, args.strobe_w_max_offset, PRIME, w)]
                read_matches = get_matches(strobes, idx, args.k, args.dont_merge_matches, ref_id_to_accession, acc, args.selfalign)

            elif args.n == 3: 
                strobes = [(p1, p2, p3, h) for p1, p2, p3, h in seq_to_minstrobes3_iter(seq, args.k, args.strobe_w_min_offset, args.strobe_w_max_offset, PRIME, w)]
                read_matches = get_matches3(strobes, idx, args.k, args.dont_merge_matches, ref_id_to_accession, acc, args.selfalign)


        else:
            if args.n == 2: 
                strobes = [(p1, p2, h) for p1, p2, h in seq_to_randstrobes2_iter(seq, args.k, args.strobe_w_min_offset, args.strobe_w_max_offset, PRIME, w)]
                read_matches = get_matches(strobes, idx, args.k, args.dont_merge_matches, ref_id_to_accession, acc, args.selfalign)

            elif args.n == 3: 
                strobes = [(p1, p2, p3, h) for p1, p2, p3, h in seq_to_randstrobes3_iter(seq, args.k, args.strobe_w_min_offset, args.strobe_w_max_offset, PRIME, w)]
                read_matches = get_matches3(strobes, idx, args.k, args.dont_merge_matches, ref_id_to_accession, acc, args.selfalign)

            else:
                raise NotImplementedError

        # print(strobes)
        query_matches.append( (acc, read_matches) )

        if i % 1000 == 0:
            print("Finished processing {0} query sequences.".format(i))
            print_matches_to_file(query_matches, ref_id_to_accession, outfile, False)
            query_matches = []

        if args.rev_comp:
            if args.kmer_index:
                strobes_rc = [(p1, p2, h) for p1, p2, h in seq_to_kmer_iter(rc(seq), args.k, w)] 
                read_matches_rc = get_matches(strobes_rc, idx, args.k, args.dont_merge_matches,ref_id_to_accession, acc, args.selfalign)

            elif args.minstrobe_index:
                if args.n == 2: 
                    strobes = [(p1, p2, h) for p1, p2, h in seq_to_minstrobes2_iter(rc(seq), args.k, args.strobe_w_min_offset, args.strobe_w_max_offset, PRIME, w)]
                    read_matches = get_matches(strobes, idx, args.k, args.dont_merge_matches, ref_id_to_accession, acc, args.selfalign)

                elif args.n == 3: 
                    strobes = [(p1, p2, p3, h) for p1, p2, p3, h in seq_to_minstrobes3_iter(rc(seq), args.k, args.strobe_w_min_offset, args.strobe_w_max_offset, PRIME, w)]
                    read_matches = get_matches3(strobes, idx, args.k, args.dont_merge_matches, ref_id_to_accession, acc, args.selfalign)

                else:
                    raise NotImplementedError

            else:
                if args.n == 2: 
                    strobes_rc = [(p1, p2, h) for p1, p2, h in seq_to_randstrobes2_iter(rc(seq), args.k, args.strobe_w_min_offset, args.strobe_w_max_offset, PRIME, w)]
                    read_matches_rc = get_matches(strobes_rc, idx, args.k, args.dont_merge_matches,ref_id_to_accession, acc, args.selfalign)

                elif args.n == 3: 
                    strobes_rc = [(p1, p2, p3, h) for p1, p2, p3, h in seq_to_randstrobes3_iter(rc(seq), args.k, args.strobe_w_min_offset, args.strobe_w_max_offset, PRIME, w)]
                    read_matches_rc = get_matches3(strobes_rc, idx, args.k, args.dont_merge_matches,ref_id_to_accession, acc, args.selfalign)

                else:
                    raise NotImplementedError

            matches_rc.append((acc, read_matches_rc))
            if i % 1000 == 0:
                print_matches_to_file(matches_rc, ref_id_to_accession, outfile, True)
                matches_rc = []

        # sys.exit()
    print_matches_to_file(query_matches, ref_id_to_accession, outfile, False)
    
    if args.rev_comp:
        print_matches_to_file(matches_rc, ref_id_to_accession, outfile, True)
        # outfile_rc.close()
    outfile.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calc identity", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--queries', type=str,  default=False, help='Path to query fasta or fastq file')
    parser.add_argument('--references', type=str,  default=False, help='Path to reference fasta or fastq file')
    parser.add_argument('--k', type=int, default=15, help='Strobe size')
    parser.add_argument('--strobe_w_min_offset', type=int, default=20, help='Strobemer window start offset from first k-mer. If kmer start at pos i, first\
                                                                            window will start at i+k+strobe_w_min_offset. Default: 20nt donwstream from end of first kmer.')
    parser.add_argument('--strobe_w_max_offset', type=int, default=70, help='Strobemer window end. If kmer start at pos i, first\
                                                                            window will stop at i+k+strobe_w_max_offset. Default: 70nt donwstream from end of first kmer.')
    parser.add_argument('--w', type=int, default=1, help='Thinning window size (default = 1, i.e., no thinning)')
    parser.add_argument('--n', type=int, default=2, help='Order on strobes')
    parser.add_argument('--dont_merge_matches', action="store_true",  help='Do not merge matches with this option. It is seriously advised to\
                                                                     merge matches as the files can become huge otherwise and fill up all diskspace.\
                                                                     Do not specify this option unless you know what you are doing! Moslty here for\
                                                                     development/bugchecking purposas. The default option is to merge matches if they\
                                                                     are consectutive on both query and reference to create MAM-like matches \
                                                                     (maximal approximate matches) of various lengths, much like the output of MUMmer. This is\
                                                                     disk space frendilier, although these files can get large too.')
    parser.add_argument('--outfolder', type=str,  default=None, help='Folder to output TSV match file.')
    parser.add_argument('--prefix', type=str,  default="matches", help='Filename prefix (default "matches").')
    parser.add_argument('--kmer_index', action="store_true",  help='Kmers can be used instead of randstrobes (default), used for performance comparison')
    parser.add_argument('--minstrobe_index', action="store_true",  help='Kmers can be used instead of randstrobes (default), used for performance comparison')
    parser.add_argument('--selfalign', action="store_true",  help='Aligns sequences to itself (mainly used for bugfixing). Default is not align\
                                                                    sequneces to themselves if the same file is given as references and queries.')
    # parser.add_argument('--compress', type=str,  default=None, help='Compress output')
    parser.add_argument('--rev_comp', action="store_true",  help='Match reverse complement of reads (output to separate file)')
    # parser.add_argument('--pickled_subreads', type=str, help='Path to an already parsed subreads file in pickle format')
    # parser.set_defaults(which='main')
    args = parser.parse_args()



    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()

    if args.outfolder and not os.path.exists(args.outfolder):
        os.makedirs(args.outfolder)

    # if args.w != 1:
    #     raise NotImplementedError("Currently only w=1 is allowed, i.e., no thinning is implemented")

    main(args)

