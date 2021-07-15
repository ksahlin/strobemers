#!/usr/bin/env python

import argparse
import sys, os
import random
import pysam
import re

from collections import defaultdict
from collections import namedtuple

import RMQT as RMaxQST

NAM = namedtuple('NAM', ['x', 'y', 'c', 'd', 'val', 'j', "chr_id"])


'''
    Below fast[a/q] reader function taken from https://github.com/lh3/readfq/blob/master/readfq.py
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
        name, seqs, last = last[1:].split()[0], [], None
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

def get_NAM_records(NAM_path, reads):
    '''
        Reads contains all the relevant reads in the batch to read mems from 
    '''
    for i, line in enumerate(open(NAM_path, 'r')):
        if line[0] == '>':
            acc = line[1:].strip()
            if i == 0:
                read_acc = acc  
            else:

                for chr_id in list(read_mems_tmp.keys()):
                    coordinate_sorted_tuples = sorted(read_mems_tmp[chr_id], key = lambda x: x[1])
                    sorted_mems = [ NAM(x,y,c,d,val,j,e_id) for j, (x, y, c, d, val, e_id) in enumerate(coordinate_sorted_tuples) ]
                    read_mems_tmp[chr_id] = sorted_mems

                yield read_acc, read_mems_tmp
                read_acc = acc 
            
            read_mems_tmp = defaultdict(list)

        else:
                vals =  line.split() #11404_11606           1     11405       202
                chr_id = vals[0]
                mem_len = int(vals[3])
                # convert to 0-indexed reference as in python
                # however, for NAM length last coordinate is inclusive of the hit in NAM solvers, not as in python end-indexing
                mem_read_start = int(vals[2]) - 1
                mem_genome_start = int(vals[1]) - 1
                                
                info_tuple = ( mem_genome_start, mem_genome_start + mem_len - 1,
                                mem_read_start, mem_read_start + mem_len - 1, 
                                mem_len, chr_id)
                read_mems_tmp[chr_id].append( info_tuple )


    for chr_id in list(read_mems_tmp.keys()):
        coordinate_sorted_tuples = sorted(read_mems_tmp[chr_id], key = lambda x: x[1])
        sorted_mems = [ NAM(x,y,c,d,val,j,e_id) for j, (x, y, c, d, val, e_id) in enumerate(coordinate_sorted_tuples) ]
        read_mems_tmp[chr_id] = sorted_mems
    # print("READ {0} LINES.".format(i))
    yield read_acc, read_mems_tmp



def argmax(iterable):
    return max(enumerate(iterable), key=lambda x: x[1])[0]

def max_both(iterable):
    return max(enumerate(iterable), key=lambda x: x[1])

def all_solutions_c_max_indicies(C, C_max):
    return [i for i, c in enumerate(C) if c == C_max] 

def reconstruct_all_solutions(mems, all_C_max_indicies, trace_vector, C, mam_mode = False):
    # solution_index = argmax(C)
    solutions = []
    for solution_index in all_C_max_indicies:
        value = C[solution_index]
        solution = []
        while solution_index > 0:
            if mam_mode:
                solution.append(mems[solution_index])  # trace vector is shifted on to the right so need to remove 1 from vectore to get j_index 
            else:
                solution.append(mems[solution_index - 1])  # trace vector is shifted on to the right so need to remove 1 from vectore to get j_index 
            solution_index = trace_vector[solution_index]
        solutions.append( solution[::-1] )
    return value, solutions

def colinear_solver_read_coverage(mems, max_intron):
    """
        Algorithm 15.1 in Genome scale algorithmic design, Makinen et al.

        Using lists instead of Binary search trees for range max queries,
        so n^2 time complexity instead of O(n log n).

        each NAM is an Namedtuple. python object

    """
    # assert mems == sorted(mems, key = lambda x: x.y )

    if len(mems) > 1000:
        print('NAM',len(mems))

    # print("Going in to NAM chaining:", mems)
    T = [ (v.d, v.val)  for v in mems]
    I = [ (v.d, v.val)  for v in mems]
    
    C = [0]*(len(T)+1)
    traceback_vector = [None]*(len(T)+1)

    for j in range(len(T)):
        v =  mems[j]

        # linear scan -- replace with range max Q tree
        T_values = [(j_prime, c_val) for j_prime, c_val in enumerate(C[1:]) if  mems[j_prime].d < v.c and j_prime < j and v.y - mems[j_prime].y < max_intron ]
        if T_values:
            # print(T_values)
            T_traceback_index, max_c_value_case_a = max(reversed(T_values), key=lambda x: x[1])
        else:
            max_c_value_case_a = 0
            T_traceback_index = -1

        I_values = [(j_prime, c_val) for j_prime, c_val in enumerate(C[1:]) if v.c <= mems[j_prime].d  <= v.d and j_prime < j and v.y - mems[j_prime].y < max_intron ]
        if I_values:
            # print(I_values)
            I_values_plus_chord_diff = [ (j_prime, c_val + (v.d - mems[j_prime].d)) for j_prime, c_val in I_values]
            I_traceback_index, max_c_value_case_b = max(reversed(I_values_plus_chord_diff), key=lambda x: x[1])
            # I_v_prev_coord = mems[I_traceback_index].d
            # C_b[j] = (v.d - I_v_prev_coord) + max_c_value_case_b # shouldnt it be v.d - v_tmp.d
            C_b = max_c_value_case_b # shouldnt it be v.d - v_tmp.d

        else:
            I_v_prev_coord = v.c - 1
            I_traceback_index = -1
            max_c_value_case_b = 0
            C_b = 0


        C_a = (v.d - v.c + 1) +  max_c_value_case_a

        index, value = max_both([C_a, C_b])
        C[j+1] = value
        if index == 0: # Updating with C_a
            j_prime = T_traceback_index
        else: # Updating with C_b
            j_prime = I_traceback_index

        if j_prime < 0: # first j (i.e. j=0) 
            traceback_vector[j+1]= 0
        else:
            traceback_vector[j+1]= j_prime + 1

    solution_index = argmax(C)
    # print(C)
    # print(traceback_vector)
    C_max = C[solution_index]
    all_C_max_indicies = all_solutions_c_max_indicies(C, C_max)
    # print(all_C_max_indicies)
    # print("number solutions with the same score:", all_solutions_c_max_indicies(C, C_max))
    C_max, solutions = reconstruct_all_solutions(mems, all_C_max_indicies, traceback_vector, C)
    # solutions = []
    # print("NAM Solution:", solution[::-1])
    return solutions, C_max #, is_unique_solution(C)
    # traceback(C, best_solution_index)

def make_leafs_power_of_2(mems):
    nodes = []
    nodes.append( RMaxQST.Node(-1, -1, -2**32, -1) ) # add an start node in case a search is smaller than any d coord in tree
    for i, mem in enumerate(mems):
        # if i > 3: break
        m = RMaxQST.Node(mem.d, mem.j, -2**32, mem.j)
        nodes.append(m)

    for i in range(20):
        if len(nodes) == 2**i or len(nodes) == 2**(i+1):
            break
        elif 2**i < len(nodes) < 2**(i+1):
            remainder = 2**(i+1) - len(nodes) 
            for i in range(remainder):
                nodes.append( RMaxQST.Node(-1, -i - 2, -2**32, -i - 2) ) # fill up nodes to have leaves a power of 2
            break

    leafs = sorted(nodes, key= lambda x: x.d)
    # n = len(leafs)
    return leafs


def n_logn_read_coverage(mems):
    """
        Algorithm 15.1 in Genome scale algorithmic design, Makinen et al.

        Using Binary search trees for range max queries,
        so n log n time complexity. Each mem is an Namedtuple. python object

    """
    # assert mems == sorted(mems, key=lambda x: x.y)

    # if len(mems) > 1000:
    #     print('NAM',len(mems))

    T_leafs = make_leafs_power_of_2(mems)
    I_leafs = make_leafs_power_of_2(mems)
    n = len(T_leafs)
    T = [0 for i in range(2 * n) ]  
    I = [0 for i in range(2 * n) ]  
    # T_leafs = copy.deepcopy(leafs)
    RMaxQST.construct_tree(T, T_leafs, n)
    # I_leafs = copy.deepcopy(leafs)
    RMaxQST.construct_tree(I, I_leafs, n)

    mem_to_leaf_index = {l.j : i for i,l in enumerate(T_leafs)}

    C = [0]* (len(mems) + 1) #(len(leafs))
    trace_vector = [None]*(len(mems) + 1)

    RMaxQST.update(T, 0, 0, n) # point update 
    RMaxQST.update(I, 0, 0, n) # point update 

    for j, mem in enumerate(mems):
        leaf_to_update = mem_to_leaf_index[j]
        c = mem.c
        T_max, j_prime_a, node_pos  = RMaxQST.range_query(T, -1, c-1, len(T_leafs)) 
        # print("C_a:",  T_max +  mem.d - mem.c + 1, j_prime_a, node_pos, leaf_to_update )
        # print("T TREE:", [(s, zz.j, zz.d, zz.Cj, zz.j_max) for s, zz in enumerate(T) if type(zz) != int])
        C_a =  T_max +  mem.d - mem.c + 1  # add the mem_length to T since disjoint

        if T_max < 0:
            print("BUG", T_max)
            sys.exit()

        
        d = mem.d
        I_max, j_prime_b, node_pos  = RMaxQST.range_query(I, c, d, len(I_leafs))         
        # print("C_b:", I_max +  mem.d, I_max, j_prime_b, node_pos, leaf_to_update )
        # print( I_max, mem.d, mems[j_prime_b].d, mems[j_prime_b])
        # print("I TREE:", [(s, zz.j, zz.d, zz.Cj, zz.j_max) for s, zz in enumerate(I) if type(zz) != int])
        C_b =  I_max +  mem.d #- mems[j_prime_b].d   # add the part of the mem that is not overlapping

        # if C_b < 0:
        #     print("BUG")
        #     sys.exit()

        index, value = max_both([C_a, C_b])
        C[j+1] = value
        if index == 0: # Updating with C_a
            j_prime = j_prime_a
        else: # Updating with C_b
            j_prime = j_prime_b


        if j_prime < 0: # any of the additional leaf nodes (with negative index) we add to make number of leafs 2^n
            trace_vector[j+1] = 0
        elif value == 0: # first j (i.e. j=0) 
            trace_vector[j+1]= 0
        else:
            trace_vector[j+1] = j_prime +1

        RMaxQST.update(T, leaf_to_update, value, n) # point update 
        RMaxQST.update(I, leaf_to_update, value - mem.d, n) # point update 


    # C_max, solution = reconstruct_solution(mems, C, trace_vector)
    # print("C", C)
    # print(trace_vector)

    solution_index = argmax(C)
    C_max = C[solution_index]
    all_C_max_indicies = all_solutions_c_max_indicies(C, C_max)
    # print("number solutions with the same score:", all_solutions_c_max_indicies(C, C_max))
    C_max, solutions = reconstruct_all_solutions(mems, all_C_max_indicies, trace_vector, C)

    return solutions, C_max #, is_unique_solution(C)


def read_sam(sam_file):
    SAM_file = pysam.AlignmentFile(sam_file, "r", check_sq=False)
    read_positions = {} # acc -> [ref_id, ref_start, refstop]


    for read in SAM_file.fetch(until_eof=True):
        if read.flag == 0:
            read_positions[read.query_name] = (read.reference_name, read.reference_start, read.reference_end, False)
        elif read.flag == 16:
            # print(read.query_name, len(read_positions))
            read_positions[read.query_name] = (read.reference_name, read.reference_start, read.reference_end, True)
        elif read.flag == 4:
            read_positions[read.query_name] = False

    return read_positions  


def compute_overlap_with_truth_genome(true_chr_id, t_start, t_stop, nams):
    tot_overlap = 0
    tot_nam_length = 0
    for n in nams:
        # print(n.chr_id)
        # print(true_chr_id, n.chr_id, true_chr_id == n.chr_id, t_start, t_stop, n.x, n.y )
        if n.chr_id == true_chr_id:
            if t_start <= n.x <=  n.y <= t_stop:
                tot_overlap += n.y - n.x
            elif t_stop < n.x:
                pass
            elif n.y < t_start:
                pass
            elif t_start <= n.x <= t_stop < n.y:
                tot_overlap += t_stop - n.x
            elif  n.x < t_start <= n.y <= t_stop:
                tot_overlap += n.y - t_start
            elif n.x < t_start <  t_stop < n.y :
                tot_overlap += t_stop - t_start
            else:
                print("BUG!", true_chr_id, t_start, t_stop, n.x, n.y )
                 

        tot_nam_length += n.y - n.x
    # print(tot_overlap, tot_nam_length)
    # print(tot_overlap/tot_nam_length)
    return tot_overlap/tot_nam_length


def read_paf(paf_file):
    read_positions = {} # acc -> [ref_id, ref_start, refstop]
    for line in open(paf_file, 'r'):
        vals = line.split()
        read_acc, ref_name, reference_start, reference_end  = vals[0], vals[5], int(vals[7]), int(vals[8])
        if read_acc not in read_positions:
            read_positions[read_acc] = {}

        read_positions[read_acc][ref_name] = (reference_start, reference_end)
    return read_positions

def overlap(q_a, q_b, p_a, p_b):
    assert q_a <= q_b and p_a <= p_b
    # if (q_a == q_b) or (p_a == p_b):
    #     print("Cigar bug")
    return  (p_a <= q_a <= p_b) or (p_a <= q_b <= p_b) or (q_a <= p_a <= q_b) or (q_a <= p_b <= q_b)

def compute_overlap_with_truth_read(true_locations, ref_id, read_id, nams, read_lengths):
    # get overlap span coordinates on 'reference read' from genome alignment coordinates
    ref_chr_id, ref_t_start, ref_t_end, ref_is_rc =  true_locations[ref_id]
    query_chr_id,  query_t_start,  query_t_end, query_is_rc =  true_locations[read_id]
    start_offset, end_offset = 0,0
    if ref_chr_id == query_chr_id:
        if overlap(ref_t_start,ref_t_end,query_t_start,query_t_end):
            # t_start = 0
            # t_stop = read_lengths[ref_id]
            # Reference is fwd
            if (not ref_is_rc):
                if query_t_start < ref_t_start:
                    t_start = 0
                else:
                    start_offset = query_t_start - ref_t_start
                    t_start = start_offset
                if ref_t_end < query_t_end:
                    t_stop = read_lengths[ref_id]
                else:
                    end_offset = ref_t_end - query_t_end
                    t_stop = read_lengths[ref_id] - end_offset

            # Reference is rev comp
            else:
                if query_t_start < ref_t_start:
                    t_stop = read_lengths[ref_id]
                else:
                    end_offset = query_t_start - ref_t_start
                    t_stop = read_lengths[ref_id] - end_offset
                if ref_t_end < query_t_end:
                    t_start = 0
                else:
                    start_offset = ref_t_end - query_t_end
                    t_start = start_offset

        else:
            return 0.0
    else:
        return 0.0

    # print(nams)
    tot_overlap = 0
    tot_nam_length = 0
    for n in nams:
        if t_start <= n.x <= n.y <= t_stop:
            tot_overlap += n.y - n.x
        elif t_stop < n.x:
            pass
        elif n.y < t_start:
            pass
        elif t_start <= n.x <= t_stop < n.y:
            tot_overlap += t_stop - n.x
        elif  n.x < t_start <= n.y <= t_stop:
            tot_overlap += n.y - t_start
        elif n.x < t_start <  t_stop < n.y :
            tot_overlap += t_stop - t_start
        else:
            print("BUG!", true_chr_id, t_start, t_stop, n.x, n.y )
                 

        tot_nam_length += n.y - n.x
    # print(tot_overlap, tot_nam_length)
    # print(tot_overlap/tot_nam_length)
    return tot_overlap/tot_nam_length


def e_size(collinear_chain_nam_sizes, genome_size):
    sum_of_squares = sum([x**2 for x in collinear_chain_nam_sizes])
    # print("sum_of_squares", sum_of_squares, genome_size)
    return sum_of_squares/genome_size

def get_time_and_mem(runtime_file):
    time, mem = '', ''
    for l in runtime_file:
        if re.search('[\d.]+ real', l):
            time = re.search('[\d.]+ real', l)
            time = float(time.group().split(" ")[0])
        if re.search('[\d.]+  maximum resident set size', l):
            mem =  re.search('[\d.]+  maximum resident set size', l)
            mem = float(mem.group().split(" ")[0])
    # print(time, mem)
    return time, mem/1000000

def main(args):

    time_in_sec, mem_in_mb = get_time_and_mem(open(args.runtime_file, 'r'))
    if args.refs:
        ref_lengths = { acc : len(seq) for acc, (seq,qual) in readfq(open(args.refs, 'r'))}
        # print(ref_lengths)
    read_coverage_solution = {}
    total_disjoint_matches = 0
    tot_genome_length = 0

    for (query_acc, nams) in get_NAM_records(args.matchfile, ref_lengths):
        q_acc = query_acc.split()[0]
        tot_genome_length += ref_lengths[q_acc]
        # print(t_chr_id, t_start, t_end)
        for ref_id in nams:
            chrom_nams = nams[ref_id]
            total_disjoint_matches += len(chrom_nams)
            # solutions, opt_cov = colinear_solver_read_coverage(chrom_nams, 10000000)
            solutions, opt_cov = n_logn_read_coverage(chrom_nams)
            # pick best from forward and reverse strand
            if q_acc in read_coverage_solution:
                c_prev, _ = read_coverage_solution[q_acc]
                if c_prev < opt_cov:
                    # print(ref_id, opt_cov)
                    read_coverage_solution[q_acc] = (opt_cov, solutions[0])                    
            else:
                # print(ref_id, opt_cov)
                read_coverage_solution[q_acc] = (opt_cov, solutions[0])  


    collinear_chain_nam_sizes = [] 
    total_bp_covered = 0
    if args.collinear_matches_out:
        collinear_outfile = open(args.collinear_matches_out, "w")
    for q_acc in read_coverage_solution:
        opt_cov, solution = read_coverage_solution[q_acc]
        total_bp_covered += opt_cov
        for n in solution:
            collinear_chain_nam_sizes.append(n.y - n.x)

        if args.collinear_matches_out:
            collinear_outfile.write("> {0}\n".format(q_acc))
            for n in solution:
                 collinear_outfile.write("  {0} {1} {2} {3}\n".format(n.chr_id, n.x, n.c, n.val))


    coll_esize = e_size(collinear_chain_nam_sizes, tot_genome_length)
    print("{0} & {1} & {2} & {3} & {4} ".format(total_disjoint_matches, round(total_bp_covered/tot_genome_length,2), round(coll_esize,1), time_in_sec, round(mem_in_mb, 1)))

if __name__ == '__main__':

    parser = argparse.ArgumentParser("Parses alignments and subsamples a fastq file based on alignments.")
    parser.add_argument('matchfile', type=str, help='TSV ')
    parser.add_argument('runtime_file', type=str, help='TSV ')
    parser.add_argument('--refs', type=str, help='Used for ref lengths ')
    parser.add_argument('--collinear_matches_out', type=str, default= '', help='Output mummerlike TSV file with only matches present in collinear chaining solution. ')
    # parser.add_argument('--true_locations', type=str, help='Path to minimap2 SAM alignments to approximate ground truth.')
    # parser.add_argument('--read_vs_read_mode', action="store_true", help='Path to minimap2 genomme SAM alignments to approximate ground truth overlap with read vs read overlaps.')
    # parser.add_argument('--method', type=str, help='Method ')
    # parser.add_argument('--outfile', type=str, help='outfile ')

    args = parser.parse_args()

    main(args)

