#!/usr/bin/env python

import argparse
import sys, os
import random
import pysam

from collections import defaultdict
from collections import namedtuple


NAM = namedtuple('NAM', ['x', 'y', 'c', 'd', 'val', 'j', "exon_part_id"])


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
    print("READ {0} LINES.".format(i))
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



def main(args):
    if args.reads:
        read_lengths = { acc : len(seq) for acc, (seq,qual) in readfq(open(args.reads, 'r'))}

    read_coverage_solution = {}
    read_tot_hits = {}
    for (read_acc, nams) in get_NAM_records(args.infile, read_lengths):
        r_acc = read_acc.split()[0]
        # print(read_acc, r_acc)

        for chrom in nams:
            chrom_nams = nams[chrom]
            solutions, opt_cov = colinear_solver_read_coverage(chrom_nams, 10000000)

            # pick best from forward and reverse strand
            if r_acc in read_coverage_solution:
                c_prev = read_coverage_solution[r_acc]
                if c_prev < opt_cov:
                    read_coverage_solution[r_acc] = opt_cov
            else:
                read_coverage_solution[r_acc] = opt_cov


            # pick largest nr of hits
            if r_acc in read_tot_hits:
                hits_prev = read_tot_hits[r_acc]
                if hits_prev < len(chrom_nams):
                    read_tot_hits[r_acc] = len(chrom_nams)
            else:
                read_tot_hits[r_acc] = len(chrom_nams)


    out = open(args.outfile, "w")
    out.write("method,r_acc,r_len,r_frac_cov,r_hits\n".format())
    for r_acc in read_coverage_solution:
        r_len = read_lengths[r_acc]
        r_cov = read_coverage_solution[r_acc]
        r_hits = read_tot_hits[r_acc]
        r_frac_cov = round(float(r_cov)/r_len, 2)
        out.write("{0},{1},{2},{3},{4}\n".format(args.method, r_acc, r_len, r_frac_cov, r_hits))
    out.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser("Parses alignments and subsamples a fastq file based on alignments.")
    parser.add_argument('infile', type=str, help='TSV ')
    parser.add_argument('--reads', type=str, help='Used for read lengths ')
    parser.add_argument('--method', type=str, help='Method ')
    parser.add_argument('--outfile', type=str, help='outfile ')

    args = parser.parse_args()

    main(args)

