#!/usr/bin/env python

import argparse
import sys, os
import random
import pysam

from collections import defaultdict



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


# results_file has format: method,read_id,ref_id,hit_len, transcript_length? if we want to normalize

def main(args):
    if args.refs:
        ref_lengths = [ (acc, len(seq)) for acc, (seq,qual) in readfq(open(args.refs, 'r'))]
        assert len(ref_lengths) == 1
        ref_len = ref_lengths[0][1]
        ref_acc = ref_lengths[0][0]

    coverage_file = open( os.path.join(args.outfolder, 'coverage.csv'), "a+")
    normalized_hit_length_file = open( os.path.join(args.outfolder, 'normalized_hit_length.csv'), "a+")
    nr_hits_file = open( os.path.join(args.outfolder, 'nr_hits.csv'), "a+")

    for i, line in enumerate(open(args.infile, 'r')):
        if line[0] == ">":
            if i == 0:
                nr_hits = 0
                coverage = 0
            else:
                nr_hits_file.write(",".join([str(x) for x in [args.method,read_acc, ref_acc, ref_len, nr_hits ]]) + "\n")
                coverage_file.write(",".join([str(x) for x in [args.method,read_acc, ref_acc, ref_len, coverage ]]) + "\n")
                nr_hits = 0
                coverage = 0

            read_acc = line[1:].strip()
            continue
        else: 
            ref, r_pos, q_pos, match_length = line.split()
            normalized_hit_length_file.write(",".join([str(x) for x in [args.method,read_acc, ref_acc, ref_len, match_length, int(match_length)/ref_len ]]) + "\n")
            nr_hits += 1
            coverage += int(match_length)/ref_len





if __name__ == '__main__':

    parser = argparse.ArgumentParser("Parses alignments and subsamples a fastq file based on alignments.")
    parser.add_argument('infile', type=str, help='TSV ')
    parser.add_argument('--refs', type=str, help='Used for lengths ')
    parser.add_argument('--method', type=str, help='Method ')
    parser.add_argument('--outfolder', type=str, help='outfolder ')

    args = parser.parse_args()

    main(args)

