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


def get_subsamples(transcript_cov, coverage):
    subsamples = {}

    large_enough_cov = [tr_id for tr_id in transcript_cov.keys() if len(transcript_cov[tr_id]) >= coverage]
    print("Transcripts with more than sample size coverage:", len(large_enough_cov))
    for tr_id in large_enough_cov:
        read_subset = random.sample(transcript_cov[tr_id], coverage)
        subsamples[tr_id] = read_subset

    return subsamples


def get_abundance_aligned_reads(sam_file, sample_size):
    SAM_file = pysam.AlignmentFile(sam_file, "r", check_sq=False)
    references = SAM_file.references
    
    transcript_cov = defaultdict(set)
    # amgiguous_primary = defaultdict(set)
    for read in SAM_file.fetch(until_eof=True):
        if (read.flag == 0 or read.flag == 16) and read.mapping_quality >= 10:
            transcript_id = read.reference_name
            transcript_cov[transcript_id].add(read.query_name)

        # elif (read.flag == 0 or read.flag == 16) and read.mapping_quality < 10:
        #     transcript_id = read.reference_name
        #     amgiguous_primary[transcript_id].add(read.query_name)

    return transcript_cov #, amgiguous_primary

def main(args):
    transcript_cov = get_abundance_aligned_reads(args.alignments, args.sample_size)
    min_cov = min([len(v) for v in transcript_cov.values()])
    print("Minimum SIRV coverage is:", min_cov)

    subsamples = get_subsamples(transcript_cov,  args.sample_size)

    fastq = { acc : (seq,qual) for acc, (seq,qual) in readfq(open(args.fastq, 'r'))}
    ref_fasta = { acc : seq for acc, (seq,qual) in readfq(open(args.ref_fasta, 'r'))}
    for tr_id in subsamples:
        ref_outfile = open(os.path.join(args.outfolder, str(tr_id) + '_ref.fastq'), "w")
        ref_outfile.write(">{0}\n{1}\n".format(tr_id, ref_fasta[tr_id]))

        reads_outfile = open(os.path.join(args.outfolder, str(tr_id) + '.fastq'), "w")
        sampled_reads = subsamples[tr_id]
        for acc in sampled_reads:
            seq, qual = fastq[acc]
            reads_outfile.write("@{0}\n{1}\n+\n{2}\n".format(acc, seq, qual))
        reads_outfile.close()
        ref_outfile.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser("Parses alignments and subsamples a fastq file based on alignments.")
    parser.add_argument('fastq', type=str, help='fastq file. ')
    parser.add_argument('ref_fasta', type=str, help='fastq file. ')
    parser.add_argument('alignments', type=str, help='fastq file. ')
    parser.add_argument('outfolder', type=str, help='Fastq file. ')
    parser.add_argument('sample_size', type=int, help='sample_size ')
    # parser.add_argument('max_size', type=int, help='Min size of reads.')

    args = parser.parse_args()

    main(args)

