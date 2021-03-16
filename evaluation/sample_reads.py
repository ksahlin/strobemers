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


def reverse_complement(string):
    rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'t', 'c':'g', 'g':'c', 't':'a', 'N':'N', 'X':'X', 'n':'n', 'Y':'R', 'R':'Y', 'K':'M', 'M':'K', 'S':'S', 'W':'W', 'B':'V', 'V':'B', 'H':'D', 'D':'H', 'y':'r', 'r':'y', 'k':'m', 'm':'k', 's':'s', 'w':'w', 'b':'v', 'v':'b', 'h':'d', 'd':'h'}
    rev_comp = ''.join([rev_nuc[nucl] for nucl in reversed(string)])
    return(rev_comp)

def get_subsamples(transcript_cov, coverage):
    subsamples = {}

    large_enough_cov = [tr_id for tr_id in transcript_cov.keys() if len(transcript_cov[tr_id]) >= coverage]
    print("Transcripts with more than sample size coverage:", len(large_enough_cov))
    for tr_id in large_enough_cov:
        read_subset = random.sample(transcript_cov[tr_id], coverage)
        subsamples[tr_id] = read_subset

    return subsamples


def get_abundance_aligned_full_length_reads(sam_file, sample_size, ref_fasta):
    SAM_file = pysam.AlignmentFile(sam_file, "r", check_sq=False)
    references = SAM_file.references
    
    transcript_cov = defaultdict(set)
    # amgiguous_primary = defaultdict(set)
    for read in SAM_file.fetch(until_eof=True):
        ref_acc = read.reference_name
        # print(ref_acc)
        if ref_acc:
            ref_len = len(ref_fasta[ref_acc])
        if read.flag == 0 and read.mapping_quality > 0 and read.reference_start < 10 and (ref_len - read.reference_end) < 10 :
            transcript_id = read.reference_name
            transcript_cov[transcript_id].add(read.query_name)
        elif read.flag == 16 and read.mapping_quality > 0 and read.reference_start < 10 and (ref_len - read.reference_end) < 10: 
            transcript_id = read.reference_name
            transcript_cov[transcript_id + '_rev'].add(read.query_name)

        # elif (read.flag == 0 or read.flag == 16) and read.mapping_quality < 10:
        #     transcript_id = read.reference_name
        #     amgiguous_primary[transcript_id].add(read.query_name)

    return transcript_cov #, amgiguous_primary

def main(args):

    ref_fasta = { acc : seq for acc, (seq,qual) in readfq(open(args.ref_fasta, 'r'))}
    transcript_cov = get_abundance_aligned_full_length_reads(args.alignments, args.sample_size, ref_fasta)
    min_cov = min([len(v) for v in transcript_cov.values()])
    for k, v in transcript_cov.items():
        print(k, len(v))
    print("Minimum SIRV coverage is:", min_cov)

    subsamples = get_subsamples(transcript_cov,  args.sample_size)

    fastq = { acc : (seq,qual) for acc, (seq,qual) in readfq(open(args.fastq, 'r'))}
    is_rc = False
    for tr_id in subsamples:
        if tr_id[-4:] == "_rev":
            continue
            # is_rc = True
            # ref_seq = reverse_complement(ref_fasta[tr_id])
        else:
            is_rc = False
            ref_seq = ref_fasta[tr_id]

        # if is_rc:
        #     ref_acc = tr_id[-4:]
        # else:
        ref_acc = tr_id
        ref_outfile = open(os.path.join(args.outfolder, str(ref_acc) + '_ref.fastq'), "w")
        ref_outfile.write(">{0}\n{1}\n".format(ref_acc, ref_seq))

        reads_outfile = open(os.path.join(args.outfolder, str(ref_acc) + '.fastq'), "w")
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

