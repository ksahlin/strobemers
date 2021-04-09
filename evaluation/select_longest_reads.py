#!/usr/bin/env python
import sys
import argparse
import pysam


'''
    Below fastq reader function taken from https://github.com/lh3/readfq/blob/master/readfq.py
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


def select_aligned_reads(sam_file, reads_fastq):
    SAM_file = pysam.AlignmentFile(sam_file, "r", check_sq=False)
    references = SAM_file.references
    error_rates = {}
    aligned = {}
    for read in SAM_file.fetch(until_eof=True):
        if read.flag == 0 or read.flag == 16:
            read_acc = read.query_name
            read_len = len(reads_fastq[read_acc][0])
            if float(read.query_alignment_length) / read_len > 0.95:
                aligned[read_acc] = reads_fastq[read_acc]

                ins = sum([length for type_, length in read.cigartuples if type_ == 1])
                del_ = sum([length for type_, length in read.cigartuples if type_ == 2])
                subs = sum([length for type_, length in read.cigartuples if type_ == 8])
                # matches = sum([length for type_, length in read.cigartuples if type_ == 7])
                # tot_align = ins + del_ + subs + matches
                # identity = matches/float(tot_align)
                # print(float(ins + del_ + subs), read.query_alignment_length)
                error_rates[read_acc] = float(ins + del_ + subs)/ read.query_alignment_length

    print(sum(error_rates.values())/len(error_rates), min(error_rates.values()), max(error_rates.values()) )
    return aligned

def main(args):
    fastq = { acc : (seq, qual) for acc, (seq,qual) in readfq(open(args.fastq, 'r'))}
    aligned = select_aligned_reads(args.alignments, fastq)

    outfile = open(args.outfile, "w")
    i = 0
    tot_len = 0
    for acc, (seq, qual) in sorted(aligned.items(), key = lambda x: len(x[1][0]), reverse = True):
        outfile.write("@{0}\n{1}\n+\n{2}\n".format(acc, seq, qual))
        tot_len += len(seq)
        if i > args.sample_size:
            break
        i += 1
    outfile.close()
    print("tot len:", tot_len)

if __name__ == '__main__':

    parser = argparse.ArgumentParser("Parses alignments and subsamples a fastq file based on alignments.")
    parser.add_argument('fastq', type=str, help='fastq file. ')
    parser.add_argument('alignments', type=str, help='SAM file. ')
    parser.add_argument('sample_size', type=int, help='sample_size ')
    parser.add_argument('outfile', type=str, help='Fastq file. ')

    args = parser.parse_args()

    main(args)
