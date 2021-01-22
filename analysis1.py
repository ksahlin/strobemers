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

def statistics(ivls, seq, k):
    if not ivls:
        return 0, [0], 0
    seq_covered = 0 # lower estimate for strobemers and randomers
    nr_islands = 0
    island_lengths = []
    if ivls[0][0] > 0:
        nr_islands +=1
        island_lengths.append(ivls[0][0])

    prev_stop = ivls[0][1]
    for start, stop in ivls:
        seq_covered += (stop-start+1) + min(k, k - (prev_stop - start))
        nr_islands += 1
        island_lengths.append(start - prev_stop)
        prev_stop = stop

    if ivls[-1][0] < len(seq):
        nr_islands +=1
        island_lengths.append(len(seq) - prev_stop)

    return nr_islands, island_lengths, seq_covered

def positions_matching_kmers(seq1, seq2, k_size):

    # print("length seq1:", len(seq1))
    # print("length seq2:", len(seq2))

    #kmers
    kmers_pos1 = {p : seq1[i:i+k_size] for p, i in enumerate(range(len(seq1) - k_size +1))}
    kmers_seq1 = set([seq1[i:i+k_size] for i in range(len(seq1) - k_size +1)])
    kmers_seq2 = set([seq2[i:i+k_size] for i in range(len(seq2) - k_size +1)])
    matches  = kmers_seq1 & kmers_seq2
    ivls = get_intervals(kmers_pos1, matches)
    nr_islands, island_lengths, seq_covered = statistics(ivls, seq1, k_size)
    print("kmers nr_matches:", len(matches))    
    print("Number of islands (gaps):", nr_islands)
    print("Avg island size (gaps):", sum(island_lengths)/len(island_lengths))
    print("Seq covered:", seq_covered)
    # print("kmer intervals:", ivls)

    # randomers
    assert k_size % 2 == 0, "Not even kmer length, results will be different"
    randomers1 = indexing.randomers(seq1, k_size, order = 2, N_1 = 50 )
    randomers2 = indexing.randomers(seq2, k_size, order = 2, N_1 = 50 )    
    matches2 = set(randomers1.values()) & set(randomers2.values())
    ivls = get_intervals(randomers1, matches2)
    nr_islands, island_lengths, seq_covered = statistics(ivls, seq1, k_size//2)
    print("2-spaced randomers nr_matches:", len(matches2))  
    print("Number of islands (gaps):", nr_islands)
    print("Avg island size (gaps):", sum(island_lengths)/len(island_lengths))
    print("Seq covered:", seq_covered)
    # print("2-spaced randomers intervals:", ivls)


    assert k_size % 3 == 0, "Not div by 3 kmer length, results will be different"
    randomers1 = indexing.randomers(seq1, k_size, order = 3, N_1 = 25, N_2 = 25)
    randomers2 = indexing.randomers(seq2, k_size, order = 3, N_1 = 25, N_2 = 25 )
    matches3 = set(randomers1.values()) & set(randomers2.values())
    ivls = get_intervals(randomers1, matches3)
    nr_islands, island_lengths, seq_covered = statistics(ivls, seq1, k_size//3)
    print("3-spaced randomers nr_matches:", len(matches3))  
    print("Number of islands (gaps):", nr_islands)
    print("Avg island size (gaps):", sum(island_lengths)/len(island_lengths))
    print("Seq covered:", seq_covered)
    # print("3-spaced randomers intervals:", ivls)


    # strobemers
    assert k_size % 2 == 0, "Not even kmer length, results will be different"
    strobemers1 = indexing.strobemers(seq1, k_size, order = 2, N_1 = 50 )
    strobemers2 = indexing.strobemers(seq2, k_size, order = 2, N_1 = 50 )    
    matches2 = set(strobemers1.values()) & set(strobemers2.values())
    ivls = get_intervals(strobemers1, matches2)
    nr_islands, island_lengths, seq_covered = statistics(ivls, seq1, k_size//2)
    print("2-spaced strobemers nr_matches:", len(matches2))  
    print("Number of islands (gaps):", nr_islands)
    print("Avg island size (gaps):", sum(island_lengths)/len(island_lengths))
    print("Seq covered:", seq_covered)
    # print("2-spaced strobemers intervals:", ivls)


    assert k_size % 3 == 0, "Not div by 3 kmer length, results will be different"
    strobemers1 = indexing.strobemers(seq1, k_size, order = 3, N_1 = 25, N_2 = 25)
    strobemers2 = indexing.strobemers(seq2, k_size, order = 3, N_1 = 25, N_2 = 25 )
    matches3 = set(strobemers1.values()) & set(strobemers2.values())
    ivls = get_intervals(strobemers1, matches3)
    nr_islands, island_lengths, seq_covered = statistics(ivls, seq1, k_size//3)
    print("3-spaced strobemers nr_matches:", len(matches3)) 
    print("Number of islands (gaps):", nr_islands)
    print("Avg island size (gaps):", sum(island_lengths)/len(island_lengths))
    print("Seq covered:", seq_covered)
    # print("3-spaced strobemers intervals:", ivls)

def main(args):
    # reads = { acc: seq for i, (acc, (seq, _)) in enumerate(readfq(open(args.fasta, 'r')))}
    # seq1, seq2 = list(reads.values())
    L = 1000
    k = 30

    # # controlled experiment
    # for exp_id in range(10):
    #     print()
    #     print()
    #     seq1 = "".join([random.choice("ACGT") for i in range(L)])
    #     muts = set(range(20,1000,20)) 
    #     seq2 = "".join([seq1[i] if i not in muts else random.choice(['', help_functions.reverse_complement(seq1[i]), seq1[i] + random.choice("ACGT")]) for i in range(len(seq1))])
    #     positions_matching_kmers(seq1, seq2, k)

    # random mutation mositions
    for mut_freq in [0.01, 0.05, 0.1]:
        for exp_id in range(10):
            seq1 = "".join([random.choice("ACGT") for i in range(L)])
            muts = set(random.sample(range(len(seq1)),int(L*mut_freq)))
            # muts = set(range(20,1000,20)) #set([20,40,60,80])
            # print(muts)
            seq2 = "".join([seq1[i] if i not in muts else random.choice(['', help_functions.reverse_complement(seq1[i]), seq1[i] + random.choice("ACGT")]) for i in range(len(seq1))])
            print()
            print("MUT FREQ:", mut_freq)
            positions_matching_kmers(seq1, seq2, k)



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