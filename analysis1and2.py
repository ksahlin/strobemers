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
        else:
            ivls.append((iv_start, p_prev))
            p_prev = p
            iv_start = p
    ivls.append((iv_start, p_prev))
    return ivls

def statistics(ivls, seq, k):
    if not ivls:
        return 0, [len(seq)], 0
    seq_covered = 0 # lower estimate for minstrobes and randstrobes
    nr_islands = 0
    island_lengths = []

    prev_stop = 0 #ivls[0][1] + k
    for start, stop in ivls:
        if start > prev_stop:
            seq_covered += (stop-start) + k
            nr_islands += 1
            island_lengths.append(start - prev_stop)
            assert start - prev_stop >= 0, "found: {0},{1},{2}: {3}".format(start, prev_stop, k, ivls)
        else: # overlapping previous hit 
            seq_covered += (stop-start) + k - (prev_stop - start)

        prev_stop = stop + k

    if ivls[-1][0] + k < len(seq):
        nr_islands +=1
        island_lengths.append(len(seq) - prev_stop)

    return nr_islands, island_lengths, seq_covered

# def positions_matching_kmers(seq1, seq2, k_size):

#     # randstrobes
#     assert k_size % 2 == 0, "Not even kmer length, results will be different"
#     randomers1 = indexing.randstrobes(seq1, k_size, order = 2, N_1 = 50 )
#     randomers2 = indexing.randstrobes(seq2, k_size, order = 2, N_1 = 50 )    
#     matches2 = set(randomers1.values()) & set(randomers2.values())
#     ivls = get_intervals(randomers1, matches2)
#     nr_islands, island_lengths, seq_covered = statistics(ivls, seq1, k_size//2)
#     print("2-spaced randstrobes nr_matches:", len(matches2))  
#     print("Number of islands (gaps):", nr_islands)
#     print("Avg island size (gaps):", sum(island_lengths)/len(island_lengths))
#     print("Seq covered:", seq_covered)
#     # print("2-spaced randstrobes intervals:", ivls)


#     assert k_size % 3 == 0, "Not div by 3 kmer length, results will be different"
#     randomers1 = indexing.randstrobes(seq1, k_size, order = 3, N_1 = 25, N_2 = 25)
#     randomers2 = indexing.randstrobes(seq2, k_size, order = 3, N_1 = 25, N_2 = 25 )
#     matches3 = set(randomers1.values()) & set(randomers2.values())
#     ivls = get_intervals(randomers1, matches3)
#     nr_islands, island_lengths, seq_covered = statistics(ivls, seq1, k_size//3)
#     print("3-spaced randstrobes nr_matches:", len(matches3))  
#     print("Number of islands (gaps):", nr_islands)
#     print("Avg island size (gaps):", sum(island_lengths)/len(island_lengths))
#     print("Seq covered:", seq_covered)
#     # print("3-spaced randstrobes intervals:", ivls)

#     assert k_size % 3 == 0, "Not div by 3 kmer length, results will be different"
#     strobemers1 = indexing.minstrobes(seq1, k_size, order = 3, N_1 = 25, N_2 = 25)
#     strobemers2 = indexing.minstrobes(seq2, k_size, order = 3, N_1 = 25, N_2 = 25 )
#     matches3 = set(strobemers1.values()) & set(strobemers2.values())
#     ivls = get_intervals(strobemers1, matches3)
#     nr_islands, island_lengths, seq_covered = statistics(ivls, seq1, k_size//3)
#     print("3-spaced minstrobes nr_matches:", len(matches3)) 
#     print("Number of islands (gaps):", nr_islands)
#     print("Avg island size (gaps):", sum(island_lengths)/len(island_lengths))
#     print("Seq covered:", seq_covered)
#     # print("3-spaced minstrobes intervals:", ivls)



def analyze_strobemers(seq1, seq2, k_size, order, hash_fcn, N_1 = 50, N_2 = 50 ):
    # minstrobes
    if order == 2:
        assert k_size % 2 == 0, "Not even kmer length, results will be different"
        if hash_fcn == "randstrobes":
            strobemers1 = indexing.randstrobes(seq1, k_size, order = 2, N_1 = N_1 )
            strobemers2 = indexing.randstrobes(seq2, k_size, order = 2, N_1 = N_1 )               
        elif hash_fcn == "minstrobes":
            strobemers1 = indexing.minstrobes(seq1, k_size, order = 2, N_1 = N_1 )
            strobemers2 = indexing.minstrobes(seq2, k_size, order = 2, N_1 = N_1 )    
    elif order == 3:
        assert k_size % 3 == 0, "Not div by 3 kmer length, results will be different"
        if hash_fcn == "randstrobes":
            strobemers1 = indexing.randstrobes(seq1, k_size, order = 3, N_1 = N_1, N_2 = N_2)
            strobemers2 = indexing.randstrobes(seq2, k_size, order = 3, N_1 = N_1, N_2 = N_2 )
        elif hash_fcn == "minstrobes":
            strobemers1 = indexing.minstrobes(seq1, k_size, order = 3, N_1 = N_1, N_2 = N_2)
            strobemers2 = indexing.minstrobes(seq2, k_size, order = 3, N_1 = N_1, N_2 = N_2 )

    matches = set(strobemers1.values()) & set(strobemers2.values())
    m = len(matches)
    ivls = get_intervals(strobemers1, matches)
    nr_islands, island_lengths, c = statistics(ivls, seq1, k_size//order)
    # print("2-spaced minstrobes nr_matches:", len(matches2))  
    # print("Number of islands (gaps):", nr_islands)
    # print("Avg island size (gaps):", sum(island_lengths)/len(island_lengths))
    # print("Seq covered:", seq_covered)
    # print("2-spaced minstrobes intervals:", ivls)
    return m, c, island_lengths


def analyze_kmers(seq1, seq2, k_size):
    #kmers
    kmers_pos1 = {p : seq1[i:i+k_size] for p, i in enumerate(range(len(seq1) - k_size +1))}
    kmers_seq1 = set([seq1[i:i+k_size] for i in range(len(seq1) - k_size +1)])
    kmers_seq2 = set([seq2[i:i+k_size] for i in range(len(seq2) - k_size +1)])
    matches  = kmers_seq1 & kmers_seq2
    m = len(matches)
    ivls = get_intervals(kmers_pos1, matches)
    nr_islands, island_lengths, c = statistics(ivls, seq1, k_size)
    # print("kmers nr_matches:", len(matches))    
    # print("Number of islands (gaps):", nr_islands)
    # print("Avg island size (gaps):", sum(island_lengths)/len(island_lengths))
    # print("Seq covered:", seq_covered)
    # print("kmer intervals:", ivls)
    return m, c, island_lengths


def main(args):
    # reads = { acc: seq for i, (acc, (seq, _)) in enumerate(readfq(open(args.fasta, 'r')))}
    # seq1, seq2 = list(reads.values())
    L = 10000
    k_size = 30
    nr_exp = 1000
    # mut_freq = 0.5 #0.01 #, 0.05, 0.1]


    for mut_freq in [0.01, 0.05, 0.1]:
        print("MUTATION RATE:", mut_freq)
        results = {"kmers" : {"m": 0, "c": 0, "islands": []},
                    "minstrobes" : { (2,15,50): {"m": 0, "c": 0, "islands": []}, (3,10,25): {"m": 0, "c": 0, "islands": []} },
                    "randstrobes" : { (2,15,50): {"m": 0, "c": 0, "islands": []}, (3,10,25): {"m": 0, "c": 0, "islands": []} }
                   }
        for exp_id in range(nr_exp):
            seq1 = "".join([random.choice("ACGT") for i in range(L)])
            
            # controlled or random experiment
            # muts = set(range(20,10000,20)) 
            muts = set(random.sample(range(len(seq1)),int(L*mut_freq)))

            seq2 = "".join([seq1[i] if i not in muts else random.choice(['', help_functions.reverse_complement(seq1[i]), seq1[i] + random.choice("ACGT")]) for i in range(len(seq1))])
            
            #kmers
            m,c,islands = analyze_kmers(seq1, seq2, k_size)
            results["kmers"]["m"] += m 
            results["kmers"]["c"] += c 
            results["kmers"]["islands"].append(islands) 

            m,c,islands = analyze_strobemers(seq1, seq2, k_size, 2, "minstrobes", N_1 = 50 )
            results["minstrobes"][(2,15,50)]["m"] += m 
            results["minstrobes"][(2,15,50)]["c"] += c 
            results["minstrobes"][(2,15,50)]["islands"].append(islands) 

            m,c,islands = analyze_strobemers(seq1, seq2, k_size, 3, "minstrobes", N_1 = 25, N_2 = 25 )
            results["minstrobes"][(3,10,25)]["m"] += m 
            results["minstrobes"][(3,10,25)]["c"] += c 
            results["minstrobes"][(3,10,25)]["islands"].append(islands) 

            m,c,islands = analyze_strobemers(seq1, seq2, k_size, 2, "randstrobes", N_1 = 50 )
            results["randstrobes"][(2,15,50)]["m"] += m 
            results["randstrobes"][(2,15,50)]["c"] += c 
            results["randstrobes"][(2,15,50)]["islands"].append(islands) 

            m,c,islands = analyze_strobemers(seq1, seq2, k_size, 3, "randstrobes", N_1 = 25, N_2 = 25 )
            results["randstrobes"][(3,10,25)]["m"] += m 
            results["randstrobes"][(3,10,25)]["c"] += c 
            results["randstrobes"][(3,10,25)]["islands"].append(islands) 

            # positions_matching_kmers(seq1, seq2, k)
        for protocol in results:
            if protocol == "kmers":
                flat = [g for l in results[protocol]["islands"] for g in l]
                avg_island_len = sum(flat)/len(flat)
                res = [results[protocol]["m"]/nr_exp, 100*results[protocol]["c"]/(L*nr_exp), avg_island_len]
                print(protocol, " & ".join([ str(round(r, 1)) for r in res]) )
            else:
                for params in results[protocol]:
                    flat = [g for l in results[protocol][params]["islands"] for g in l]
                    avg_island_len = sum(flat)/len(flat)
                    res = [results[protocol][params]["m"]/nr_exp, 100*results[protocol][params]["c"]/(L*nr_exp), avg_island_len]
                    print(protocol, params, " & ".join([ str(round(r, 1)) for r in res]) )

    # print(results)


    # # random mutation mositions
    # for mut_freq in [0.01, 0.05, 0.1]:
    #     for exp_id in range(10):
    #         seq1 = "".join([random.choice("ACGT") for i in range(L)])
    #         muts = set(random.sample(range(len(seq1)),int(L*mut_freq)))
    #         # muts = set(range(20,1000,20)) #set([20,40,60,80])
    #         # print(muts)
    #         seq2 = "".join([seq1[i] if i not in muts else random.choice(['', help_functions.reverse_complement(seq1[i]), seq1[i] + random.choice("ACGT")]) for i in range(len(seq1))])
    #         print()
    #         print("MUT FREQ:", mut_freq)
    #         positions_matching_kmers(seq1, seq2, k)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calc identity", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # parser.add_argument('--fasta', type=str,  default=False, help='Path to consensus fastq file(s)')
    # parser.add_argument('--k', type=int, default=13, help='Kmer size')
    # parser.add_argument('--w', type=int, default=20, help='Window size')
    # parser.add_argument('--outfolder', type=str,  default=None, help='A fasta file with transcripts that are shared between samples and have perfect illumina support.')
    # parser.add_argument('--pickled_subreads', type=str, help='Path to an already parsed subreads file in pickle format')
    # parser.set_defaults(which='main')
    args = parser.parse_args()



    # if len(sys.argv)==1:
    #     parser.print_help()
    #     sys.exit()

    main(args)