import os,sys
import argparse

import random

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except (ImportError, RuntimeError):
    print("COULD not import matplotlib")

import seaborn as sns
import pandas as pd
from matplotlib import pyplot

from modules import indexing, help_functions

def get_intervals(mers, matches):
    if not matches:
        return [], []
    all_pos_vector = [p for (p, k) in mers.items() if k in matches ]
    ivls = []
    iv_start = all_pos_vector[0]
    p_prev = all_pos_vector[0]
    for p in all_pos_vector[1:]:
        if p == p_prev + 1:
            p_prev += 1
        else:
            ivls.append((iv_start, p_prev))
            p_prev = p
            iv_start = p
    ivls.append((iv_start, p_prev))
    return ivls, all_pos_vector

def statistics(ivls, seq, k):
    if not ivls:
        return 0, [len(seq)], 0
    seq_covered = 0 # lower estimate for minstrobes and randstrobes
    nr_islands = 0
    island_lengths = []

    prev_stop = 0 #ivls[0][1] + k
    for start, stop in ivls:
        if start - k > prev_stop:
            seq_covered += (stop-start) + k
            nr_islands += 1
            island_lengths.append(start - prev_stop)
            assert start - prev_stop >= 0, "found: {0},{1},{2}: {3}".format(start, prev_stop, k, ivls)
        else: # overlapping previous hit 
            seq_covered += (stop-start) + k - (k - (start - prev_stop) )

        prev_stop = stop

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
    ivls, all_pos_vector = get_intervals(strobemers1, matches)
    nr_islands, island_lengths, c = statistics(ivls, seq1, k_size//order)
    # print("2-spaced minstrobes nr_matches:", len(matches2))  
    # print("Number of islands (gaps):", nr_islands)
    # print("Avg island size (gaps):", sum(island_lengths)/len(island_lengths))
    # print("Seq covered:", seq_covered)
    # print("2-spaced minstrobes intervals:", ivls)
    return m, c, island_lengths, all_pos_vector


def analyze_kmers(seq1, seq2, k_size):
    #kmers
    kmers_pos1 = {p : seq1[i:i+k_size] for p, i in enumerate(range(len(seq1) - k_size +1))}
    kmers_seq1 = set([seq1[i:i+k_size] for i in range(len(seq1) - k_size +1)])
    kmers_seq2 = set([seq2[i:i+k_size] for i in range(len(seq2) - k_size +1)])
    matches  = kmers_seq1 & kmers_seq2
    m = len(matches)
    ivls, all_pos_vector = get_intervals(kmers_pos1, matches)
    nr_islands, island_lengths, c = statistics(ivls, seq1, k_size)
    # print("kmers nr_matches:", len(matches))    
    # print("Number of islands (gaps):", nr_islands)
    # print("Avg island size (gaps):", sum(island_lengths)/len(island_lengths))
    # print("Seq covered:", seq_covered)
    # print("kmer intervals:", ivls)
    return m, c, island_lengths, all_pos_vector


def plot_island_distribution(results, mut_freq, outfolder):
    # pd.set_option("display.precision", 8)
    filename = os.path.join(outfolder, "{0}.pdf".format(mut_freq))
    plt.yscale('log', nonposy='clip')
    # bins = [0.1*i for i in range(300)]
    for label in results:
        if label == "kmers":
            flat = [g for l in results[label]["islands"] for g in l]
            pyplot.hist(flat, 100, range=[0, 500], alpha=0.5, label=label)
        else:
            for t in results[label]:
                flat = [g for l in results[label][t]["islands"] for g in l]
                tmp_label = label + "-{0}".format(t)
                if label == "randstrobes" and t == (3,10,25):
                    pyplot.hist(flat, 100, range=[0, 500], alpha=0.5, label=tmp_label)

    pyplot.legend(loc='upper right')
    # pyplot.xlabel("Difference to genome (%)")
    pyplot.xlabel("Island length")
    pyplot.ylabel("Count")
    plt.savefig(filename)
    plt.savefig(filename)
    plt.close()

def print_matches(all_pos_vector, method):
    s = set(all_pos_vector)
    for i in range(100):
        if i in s:
            print("X",end='')
        else:
            print(" ",end='')
    print(method)
    print()

import numpy as np
def plot_matches(all_data, method, L, k_size,outfolder):

    data = np.random.randn(5, 2)
    print(data)
    binary_matrices = []
    for all_runs in all_data:
        binary_matrix = []
        for run in all_runs:
            binary_vector = []
            s = set(run)
            for i in range(L - k_size +1):
                if i in s:
                    # print("X",end='')
                    binary_vector.append(1)
                else:
                    # print(" ",end='')
                    binary_vector.append(0)
            binary_matrix.append(binary_vector)
        binary_matrices.append( np.array(binary_matrix)  )
    # print(binary_matrix)

    # np_matrix = np.array(binary_matrix)  
    fig, ax = plt.subplots(4,sharex=True, sharey=True)
    fig.suptitle('Match distribution')
    plt.yticks([])
    id_labels = ["1", "2", "3", "4", "5"]
    ax[0].set_title('minstrobes (2,9,40)')
    mat = ax[0].imshow(binary_matrices[0], cmap='GnBu', interpolation='nearest')
    # ax[0].set_yticks(range(binary_matrices[0].shape[0]), id_labels)

    ax[1].set_title('minstrobes (3,6,20)')
    mat = ax[1].imshow(binary_matrices[1], cmap='GnBu', interpolation='nearest')
    # ax[1].set_yticks(range(binary_matrices[1].shape[0]), id_labels)

    # ax[2].set_yticks(range(binary_matrices[2].shape[0]), id_labels)
    ax[2].set_title('randstrobes (2,9,40)')
    mat = ax[2].imshow(binary_matrices[2], cmap='GnBu', interpolation='nearest')

    # ax[3].set_yticks(range(binary_matrices[3].shape[0]), id_labels)
    ax[3].set_title('randstrobes (3,6,20)')
    mat = ax[3].imshow(binary_matrices[3], cmap='GnBu', interpolation='nearest')


    # plt.xticks(range(id_matrix.shape[1]), concert_dates)
    # plt.xticks(rotation=30)
    plt.xlabel('Position')

    # # this places 0 or 1 centered in the individual squares
    # for x in range(np_matrix.shape[0]):
    #     for y in range(np_matrix.shape[1]):
    #         ax.annotate(str(np_matrix[x, y])[0], xy=(y, x), 
    #                     horizontalalignment='center', verticalalignment='center')

    # ax = sns.heatmap(data, cbar=False, xticklabels = False, yticklabels=False)
    # ax.tick_params(left=False, bottom=False)
    filename = os.path.join(outfolder, "{0}_ex.pdf".format(method))
    plt.tight_layout()
    plt.savefig(filename)


def main(args):
    # reads = { acc: seq for i, (acc, (seq, _)) in enumerate(readfq(open(args.fasta, 'r')))}
    # seq1, seq2 = list(reads.values())
    L = 100
    k_size = 18
    nr_exp = 5
    mut_freqs = [0.1] # [0.01, 0.05, 0.1]
    # mut_freq = 0.5 #0.01 #, 0.05, 0.1]
    list_for_illustration = [[],[],[],[]]

    for mut_freq in mut_freqs:
        print("MUTATION RATE:", mut_freq)
        results = {"kmers" : {"m": 0, "c": 0, "islands": []},
                    "minstrobes" : { (2,15,50): {"m": 0, "c": 0, "islands": []}, (3,10,25): {"m": 0, "c": 0, "islands": []} },
                    "randstrobes" : { (2,15,50): {"m": 0, "c": 0, "islands": []}, (3,10,25): {"m": 0, "c": 0, "islands": []} }
                   }
        for exp_id in range(nr_exp):
            seq1 = "".join([random.choice("ACGT") for i in range(L)])
            
            # controlled or random experiment
            muts = set(range(20,L,15)) 
            # muts = set(random.sample(range(len(seq1)),int(L*mut_freq)))

            seq2 = "".join([seq1[i] if i not in muts else random.choice(['', help_functions.reverse_complement(seq1[i]), seq1[i] + random.choice("ACGT")]) for i in range(len(seq1))])
            
            #kmers
            m,c,islands,all_pos_vector = analyze_kmers(seq1, seq2, k_size)
            results["kmers"]["m"] += m 
            results["kmers"]["c"] += c 
            results["kmers"]["islands"].append(islands) 
            print_matches(all_pos_vector, "kmers")
            m,c,islands,all_pos_vector = analyze_strobemers(seq1, seq2, k_size, 2, "minstrobes", N_1 = 50 )
            results["minstrobes"][(2,15,50)]["m"] += m 
            results["minstrobes"][(2,15,50)]["c"] += c 
            results["minstrobes"][(2,15,50)]["islands"].append(islands) 
            print_matches(all_pos_vector, "minstrobes2")
            list_for_illustration[0].append(all_pos_vector)
            # print(islands)

            m,c,islands,all_pos_vector = analyze_strobemers(seq1, seq2, k_size, 3, "minstrobes", N_1 = 25, N_2 = 25 )
            results["minstrobes"][(3,10,25)]["m"] += m 
            results["minstrobes"][(3,10,25)]["c"] += c 
            results["minstrobes"][(3,10,25)]["islands"].append(islands) 
            print_matches(all_pos_vector, "minstrobes3") 
            list_for_illustration[1].append(all_pos_vector)
            # print(islands)

            m,c,islands,all_pos_vector = analyze_strobemers(seq1, seq2, k_size, 2, "randstrobes", N_1 = 50 )
            results["randstrobes"][(2,15,50)]["m"] += m 
            results["randstrobes"][(2,15,50)]["c"] += c 
            results["randstrobes"][(2,15,50)]["islands"].append(islands) 
            print_matches(all_pos_vector, "randstrobes2") 
            list_for_illustration[2].append(all_pos_vector)
            # print(islands)

            m,c,islands,all_pos_vector = analyze_strobemers(seq1, seq2, k_size, 3, "randstrobes", N_1 = 25, N_2 = 25 )
            results["randstrobes"][(3,10,25)]["m"] += m 
            results["randstrobes"][(3,10,25)]["c"] += c 
            results["randstrobes"][(3,10,25)]["islands"].append(islands) 
            print_matches(all_pos_vector, "randstrobes3") 
            list_for_illustration[3].append(all_pos_vector)
            # print(islands)
            print(len(list_for_illustration))
        
        plot_matches(list_for_illustration, "m", L, k_size, args.outfolder)

        plot_island_distribution(results, mut_freq, args.outfolder)

        
        for protocol in results:
            if protocol == "kmers":
                flat = [g for l in results[protocol]["islands"] for g in l]
                avg_island_len = sum(flat)/len(flat)
                res = [round(100*results[protocol]["m"]/(L*nr_exp), 1), 100*results[protocol]["c"]/(L*nr_exp), avg_island_len]
                print(protocol, " & ".join([ str(round(r, 1)) for r in res]) )
            else:
                for params in results[protocol]:
                    flat = [g for l in results[protocol][params]["islands"] for g in l]
                    avg_island_len = sum(flat)/len(flat)
                    res = [round(100*results[protocol][params]["m"]/(L*nr_exp), 1), 100*results[protocol][params]["c"]/(L*nr_exp), avg_island_len]
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
    parser.add_argument('--outfolder', type=str,  default=None, help='A fasta file with transcripts that are shared between samples and have perfect illumina support.')
    # parser.add_argument('--pickled_subreads', type=str, help='Path to an already parsed subreads file in pickle format')
    # parser.set_defaults(which='main')
    args = parser.parse_args()



    # if len(sys.argv)==1:
    #     parser.print_help()
    #     sys.exit()

    main(args)