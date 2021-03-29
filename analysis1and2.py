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

from modules import indexing, indexing2, help_functions


def get_match_coverage(seq_len, mers, matches, order, span):
    """
        Span is k if kmer, length between first and last position if spaced kmer,
        and span between the start of the first strobe and the last nucleotide of the last strobe.
    """
    if not matches:
        return 0

    if order == 1:
        match_list = sorted([ (p, p+span) for (p, h) in mers.items() if h in matches ])
    elif order == 2:
        match_list = sorted([ (p1,p2+span) for ((p1,p2), h) in mers.items() if h in matches ])
    elif order == 3:
        match_list = sorted([ (p1, p3+span) for ((p1,p2, p3), h) in mers.items() if h in matches ])

    covered_bases = match_list[0][1] - match_list[0][0]
    max_stop = match_list[0][1]
    for start,stop in match_list:
        if start < max_stop and stop > max_stop:
            covered_bases += stop - max_stop
            max_stop = stop
        elif start >= max_stop:
            covered_bases += stop - start
            max_stop = stop
    # print(covered_bases)
    return covered_bases #round(100*covered_bases/seq_len, 1)


def seq_covered_spaced_kmers(mers, matches, seq, positions):
    """
        Function specific to calculate the coverage of spaced k-mers
        since we have all the "sampled" positions in a spaced k-mer
        we can keep an active set of covered positions as we iterate 
        through the string.
    """
    seq_covered = 0
    if not matches:
        return seq_covered

    spaced_kmer_sampled_pos = sorted(positions)
    covered_pos = [0]*len(seq)
    all_match_pos_vector = sorted([p for (p, k) in mers.items() if k in matches ])

    for p in all_match_pos_vector:
        for j in spaced_kmer_sampled_pos:
            covered_pos[p + j] +=1

    # all positions covered at least once
    c = sum([1 for p in covered_pos if p > 0])
    return c


def get_intervals(mers, matches, order):
    if not matches:
        return [], []

    if order == 1:
        all_pos_vector = sorted([p for (p, k) in mers.items() if k in matches ])
    elif order == 2:
        match_list = [ [p1,p2] for ((p1,p2), k) in mers.items() if k in matches ]
        all_pos_vector = sorted([p for sublist in match_list for p in sublist])
    elif order == 3:
        match_list = [ [p1,p2, p3] for ((p1,p2, p3), k) in mers.items() if k in matches ]
        all_pos_vector = sorted([p for sublist in match_list for p in sublist])


    # ivls = []
    # iv_start = all_pos_vector[0]
    # p_prev = all_pos_vector[0]
    # for p in all_pos_vector[1:]:
    #     if p == p_prev + 1:
    #         p_prev += 1
    #     else:
    #         ivls.append((iv_start, p_prev))
    #         p_prev = p
    #         iv_start = p
    # ivls.append((iv_start, p_prev))

    ivls = []
    iv_start = all_pos_vector[0]
    length = 0
    for p1,p2 in zip(all_pos_vector[:-1], all_pos_vector[1:]):
        if p2 == p1 + 1:
            length += 1
        # elif p2 == p1: # for the strobes
        #     pass
        elif p2 > p1 + 1:
            ivls.append((iv_start, iv_start+length))
            length = 0
            iv_start = p2

    if len(all_pos_vector) > 1:
        if p2 <= p1 + 1:
            ivls.append((iv_start, iv_start+length))
    elif len(all_pos_vector) == 1:
        ivls.append((iv_start, iv_start))
    # print(ivls)
    return ivls, all_pos_vector


def statistics(ivls, seq, k):
    if not ivls:
        return 1, [len(seq)], 0
    seq_covered = 0 # lower estimate for minstrobes and randstrobes
    nr_islands = 0
    island_lengths = []

    prev_stop = 0 #ivls[0][1] + k
    for i, (start, stop) in enumerate(ivls):
        if i == 0:
            seq_covered += (stop-start) + k
            if start > 0:
                nr_islands += 1
                island_lengths.append(start - prev_stop)

        elif start - k > prev_stop:
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



def analyze_strobemers(seq1, seq2, k_size, order, hash_fcn, w, w_low = 0, w_high = 50):
    # minstrobes
    if order == 2:
        assert k_size % 2 == 0, "Not even kmer length, results will be different"
        if hash_fcn == "randstrobes":
            strobemers1 = indexing2.randstrobes(seq1, k_size, w_low, w_high, w, order = 2 )
            strobemers2 = indexing2.randstrobes(seq2, k_size, w_low, w_high, w, order = 2 )               
            # print("randstrobes2",  len(strobemers1), len(strobemers2), len(set(strobemers1.values()) & set(strobemers2.values())))
            # print(strobemers1)
            # print(sorted(set(strobemers1.values()))[:20])
            # print(sorted(set(strobemers2.values()))[:20])
        elif hash_fcn == "minstrobes":
            strobemers1 = indexing2.minstrobes(seq1, k_size, w_low, w_high, w, order = 2)
            strobemers2 = indexing2.minstrobes(seq2, k_size, w_low, w_high, w, order = 2)    
            # print("minstrobes2",  len(strobemers2))
        elif hash_fcn == "hybridstrobes":
            strobemers1 = indexing2.hybridstrobes(seq1, k_size, w_low, w_high, w, order = 2)
            strobemers2 = indexing2.hybridstrobes(seq2, k_size, w_low, w_high, w, order = 2)    
            # print("minstrobes2",  len(strobemers2))
    elif order == 3:
        assert k_size % 3 == 0, "Not div by 3 kmer length, results will be different"
        if hash_fcn == "randstrobes":
            strobemers1 = indexing2.randstrobes(seq1, k_size, w_low, w_high, w, order = 3)
            strobemers2 = indexing2.randstrobes(seq2, k_size, w_low, w_high, w, order = 3)
            # print("randstrobes3",  len(strobemers1), len(strobemers2), len(set(strobemers1.values()) & set(strobemers2.values())))
            # print(strobemers1)

        elif hash_fcn == "minstrobes":
            strobemers1 = indexing2.minstrobes(seq1, k_size, w_low, w_high, w, order = 3)
            strobemers2 = indexing2.minstrobes(seq2, k_size, w_low, w_high, w, order = 3)
            # print("minstrobes3",  len(strobemers2))
        
        elif hash_fcn == "hybridstrobes":
            strobemers1 = indexing2.hybridstrobes(seq1, k_size, w_low, w_high, w, order = 3)
            strobemers2 = indexing2.hybridstrobes(seq2, k_size, w_low, w_high, w, order = 3)    
            # print("minstrobes2",  len(strobemers2))

    # elif order == 4:
    #     assert k_size % 4 == 0, "Not div by 4 kmer length, results will be different"
    #     if hash_fcn == "randstrobes":
    #         strobemers1 = indexing2.randstrobes(seq1, k_size, order = 4, w_1 = w_1, w_2 = w_2, w_3 = w_3)
    #         strobemers2 = indexing2.randstrobes(seq2, k_size, order = 4, w_1 = w_1, w_2 = w_2, w_3 = w_3)
    #     # elif hash_fcn == "minstrobes":
    #     #     strobemers1 = indexing2.minstrobes(seq1, k_size, order = 3, w_1 = w_1, w_2 = w_2)
    #     #     strobemers2 = indexing2.minstrobes(seq2, k_size, order = 3, w_1 = w_1, w_2 = w_2 )
    # print(hash_fcn, order, len(strobemers2))
    matches = set(strobemers1.values()) & set(strobemers2.values())
    m = len(matches)
    mp = len(strobemers1.values())
    ivls, all_pos_vector = get_intervals(strobemers1, matches, order)
    nr_islands, island_lengths, c = statistics(ivls, seq1, k_size//order)
    match_coverage = get_match_coverage(len(seq1), strobemers1, matches, order, k_size//order)

    # print("2-spaced minstrobes nr_matches:", len(matches2))  
    # print("Number of islands (gaps):", nr_islands)
    # print("Avg island size (gaps):", sum(island_lengths)/len(island_lengths))
    # print("Seq covered:", seq_covered)
    # print("2-spaced minstrobes intervals:", ivls)
    return m, mp, c, island_lengths, all_pos_vector, match_coverage


def analyze_kmers(seq1, seq2, k_size, w):
    #kmers
    kmers_pos1 = indexing2.kmers(seq1, k_size, w)
    kmers_pos2 = indexing2.kmers(seq2, k_size, w)
    # print("kmers", 1, len(kmers_pos1))
    # print("kmers:",  len(kmers_pos2))
    # kmers_pos1 = {p : seq1[i:i+k_size] for p, i in enumerate(range(len(seq1) - k_size +1))}
    # kmers_seq1 = set([seq1[i:i+k_size] for i in range(len(seq1) - k_size +1)])
    # kmers_seq2 = set([seq2[i:i+k_size] for i in range(len(seq2) - k_size +1)])
    # matches  = kmers_seq1 & kmers_seq2
    matches = set(kmers_pos1.values()) & set(kmers_pos2.values())
    m = len(matches)
    mp = len(kmers_pos1.values())
    ivls, all_pos_vector = get_intervals(kmers_pos1, matches, 1)
    nr_islands, island_lengths, c = statistics(ivls, seq1, k_size)
    match_coverage = get_match_coverage(len(seq1), kmers_pos1, matches, 1, k_size)
    # print("kmers nr_matches:", len(matches))    
    # print("Number of islands (gaps):", nr_islands)
    # print("Avg island size (gaps):", sum(island_lengths)/len(island_lengths))
    # print("Seq covered:", seq_covered)
    # print("kmer intervals:", ivls)
    return m, mp, c, island_lengths, all_pos_vector, match_coverage


def analyze_spaced_kmers(seq1, seq2, k_size, span_size, w):
    positions = set(random.sample(range(1, span_size - 1 ), k_size-2)) 
    positions.add(0)
    positions.add(span_size - 1) # asserts first and last position is sampled so that we have a spaced kmer of length span size
    spaced_kmers_seq1 = indexing2.spaced_kmers(seq1, k_size, span_size, positions, w)
    spaced_kmers_seq2 = indexing2.spaced_kmers(seq2, k_size, span_size, positions, w) 
    matches  = set(spaced_kmers_seq1.values()) & set(spaced_kmers_seq2.values())
    m = len(matches)
    mp = len(spaced_kmers_seq1.values())
    ivls, all_pos_vector = get_intervals(spaced_kmers_seq1, matches, 1)
    nr_islands, island_lengths, _ = statistics(ivls, seq1, span_size)
    # we compute coverage for spaced k-mers with specific function
    c = seq_covered_spaced_kmers(spaced_kmers_seq1, matches, seq1, positions)
    match_coverage = get_match_coverage(len(seq1), spaced_kmers_seq1, matches, 1, span_size)

    # print("kmers nr_matches:", len(matches))    
    # print("Number of islands (gaps):", nr_islands)
    # print("Avg island size (gaps):", sum(island_lengths)/len(island_lengths))
    # print("Seq covered:", seq_covered)
    # print("kmer intervals:", ivls)
    return m, mp, c, island_lengths, all_pos_vector, match_coverage 


def plot_island_distribution(results, mut_freq, outfolder):
    # pd.set_option("display.precision", 8)
    filename = os.path.join(outfolder, "{0}.pdf".format(mut_freq))
    plt.yscale('log', nonposy='clip')
    # bins = [0.1*i for i in range(300)]
    for label in results:
        if label == "kmers" or label == "spaced_kmers_dense" or label == "spaced_kmers_sparse":
            if label == "kmers":
                flat = [g for l in results[label]["islands"] for g in l]
                pyplot.hist(flat, 100, range=[0, 500], alpha=0.5, label=label)
        else:
            for t in results[label]:
                flat = [g for l in results[label][t]["islands"] for g in l]
                tmp_label = label + "-{0}".format(t)
                # if label == "randstrobes" and t == (3,10,25):
                #     pyplot.hist(flat, 100, range=[0, 500], alpha=0.5, label=tmp_label)
                pyplot.hist(flat, 100, range=[0, 500], alpha=0.5, label=tmp_label)

    pyplot.legend(loc='upper right')
    # pyplot.xlabel("Difference to genome (%)")
    pyplot.xlabel("Island length")
    pyplot.ylabel("Count")
    plt.savefig(filename)
    plt.savefig(filename)
    plt.close()

def plot_island_distribution2(results, mut_freq, outfolder):
    # pd.set_option("display.precision", 8)
    filename = os.path.join(outfolder, "{0}.pdf".format(mut_freq))

    # make data correct format
    data = open(os.path.join(outfolder, "{0}.csv".format(mut_freq)), "w")
    data.write("label\tdp\tmut_freq\n")

    for label in results:
        if label == "kmers" or label == "spaced_kmers_dense" or label == "spaced_kmers_sparse":
            flat = [g for l in results[label]["islands"] for g in l]
            for dp in flat:
                data.write("{0}\t{1}\t{2}\n".format(label, dp, mut_freq))
        else:
            for t in results[label]:
                flat = [g for l in results[label][t]["islands"] for g in l]
                tmp_label = label + "-{0}".format(t)
                for dp in flat:
                    data.write("{0}\t{1}\t{2}\n".format(tmp_label, dp, mut_freq))
    data.close()
    # hue_order = ["randstrobes-(3, 10, 25, 50)", "randstrobes-(2, 15, 25, 50)", "kmers", "minstrobes-(3, 10, 25, 50)", "minstrobes-(2, 15, 25, 50)", "spaced_kmers_dense", "spaced_kmers_sparse"]
    hue_order = ["randstrobes-(3, 10, 25, 50)", "hybridstrobes-(3, 10, 25, 50)", "kmers", "minstrobes-(2, 15, 25, 50)", "spaced_kmers_dense"]
    data = pd.read_csv(data.name, sep='\t')
    # plt.yscale('log', nonposy='clip')
    # bins = [0.1*i for i in range(300)]
    sns.displot(data, x="dp", hue="label", hue_order = hue_order,
                 element="step", log_scale=(True, True)) # , kind="kde", log_scale= True, fill=True, multiple="stack"
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
    ax[0].set_title('minstrobes (2,9,10,20)')
    mat = ax[0].imshow(binary_matrices[0], cmap='GnBu', interpolation='nearest')
    # ax[0].set_yticks(range(binary_matrices[0].shape[0]), id_labels)

    ax[1].set_title('minstrobes (3,6,10,20)')
    mat = ax[1].imshow(binary_matrices[1], cmap='GnBu', interpolation='nearest')
    # ax[1].set_yticks(range(binary_matrices[1].shape[0]), id_labels)

    # ax[2].set_yticks(range(binary_matrices[2].shape[0]), id_labels)
    ax[2].set_title('randstrobes (2,9,10,20)')
    mat = ax[2].imshow(binary_matrices[2], cmap='GnBu', interpolation='nearest')

    # ax[3].set_yticks(range(binary_matrices[3].shape[0]), id_labels)
    ax[3].set_title('randstrobes (3,6,10,20)')
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


def get_e_size(all_islands, L, nr_exp):
    # print("all_islands",all_islands)
    sum_of_squares = sum([x**2 for x in all_islands])
    return sum_of_squares/(L*nr_exp) 

def main(args):
    L = 10000
    k_size = 30
    nr_exp = 1000
    w = 1 # thinning, w = 1  means no thinning
    mut_freqs = [0.01, 0.05, 0.1] #[0.1] 
    w_2low = 25
    w_3low = 25
    w_2high = 50
    w_3high = 50

    # w_3strobe = 25
    # w_4strobe = 25
    # experiment_type choose between 'only_subs', 'controlled' or 'all'
    experiment_type = "all" #"controlled" # "all" #"only_subs" # "" # for spaced kmers
    # mut_freq = 0.5 #0.01 #, 0.05, 0.1]
    list_for_illustration = [[],[],[],[]]

    for mut_freq in mut_freqs:
        print("MUTATION RATE:", mut_freq)
        results = {"kmers" : {"m": 0, "mp": 0, "c": 0, "islands": [], "mc":0}, 
                    "spaced_kmers_dense" : {"m": 0, "mp": 0, "c": 0, "islands": [], "mc":0},
                    "spaced_kmers_sparse" : {"m": 0, "mp": 0, "c": 0, "islands": [], "mc":0},
                    "minstrobes" : { (2,15,w_2low,w_2high): {"m": 0, "mp": 0, "c": 0, "islands": [], "mc":0 }, 
                                     (3,10,w_3low,w_3high): {"m": 0, "mp": 0, "c": 0, "islands": [], "mc":0} },
                    "randstrobes" : { (2,15,w_2low,w_2high): {"m": 0, "mp": 0, "c": 0, "islands": [], "mc":0}, 
                                      (3,10,w_3low,w_3high): {"m": 0, "mp": 0, "c": 0, "islands": [], "mc":0} },
                    "hybridstrobes" : { (2,15,w_2low,w_2high): {"m": 0, "mp": 0, "c": 0, "islands": [], "mc":0 },
                                        (3,10,w_3low,w_3high): {"m": 0, "mp": 0, "c": 0, "islands": [], "mc":0 }}                    
                   }
        for exp_id in range(nr_exp):
            seq1 = "".join([random.choice("ACGT") for i in range(L)])
            
            # controlled or random experiment
            if experiment_type == 'only_subs':
                muts = set(random.sample(range(len(seq1)),int(L*mut_freq)))
                seq2 = "".join([seq1[i] if i not in muts else random.choice([help_functions.reverse_complement(seq1[i])]) for i in range(len(seq1))])
            elif experiment_type == 'controlled':
                # muts = set(range(15,L,15)) # every 15th nt for figure 2 only!
                muts = set(range(20,L,20)) 
                seq2 = "".join([seq1[i] if i not in muts else random.choice(['', help_functions.reverse_complement(seq1[i]), seq1[i] + random.choice("ACGT")]) for i in range(len(seq1))])
            elif experiment_type == 'all':
                muts = set(random.sample(range(len(seq1)),int(L*mut_freq)))
                seq2 = "".join([seq1[i] if i not in muts else random.choice(['', help_functions.reverse_complement(seq1[i]), seq1[i] + random.choice("ACGT")]) for i in range(len(seq1))])
            else:
                print("Wrong experiment label specified")
                sys.exit()


            
            # kmers
            m,mp,c,islands,all_pos_vector, match_coverage = analyze_kmers(seq1, seq2, k_size, w)
            results["kmers"]["m"] += m 
            results["kmers"]["c"] += c 
            results["kmers"]["islands"].append(islands) 
            results["kmers"]["mc"] += match_coverage 
            results["kmers"]["mp"] += mp 
            # print("kmers", match_coverage)
            # print_matches(all_pos_vector, "kmers")

            # Spaced kmers dense
            m,mp,c,islands,all_pos_vector, match_coverage = analyze_spaced_kmers(seq1, seq2, k_size, k_size+k_size//2, w)
            results["spaced_kmers_dense"]["m"] += m 
            results["spaced_kmers_dense"]["c"] += c 
            results["spaced_kmers_dense"]["islands"].append(islands) 
            results["spaced_kmers_dense"]["mc"] += match_coverage 
            results["spaced_kmers_dense"]["mp"] += mp 
            # print("spaced_kmers_dense", match_coverage)

            # print_matches(all_pos_vector, "Spaced kmers")

            # Spaced kmers sparse
            m,mp,c,islands,all_pos_vector, match_coverage = analyze_spaced_kmers(seq1, seq2, k_size, 3*k_size, w)
            results["spaced_kmers_sparse"]["m"] += m 
            results["spaced_kmers_sparse"]["c"] += c 
            results["spaced_kmers_sparse"]["islands"].append(islands) 
            results["spaced_kmers_sparse"]["mc"] += match_coverage 
            results["spaced_kmers_sparse"]["mp"] += mp 
            # print("spaced_kmers_sparse", match_coverage)
            # print_matches(all_pos_vector, "Spaced kmers")


            m,mp,c,islands,all_pos_vector, match_coverage = analyze_strobemers(seq1, seq2, k_size, 2, "minstrobes", w , w_low = w_2low, w_high = w_2high)
            results["minstrobes"][(2,15,w_2low,w_2high)]["m"] += m 
            results["minstrobes"][(2,15,w_2low,w_2high)]["c"] += c 
            results["minstrobes"][(2,15,w_2low,w_2high)]["islands"].append(islands) 
            results["minstrobes"][(2,15,w_2low,w_2high)]["mc"] += match_coverage 
            results["minstrobes"][(2,15,w_2low,w_2high)]["mp"] += mp 
            # print_matches(all_pos_vector, "minstrobes2")
            # print("minstrobes2", match_coverage)
            list_for_illustration[0].append(all_pos_vector)
            # print(islands)

            m,mp,c,islands,all_pos_vector, match_coverage = analyze_strobemers(seq1, seq2, k_size, 3, "minstrobes", w ,  w_low = w_3low, w_high = w_3high )
            results["minstrobes"][(3,10,w_3low,w_3high)]["m"] += m 
            results["minstrobes"][(3,10,w_3low,w_3high)]["c"] += c 
            results["minstrobes"][(3,10,w_3low,w_3high)]["islands"].append(islands) 
            results["minstrobes"][(3,10,w_3low,w_3high)]["mc"] += match_coverage 
            results["minstrobes"][(3,10,w_3low,w_3high)]["mp"] += mp 
            # print_matches(all_pos_vector, "minstrobes3") 
            # print("minstrobes3", match_coverage)
            list_for_illustration[1].append(all_pos_vector)
            # print(islands)

            m,mp,c,islands,all_pos_vector, match_coverage = analyze_strobemers(seq1, seq2, k_size, 2, "randstrobes", w ,  w_low = w_2low, w_high = w_2high)
            results["randstrobes"][(2,15,w_2low,w_2high)]["m"] += m 
            results["randstrobes"][(2,15,w_2low,w_2high)]["c"] += c 
            results["randstrobes"][(2,15,w_2low,w_2high)]["islands"].append(islands) 
            results["randstrobes"][(2,15,w_2low,w_2high)]["mc"] += match_coverage 
            results["randstrobes"][(2,15,w_2low,w_2high)]["mp"] += mp 
            # print_matches(all_pos_vector, "randstrobes2") 
            # print("randstrobes2", match_coverage)
            list_for_illustration[2].append(all_pos_vector)
            # print(islands)

            # Tried randstrobe n=3 with w1=17 and w2=40 and it further decreases E-size of gaps over results in paper
            # for higher mutation rates 0.05 and 0.1
            m,mp,c,islands,all_pos_vector, match_coverage = analyze_strobemers(seq1, seq2, k_size, 3, "randstrobes", w ,  w_low = w_3low, w_high = w_3high )
            results["randstrobes"][(3,10,w_3low,w_3high)]["m"] += m 
            results["randstrobes"][(3,10,w_3low,w_3high)]["c"] += c 
            results["randstrobes"][(3,10,w_3low,w_3high)]["islands"].append(islands) 
            results["randstrobes"][(3,10,w_3low,w_3high)]["mc"] += match_coverage 
            results["randstrobes"][(3,10,w_3low,w_3high)]["mp"] += mp 
            # print_matches(all_pos_vector, "randstrobes3") 
            # print("randstrobes3", match_coverage)
            list_for_illustration[3].append(all_pos_vector)
            # print(islands)
            
            m,mp,c,islands,all_pos_vector, match_coverage = analyze_strobemers(seq1, seq2, k_size, 2, "hybridstrobes", w , w_low = w_2low, w_high = w_2high)
            results["hybridstrobes"][(2,15,w_2low,w_2high)]["m"] += m 
            results["hybridstrobes"][(2,15,w_2low,w_2high)]["c"] += c 
            results["hybridstrobes"][(2,15,w_2low,w_2high)]["islands"].append(islands) 
            results["hybridstrobes"][(2,15,w_2low,w_2high)]["mc"] += match_coverage 
            results["hybridstrobes"][(2,15,w_2low,w_2high)]["mp"] += mp 

            m,mp,c,islands,all_pos_vector, match_coverage = analyze_strobemers(seq1, seq2, k_size, 3, "hybridstrobes", w , w_low = w_3low, w_high = w_3high)
            results["hybridstrobes"][(3,10,w_3low,w_3high)]["m"] += m 
            results["hybridstrobes"][(3,10,w_3low,w_3high)]["c"] += c 
            results["hybridstrobes"][(3,10,w_3low,w_3high)]["islands"].append(islands) 
            results["hybridstrobes"][(3,10,w_3low,w_3high)]["mc"] += match_coverage 
            results["hybridstrobes"][(3,10,w_3low,w_3high)]["mp"] += mp 

            # m,c,islands,all_pos_vector = analyze_strobemers(seq1, seq2, 28, 4, "randstrobes", w_1 = 7, w_2 = 10, w_3 = 25 )
            # results["randstrobes"][(4,7,w_4strobe)]["m"] += m 
            # results["randstrobes"][(4,7,w_4strobe)]["c"] += c 
            # results["randstrobes"][(4,7,w_4strobe)]["islands"].append(islands) 
            # # print_matches(all_pos_vector, "randstrobes3") 
            # # list_for_illustration[4].append(all_pos_vector)
            # # print(islands)

            
            # print(len(list_for_illustration))
        
        # plot_matches(list_for_illustration, "m", L, k_size, args.outfolder)

        plot_island_distribution2(results, mut_freq, args.outfolder)

        
        for protocol in results:
            if protocol == "kmers" or protocol == "spaced_kmers_sparse" or protocol == "spaced_kmers_dense":
                flat = [g for l in results[protocol]["islands"] for g in l]
                if flat:
                    # avg_island_len = sum(flat)/len(flat)
                    # print(protocol)
                    e_size = get_e_size(flat, L, nr_exp)
                # else:
                #     avg_island_len = 0
                res = [round(100*results[protocol]["m"]/results[protocol]["mp"], 1), 100*results[protocol]["c"]/(L*nr_exp), 100*results[protocol]["mc"]/(L*nr_exp), e_size]
                print(protocol, " & ".join([ str(round(r, 1)) for r in res]) )
            else:
                for params in results[protocol]:
                    flat = [g for l in results[protocol][params]["islands"] for g in l]
                    if flat:
                        # avg_island_len = sum(flat)/len(flat)
                        # print(protocol, params)
                        e_size = get_e_size(flat, L, nr_exp)
                    # else:
                        # avg_island_len = 0
                    res = [round(100*results[protocol][params]["m"]/results[protocol][params]["mp"], 1), 100*results[protocol][params]["c"]/(L*nr_exp), 100*results[protocol][params]["mc"]/(L*nr_exp), e_size]
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