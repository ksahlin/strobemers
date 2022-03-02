#!/usr/bin/env python3.9
# -*- coding: utf-8 -*
import sys
import operator
import random
import copy
from collections import defaultdict, deque
from fractions import Fraction
from typing import Iterator

BITS = sys.hash_info.width
MAX = sys.maxsize
MAX_HASH_VALUE = int((2**BITS)/2) - 1


def Sequence(L: int) -> str:
    """
    Generate a random canonical DNA sequence.

    :param L: an integer representing the desired sequence length
    :returns: a string with a random DNA sequence of length L
    """
    DNA = "".join([random.choice("ACGT") for i in range(L)])
    return DNA


def argmin(array: list) -> tuple:
    """
    Find the value of x which minimizes f(x) over the set of candidates for x

    :param array: a list to minimize
    :returns: a tuple with the index position and the value of the lowest element
    """
    min_index = array.index(min(array))
    min_val = array[min_index]
    return min_index, min_val


def thinner(hash_list: list, w: int) -> list:
    """
    Thins out kmers/strobemers using a sliding window approach

    :param hash_list: a list with hash values
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :returns: a list with tuples (pos in original list, minimim hash value) for each window of w hashes
    """
    window_hashes = deque(hash_list[:w])
    min_index, curr_min_hash = argmin(window_hashes)
    thinned_hash_list = [(min_index, curr_min_hash)]

    for i in range(w, len(hash_list) + w-1):
        if i >= len(hash_list):
            new_hash = MAX
        else:
            new_hash = hash_list[i]

        # updating window
        discarded_hash = window_hashes.popleft()
        window_hashes.append(new_hash)

        # we have discarded previous windows minimizer, look for new minimizer brute force
        if curr_min_hash == discarded_hash:
            min_index, curr_min_hash = argmin(window_hashes)
            thinned_hash_list.append((min_index + i + 1 - w, curr_min_hash))

        # previous minimizer still in window, we only need to compare with the recently added kmer
        elif new_hash < curr_min_hash:
            curr_min_hash = new_hash
            thinned_hash_list.append((i, curr_min_hash))

    return thinned_hash_list


def kmers(seq: str, k_size: int, w: int) -> dict:
    """
    Sample a substrings of length k contained within a given sequence

    :param seq: a string with a nucleotide sequence
    :param k_size: length of the kmer
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :returns: a dictionary with positions along the string as keys and the kmers as value
    """
    if w > 1:
        hash_seq_list = [(i, hash(seq[i:i+k_size])) for i in range(len(seq) - k_size + 1)]
        # produce a subset of positions, still with same index as in full sequence
        hash_seq_list_thinned = thinner([h for i, h in hash_seq_list], w)
        kmers_pos = {i: h for i, h in hash_seq_list_thinned}
    else:
        kmers_pos = {i: hash(seq[i:i+k_size]) for i in range(len(seq) - k_size + 1)}

    return kmers_pos


def kmer_iter(seq: str, k_size: int, w: int) -> Iterator[tuple]:
    """
    Iterator for creating kmers

    :param seq: a string with a nucleotide sequence
    :param k_size: length of the kmer
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :returns: an iterator for creating kmers
    """
    if w > 1:
        hash_seq_list = [(i, hash(seq[i:i+k_size])) for i in range(len(seq) - k_size + 1)]
        # produce a subset of positions, still with samme index as in full sequence
        hash_seq_list_thinned = thinner([h for i, h in hash_seq_list], w)
        for p, h in hash_seq_list_thinned:
            yield p, h
    else:
        hash_seq_list = [(i, hash(seq[i:i+k_size])) for i in range(len(seq) - k_size + 1)]
        for p, h in hash_seq_list:
            yield p, h


def minimizers(seq: str, k_size: int, w_size: int) -> list:
    """
    Sample the smallest k-mer by hash-value in a pre-defined ordering of each k-mer in the window

    :param seq: a string with a nucleotide sequence
    :param k_size: length of the minimizers
    :param w_size: window size
    :returns: a list with minimizers
    """
    # kmers = [seq[i:i+k_size] for i in range(len(seq)-k_size) ]
    w = w_size - k_size
    window_kmers = deque([seq[i:i+k_size] for i in range(w+1)])
    curr_min = min(window_kmers)
    minimizers = [(curr_min, list(window_kmers).index(curr_min))]

    for i in range(w+1, len(seq) - k_size):
        new_kmer = seq[i:i+k_size]
        # updateing window
        discarded_kmer = window_kmers.popleft()
        window_kmers.append(new_kmer)

        # we have discarded previous windows minimizer, look for new minimizer brute force
        if curr_min == discarded_kmer:
            curr_min = min(window_kmers)
            minimizers.append((curr_min, list(window_kmers).index(curr_min) + i - w))

        # previous minimizer still in window, we only need to compare with the recently added kmer
        elif new_kmer < curr_min:
            curr_min = new_kmer
            minimizers.append((curr_min, i))

    return minimizers


def spaced_kmers(seq: str, k_size: int, span_size: int, positions: set,
                 w: int) -> dict:
    """
    Sample kmers with spaced seeds

    :param seq: a string with a nucleotide sequence
    :param k_size: length of the strobe
    :param span_size: length between first and last position
    :param positions: a set of positions to consider for the spaced k-mer
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :returns: a dictionary with positions along the string as keys and the spaced kmer as value
    """
    assert len(positions) == k_size

    if w > 1:
        hash_seq_list = [
            (i, hash("".join([seq[i + j] for j in range(span_size) if j in positions])))
            for i in range(len(seq) - span_size + 1)
        ]
        # produce a subset of positions, still with samme index as in full sequence
        hash_seq_list_thinned = thinner([h for i, h in hash_seq_list], w)
        spaced_kmers_pos = {i: h for i, h in hash_seq_list_thinned}
    else:
        spaced_kmers_pos = {
            i: hash("".join([seq[i + j] for j in range(span_size) if j in positions]))
            for i in range(len(seq) - span_size + 1)
        }
    # print(positions, len(positions), span_size)
    # well, this is not the most time efficient way to sample spaced kmers but works for simulations...
    return spaced_kmers_pos


def spaced_kmers_iter(seq: str, k_size: int, span_size: int,
                      positions: set) -> Iterator[int]:
    """
    Iterator for generating spaced kmers

    :param seq: a string with a nucleotide sequence
    :param k_size: length of the kmers
    :param span_size: length between first and last position
    :param positions: a set of positions to consider for the spaced k-mer
    :returns: an iterator for creating spaced_kmers
    """
    assert len(positions) == k_size
    # print(positions, len(positions), span_size)
    # well, this is not the most time efficient way to sample spaced kmers but works for simulations...
    for i in range(len(seq) - span_size + 1):
        yield hash("".join([seq[i + j] for j in range(span_size) if j in positions]))


def seq_to_randstrobes_iter(seq: str, k_size: int, strobe_w_min_offset: int,
                            strobe_w_max_offset: int, prime: int, w: int,
                            order: int) -> Iterator[tuple]:
    """
    Iterator for creation of randstrobes of any order

    :param seq: a string with a nucleotide sequence
    :param k_size: length of each strobe
    :param strobe_w_min_offset: minimum window offset to the previous window (wMin > 0)
    :param strobe_w_max_offset: maximum window offset to the previous window (wMin <= wMax)
    :param prime: prime number (q) in minimizing h(m)+h(mj) mod q
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :returns: an iterator for creating randstrobes
    """
    hash_seq_list = [(i, hash(seq[i:i+k_size])) for i in range(len(seq) - k_size + 1)]
    # thinning
    if w > 1:
        # produce a subset of positions, still with same index as in full sequence
        hash_seq_list_thinned = thinner([h for i, h in hash_seq_list], w)
    else:
        hash_seq_list_thinned = hash_seq_list

    for (p1, hash_m1) in hash_seq_list_thinned:  # [:-k_size]:
        if p1 >= len(hash_seq_list) - (order-1)*k_size:
            break
        # hash_m1 = hash_seq_list[p]

        if p1 + (order-1) * strobe_w_max_offset <= len(hash_seq_list):
            windows = list()
            for window_order in range(1, order):
                start = p1 + strobe_w_min_offset + (window_order-1) * strobe_w_max_offset
                end = min(p1 + window_order * strobe_w_max_offset, len(hash_seq_list))
                windows.append((start, end))

        else:
            windows = list()
            for window_order in range(1, order):
                start = (max(
                    p1+window_order*k_size,
                    len(hash_seq_list) + strobe_w_min_offset - (order - window_order) * strobe_w_max_offset
                    )
                )

                end = min(p1 + window_order * strobe_w_max_offset, len(hash_seq_list))
                windows.append((start, end))

        index = [p1, ]
        min_values = []
        min_hash_val = hash_m1
        for index_order in range(1, order):
            min_index, min_value = argmin([
                (min_hash_val + hash_seq_list[i][1]) % prime
                for i in range(*windows[index_order-1])
            ])

            min_hash_val = min_hash_val + (index_order * (-1)**index_order) * hash_seq_list[windows[index_order-1][0] + min_index][1]
            index.append(min_index+windows[index_order-1][0])

        yield index, min_hash_val


def randstrobes(seq: str, k_size: int, strobe_w_min_offset: int,
                strobe_w_max_offset: int, w: int, order: int = 2) -> dict:
    """
    Strobemer seeding protocol to sample randstrobes

    :param seq: a string with a nucleotide sequence
    :param k_size: length of all strobes (len(strobe_1) +  ... + len(strobe_n))
    :param strobe_w_min_offset: minimum window offset to the previous window (wMin > 0)
    :param strobe_w_max_offset: maximum window offset to the previous window (wMin <= wMax)
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :param return_min_value: a bool to specify whether the hash value we minimize in the function deciding how to pick the strobe should be returned
    :returns: a dictionary with positions along the string as keys and the randstrobes as value
    """
    randstrobes = dict()
    randstrobes_hash = dict()
    prime = 997
    assert strobe_w_min_offset > 0, "Minimum strobemer offset has to be greater than 0 in this implementation"

    if k_size % order != 0:
        print("WARNING: kmer size is not evenly divisible with {0}, will use {1} as kmer size: ".format(order, k_size - k_size % order))
        k_size = k_size - k_size % order
    m_size = k_size//order

    randstrobes = {tuple(index): h for index, h in seq_to_randstrobes_iter(
        seq, m_size, strobe_w_min_offset, strobe_w_max_offset, prime, w, order
    )}
    return randstrobes


def randstrobes_iter(seq: str, k_size: int, strobe_w_min_offset: int,
                     strobe_w_max_offset: int, w: int, order: int = 2,
                     buffer_size: int = 10000000) -> Iterator[tuple]:
    """
    Generator for creating randstrobes (less memory requiring)

    :param seq: a string with a nucleotide sequence
    :param k_size: length of all strobes (len(strobe_1) +  ... + len(strobe_n))
    :param strobe_w_min_offset: minimum window offset to the previous window (wMin > 0)
    :param strobe_w_max_offset: maximum window offset to the previous window (wMin <= wMax)
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :param buffer_size: size of buffer at the end of the sequence
    :returns: an iterator for creating randstrobes
    """
    for i in range(0, len(seq), buffer_size):
        substring = seq[i:i+buffer_size]
        for p, m in randstrobes(
                substring, k_size, strobe_w_min_offset, strobe_w_max_offset,
                w, order=order).items():

            yield p, m


def seq_to_mixedrandstrobes_iter(seq: str, k_size: int, strobe_w_min_offset: int,
                                 strobe_w_max_offset: int, prime: int, w: int,
                                 order: int, denominator: int,
                                 numerator: int) -> Iterator[tuple]:
    """
    Iterator for creating of mixedrandstrobes of any order

    :param seq: a string with a nucleotide sequence
    :param k_size: length of each strobe
    :param strobe_w_min_offset: minimum window offset to the previous window (wMin > 0)
    :param strobe_w_max_offset: maximum window offset to the previous window (wMin <= wMax)
    :param prime: prime number (q) in minimizing h(m)+h(mj) mod q
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :param denominator: denominator and numerator determine the fraction of sampled strobemers
    :param numerator: denominator and numerator determine the fraction of sampled strobemers
    :returns: an iterator for creating mixedrandstrobes
    """
    hash_seq_list = [(i, hash(seq[i:i+k_size])) for i in range(len(seq) - k_size + 1)]
    # thinning
    if w > 1:
        # produce a subset of positions, still with samme index as in full sequence
        hash_seq_list_thinned = thinner([h for i, h in hash_seq_list], w)
    else:
        hash_seq_list_thinned = hash_seq_list

    for (p1, hash_m1) in hash_seq_list_thinned:  # [:-k_size]:
        if p1 >= len(hash_seq_list) - (order-1)*k_size:
            break
        # hash_m1 = hash_seq_list[p]

        if hash_m1 % denominator < numerator:  # pick randstrobe
            if p1 + (order-1) * strobe_w_max_offset <= len(hash_seq_list):
                windows = list()
                for window_order in range(1, order):
                    start = p1 + strobe_w_min_offset + (window_order-1) * strobe_w_max_offset
                    end = min(p1 + window_order * strobe_w_max_offset, len(hash_seq_list))
                    windows.append((start, end))

            else:
                windows = list()
                for window_order in range(1, order):
                    start = (max(
                        p1+window_order*k_size,
                        len(hash_seq_list) + strobe_w_min_offset - (order - window_order) * strobe_w_max_offset
                        )
                    )

                    end = min(p1 + window_order * strobe_w_max_offset, len(hash_seq_list))
                    windows.append((start, end))

            index = [p1, ]
            min_hash_val = hash_m1
            for index_order in range(1, order):
                min_index, min_value = argmin([
                    (min_hash_val + hash_seq_list[i][1]) % prime
                    for i in range(*windows[index_order-1])
                ])

                min_hash_val = min_hash_val + (index_order * (-1)**index_order) * hash_seq_list[windows[index_order-1][0] + min_index][1]
                index.append(min_index+windows[index_order-1][0])
            yield index, min_hash_val

        else:  # pick k-mer
            index = tuple(p1 + (strobe_num) * k_size for strobe_num in range(order))
            yield index, hash(seq[p1: p1+k_size*order])


def mixedrandstrobes(seq: str, k_size: int, strobe_w_min_offset: int,
                     strobe_w_max_offset: int, w: int,
                     order: int = 2, strobe_fraction: float = 0.5) -> dict:
    """
    Mixed protocol to produce specified randstrobes and k-mers content fractions

    :param seq: a string with a nucleotide sequence
    :param k_size: length of all strobes (len(strobe_1) +  ... + len(strobe_n))
    :param strobe_w_min_offset: minimum window offset to the previous window (wMin > 0)
    :param strobe_w_max_offset: maximum window offset to the previous window (wMin <= wMax)
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :param strobe_fraction: a fraction of sampled strobemers, rest kmers (0 <= strobe_fraction <= 1)
    :returns: a dictionary with positions along the string as keys and the miaxed randstrobes/kmers as value
    """
    mixed_output = dict()
    fraction = Fraction(str(strobe_fraction))
    denominator = fraction.denominator
    numerator = fraction.numerator
    prime = 997
    assert strobe_w_min_offset > 0, "Minimum strobemer offset has to be greater than 0 in this implementation"

    if k_size % order != 0:
        print("WARNING: kmer size is not evenly divisible with {0}, will use {1} as kmer size: ".format(order, k_size - k_size % order))
        k_size = k_size - k_size % order
    m_size = k_size//order

    mixedrandstrobes = {
        tuple(index): h
        for index, h in seq_to_mixedrandstrobes_iter(
            seq, m_size, strobe_w_min_offset, strobe_w_max_offset, prime, w, order, denominator, numerator
        )
    }

    return mixedrandstrobes


def mixedrandstrobes_iter(seq: str, k_size: int, strobe_w_min_offset: int,
                         strobe_w_max_offset: int, w: int, order: int = 2,
                         strobe_fraction: float = 1,
                         buffer_size: int = 10000000) -> Iterator[tuple]:
    """
    Generator for creating mixedrandstrobes (less memory requiring)

    :param seq: a string with a nucleotide sequence
    :param k_size: length of all strobes (len(strobe_1) +  ... + len(strobe_n))
    :param strobe_w_min_offset: minimum window offset to the previous window (wMin > 0)
    :param strobe_w_max_offset: maximum window offset to the previous window (wMin <= wMax)
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :param buffer_size: size of buffer at the end of the sequence
    :returns: an iterator for creating mixedrandstrobes
    """
    for i in range(0, len(seq), buffer_size):
        substring = seq[i: i+buffer_size]
        for p, m in mixedrandstrobes(substring, k_size, strobe_w_min_offset, strobe_w_max_offset, w, order, strobe_fraction).items():
            yield p, m


def seq_to_minstrobes_iter(seq: str, k_size: int, strobe_w_min_offset: int,
                           strobe_w_max_offset: int, prime: int, w: int,
                           order: int) -> Iterator[tuple]:
    """
    Generator for creating minstrobes of any order

    :param seq: a string with a nucleotide sequence
    :param k_size: length of the strobe
    :param strobe_w_min_offset: minimum window offset to the previous window (wMin > 0)
    :param strobe_w_max_offset: maximum window offset to the previous window (wMin <= wMax)
    :param prime: prime number (q) in minimizing h(m)+h(mj) mod q
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :returns: an iterator for creating minstrobes
    """
    hash_seq_list = [(i, hash(seq[i: i+k_size])) for i in range(len(seq) - k_size + 1)]

    # produce a subset of positions, still with samme index as in full sequence
    strobes = deque(thinner([h for i, h in hash_seq_list], strobe_w_max_offset - strobe_w_min_offset))
    strobes_dict = {strobe_num: copy.deepcopy(strobes) for strobe_num in range(1, order)}

    if w > 1:
        # produce a subset of positions, still with same index as in full sequence
        hash_seq_list_thinned = thinner([h for i, h in hash_seq_list], w)
    else:
        hash_seq_list_thinned = hash_seq_list

    for (p1, hash_m1) in hash_seq_list_thinned:  # [:-m_size]:
        if p1 >= len(hash_seq_list) + k_size - k_size*order:
            break

        positions = [p1, ]
        hash_value = hash_m1

        for strobe_num in range(1, order):
            if p1 + k_size + strobe_w_min_offset + (strobe_num-1) * strobe_w_max_offset < len(seq):
                while strobes_dict[strobe_num][0][0] < min(p1 + k_size + strobe_w_min_offset + (strobe_num-1) * strobe_w_max_offset, len(hash_seq_list)-1):
                    l = strobes_dict[strobe_num].popleft()
            p, hash_val = strobes_dict[strobe_num][0]
            positions.append(p)
            hash_value += strobe_num * (-1)**strobe_num * hash_val

        yield positions, hash_value


def minstrobes(seq: str, k_size: int, strobe_w_min_offset: int,
               strobe_w_max_offset: int, w: int, order: int = 2) -> dict:
    """
    Strobemer seeding protocol to sample minstrobes

    :param seq: a string with a nucleotide sequence
    :param k_size: length of all strobes (len(strobe_1) +  ... + len(strobe_n))
    :param strobe_w_min_offset: minimum window offset to the previous window (wMin > 0)
    :param strobe_w_max_offset: maximum window offset to the previous window (wMin <= wMax)
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :returns: a dictionary with positions along the string as keys and the minstrobes as value
    """
    prime = 997
    assert strobe_w_min_offset > 0, "Minimum strobemer offset has to be greater than 0 in this implementation"

    if k_size % order != 0:
        print("WARNING: kmer size is not evenly divisible with {0}, will use {0} as kmer size: ".format(order, k_size - k_size % 2))
        k_size = k_size - k_size % order
    m_size = k_size//order
    assert m_size + (order-1) * strobe_w_max_offset < len(seq), "Last minstrobes window position is exceeding the sequence length, consider to use a lower order"

    minstrobes = {
        tuple(positions): h
        for positions, h in seq_to_minstrobes_iter(
            seq, m_size, strobe_w_min_offset, strobe_w_max_offset, prime,
            w, order
        )
    }

    return minstrobes


def minstrobes_iter(seq: str, k_size: int, strobe_w_min_offset: int,
                    strobe_w_max_offset: int, w: int, order: int = 2,
                    buffer_size: int = 10000000) -> Iterator[tuple]:
    """
    Generator for creating minstrobes (less memory requiring)

    :param seq: a string with a nucleotide sequence
    :param k_size: length of all strobes (len(strobe_1) +  ... + len(strobe_n))
    :param strobe_w_min_offset: minimum window offset to the previous window (wMin > 0)
    :param strobe_w_max_offset: maximum window offset to the previous window (wMin <= wMax)
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :param buffer_size: size of buffer at the end of the sequence
    :returns: an iterator for creating minstrobes
    """
    for i in range(0, len(seq), buffer_size):
        substring = seq[i: i+buffer_size]
        for p, m in minstrobes(substring, k_size, strobe_w_min_offset, strobe_w_max_offset, w, order=order).items():
            yield p, m


def seq_to_mixedminstrobes_iter(seq: str, k_size: int, strobe_w_min_offset: int,
                                strobe_w_max_offset: int, prime: int, w: int,
                                order: int, denominator: int,
                                numerator: int) -> Iterator[tuple]:
    """
    Generator for creating mixedminstrobes of any order

    :param seq: a string with a nucleotide sequence
    :param k_size: length of the strobe
    :param strobe_w_min_offset: minimum window offset to the previous window (wMin > 0)
    :param strobe_w_max_offset: maximum window offset to the previous window (wMin <= wMax)
    :param prime: prime number (q) in minimizing h(m)+h(mj) mod q
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :param denominator: denominator and numerator determine the fraction of sampled strobemers
    :param numerator: denominator and numerator determine the fraction of sampled strobemers
    :returns: a tuple with positions as first element and hash_value as second element.
    """
    hash_seq_list = [(i, hash(seq[i:i+k_size])) for i in range(len(seq) - k_size + 1)]
    # produce a subset of positions, still with samme index as in full sequence
    strobes = deque(thinner([h for i, h in hash_seq_list], strobe_w_max_offset - strobe_w_min_offset))
    strobes_dict = {strobe_num: copy.deepcopy(strobes) for strobe_num in range(1, order)}

    if w > 1:
        # produce a subset of positions, still with same index as in full sequence
        hash_seq_list_thinned = thinner([h for i, h in hash_seq_list], w)
    else:
        hash_seq_list_thinned = hash_seq_list

    # Decision based on hash values
    for (p1, hash_m1) in hash_seq_list_thinned:  # [:-k_size]:
        if p1 >= len(hash_seq_list) + k_size - order*k_size:
            break
        if hash_m1 % denominator < numerator:  # pick minstrobe
            positions = [p1, ]
            hash_value = hash_m1

            for strobe_num in range(1, order):
                if p1 + k_size + strobe_w_min_offset + (strobe_num-1) * strobe_w_max_offset < len(seq):
                    while strobes_dict[strobe_num][0][0] < min(p1 + k_size + strobe_w_min_offset + (strobe_num-1) * strobe_w_max_offset, len(hash_seq_list)-1):
                        l = strobes_dict[strobe_num].popleft()
                p, hash_val = strobes_dict[strobe_num][0]
                positions.append(p)
                hash_value += strobe_num * (-1)**strobe_num * hash_val
            yield positions, hash_value

        else:  # pick k-mer
            index = tuple(p1 + (strobe_num) * k_size for strobe_num in range(order))
            yield index, hash(seq[p1:p1+k_size*order])


def mixedminstrobes(seq: str, k_size: int, strobe_w_min_offset: int,
                    strobe_w_max_offset: int, w: int,
                    order: int = 2, strobe_fraction: float = 0.5) -> dict:
    """
    Mixed protocol to produce specified minstrobes and k-mers content fractions

    :param seq: a string with a nucleotide sequence
    :param k_size: length of all strobes (len(strobe_1) +  ... + len(strobe_n))
    :param strobe_w_min_offset: minimum window offset to the previous window (wMin > 0)
    :param strobe_w_max_offset: maximum window offset to the previous window (wMin <= wMax)
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :param strobe_fraction: a fraction of sampled strobemers, rest kmers (0 <= strobe_fraction <= 1)
    :returns: a dictionary with positions along the string as keys and the mixed minstrobes/kmers as value
    """
    prime = 997
    fraction = Fraction(str(strobe_fraction))
    denominator = fraction.denominator
    numerator = fraction.numerator
    assert strobe_w_min_offset > 0, "Minimum strobemer offset has to be greater than 0 in this implementation"

    if k_size % order != 0:
        print("WARNING: kmer size is not evenly divisible with {0}, will use {0} as kmer size: ".format(order, k_size - k_size % 2))
        k_size = k_size - k_size % order
    m_size = k_size//order
    assert m_size + (order-1) * strobe_w_max_offset < len(seq), "Last minstrobes window position is exceeding the sequence length, consider to use a lower order"

    mixedminstrobes = {
        tuple(positions): h
        for positions, h in seq_to_mixedminstrobes_iter(
                seq, m_size, strobe_w_min_offset, strobe_w_max_offset, prime, w,
                order, denominator, numerator
            )
        }
    return mixedminstrobes


def mixedminstrobes_iter(seq: str, k_size: int, strobe_w_min_offset: int,
                         strobe_w_max_offset: int, w: int, order: int = 2,
                         strobe_fraction: float = 1,
                         buffer_size: int = 10000000) -> Iterator[tuple]:
    """
    Generator for creating mixedminstrobes (less memory requiring)

    :param seq: a string with a nucleotide sequence
    :param k_size: length of all strobes (len(strobe_1) +  ... + len(strobe_n))
    :param strobe_w_min_offset: minimum window offset to the previous window (wMin > 0)
    :param strobe_w_max_offset: maximum window offset to the previous window (wMin <= wMax)
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :param buffer_size: size of buffer at the end of the sequence
    :returns: an iterator for creating minstrobes
    """
    for i in range(0, len(seq), buffer_size):
        substring = seq[i: i+buffer_size]
        for p, m in mixedminstrobes(substring, k_size, strobe_w_min_offset, strobe_w_max_offset, w, order, strobe_fraction).items():
            yield p, m


def update_queue(q: list, curr_min: int, min_index: int, new_hash: int, i: int,
                 start_offset: int, end_offset: int) -> tuple:
    """
    Updates windows

    :param q: a list with strobe_windows
    :param curr_min: an integer with the current minimum value
    :param min_index: an integer with the index position of the minimum value
    :param new_hash: an integer with the new hash value
    :param i: an integer with the position of the first strobe
    :param start_offset: minimum window offset
    :param end_offset: maximum window offset
    :returns: a tuple with the index position of the minimum value and the value
    """
    old_h = q.popleft()
    q.append(new_hash)

    # we have discarded previous windows minimizer, look for new minimizer brute force
    if curr_min == old_h:
        min_index, curr_min = argmin(q)
        min_index = i + start_offset + min_index

    # Previous minimizer still in window, we only need to compare with the recently added kmer
    elif new_hash < curr_min:
        curr_min = new_hash
        min_index = i + end_offset

    return min_index, curr_min


def seq_to_hybridstrobes_iter(seq: str, k_size: int, w_min, w_max, w: int,
                              order: int) -> Iterator[tuple]:
    """
    Generator for creating hybridstrobes of any orders

    :param seq: a string with a nucleotide sequence
    :param k_size: length of each strobe
    :param w_min: minimum window offset to the previous window (wMin > 0)
    :param w_max: maximum window offset to the previous window (wMin <= wMax)
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :returns: an iterator for creating hybridstrobes
    """
    hash_list = [hash(seq[i:i+k_size]) for i in range(len(seq) - k_size + 1)]
    n_partition = 3
    w_p = (w_max - w_min) // n_partition

    tmp_index_dict = dict()
    for strobe_num in range(0, order-1):
        tmp_index = []
        for partition in range(0, n_partition):
            start = w_min + w_max*strobe_num + w_p*partition
            end = (
                w_max + w_max*strobe_num if partition + 1 == n_partition
                else w_min + w_max*strobe_num + w_p + w_p*partition
            )

            strobe_window = deque(hash_list[start: end])
            min_index, min_w = argmin(strobe_window)
            min_index = min_index + w_min + w_max*strobe_num + w_p*partition
            tmp_index.append(
                (
                    strobe_window,
                    min_w,
                    min_index,
                    start, end
                )
            )
        tmp_index_dict[strobe_num] = tmp_index

    for i in range(len(hash_list) - w_max*order):  # temporary iteration
        index_hash = hash_list[i]
        positions = [i, ]

        for strobe_num in range(0, order-1):
            tmp_index = []
            for window_numer, window in enumerate(tmp_index_dict[strobe_num]):
                # updating windows
                strobe_window, min_w, min_index, start, end = window
                new_w = hash_list[i + end]
                min_index, min_w = update_queue(
                    strobe_window, min_w, min_index, new_w, i, start, end
                )

                # update tmp_index_dict
                tmp_index_dict[strobe_num][window_numer] = (
                    strobe_window, min_w, min_index, start, end
                )
                tmp_index.append((min_index, min_w))

            next_i, next_index_hash = tmp_index[index_hash % n_partition]
            positions.append(next_i)
            index_hash = index_hash + (strobe_num+1) * (-1)**(strobe_num+1) * next_index_hash

        yield positions, index_hash


def hybridstrobes(seq: str, k_size: int, strobe_w_min_offset: int,
                  strobe_w_max_offset: int, w: int, order: int = 2) -> dict:
    """
    Hybrid between minstrobes and randstrobes that uses both independent minimizers and a conditional dependence between strobes

    :param seq: a string with a nucleotide sequence
    :param k_size: length of all strobes (len(strobe_1) +  ... + len(strobe_n))
    :param strobe_w_min_offset: minimum window offset to the previous window (wMin > 0)
    :param strobe_w_max_offset: maximum window offset to the previous window (wMin <= wMax)
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :returns: a dictionary with positions along the string as keys and the hybridstrobes as value
    """
    assert strobe_w_min_offset > 0, "Minimum strobemer offset has to be greater than 0 in this implementation"
    if k_size % order != 0:
        print("WARNING: kmer size is not evenly divisible with 2, will use {0} as kmer size: ".format(k_size - k_size % 2))
        k_size = k_size - k_size % order
    m_size = k_size//order

    if w == 1:
        hybridstrobes = {
            tuple(positions): index_hash
            for positions, index_hash in seq_to_hybridstrobes_iter(
                seq, m_size, strobe_w_min_offset, strobe_w_max_offset, w, order
            )
        }

    else:  # thin out hybridstrobes
        hybridstrobes_tmp = [
            (positions, index_hash)
            for positions, index_hash in seq_to_hybridstrobes_iter(
                seq, m_size, strobe_w_min_offset, strobe_w_max_offset, w, order
            )
        ]

        thinned_hybridstrobes = thinner(
            [index_hash for (positions, index_hash) in hybridstrobes_tmp],
            w
        )

        hybridstrobes = {}
        for p1, index_hash in thinned_hybridstrobes:
            if p1 < len(hybridstrobes_tmp):
                (positions, index_hash) = hybridstrobes_tmp[p1]
                hybridstrobes[tuple(positions)] = index_hash

    return hybridstrobes


def hybridstrobes_iter(seq: str, k_size: int, strobe_w_min_offset: int,
                       strobe_w_max_offset: int, w: int, order: int = 2,
                       buffer_size: int = 10000000) -> Iterator[tuple]:
    """
    Generator for creating hybridstrobes (less memory requiring)

    :param seq: a string with a nucleotide sequence
    :param k_size: length of all strobes (len(strobe_1) +  ... + len(strobe_n))
    :param strobe_w_min_offset: minimum window offset to the previous window (wMin > 0)
    :param strobe_w_max_offset: maximum window offset to the previous window (wMin <= wMax)
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :param buffer_size: size of buffer at the end of the sequence
    :returns: an iterator for creating hybridstrobes
    """
    for i in range(0, len(seq), buffer_size):
        substring = seq[i:i+buffer_size]
        # print(substring, len(substring))
        for p, m in hybridstrobes(substring, k_size, strobe_w_min_offset, strobe_w_max_offset, w, order).items():
            yield p, m


def seq_to_mixedhybridstrobes_iter(seq: str, k_size: int, w_min: int, w_max: int,
                                   w: int, order: int, denominator: int,
                                   numerator: int) -> Iterator[tuple]:
    """
    Generator for creating mixed hybridstrobes of any orders

    :param seq: a string with a nucleotide sequence
    :param k_size: length of each strobe
    :param w_min: minimum window offset to the previous window (wMin > 0)
    :param w_max: maximum window offset to the previous window (wMin <= wMax)
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :param denominator: denominator and numerator determine the fraction of sampled strobemers
    :param numerator: denominator and numerator determine the fraction of sampled strobemers
    :returns: an iterator for creating mixedhybridstrobes
    """
    hash_list = [hash(seq[i:i+k_size]) for i in range(len(seq) - k_size + 1)]
    n_partition = 3
    w_p = (w_max - w_min) // n_partition

    tmp_index_dict = dict()
    for strobe_num in range(0, order-1):
        tmp_index = []
        for partition in range(0, n_partition):
            start = w_min + w_max*strobe_num + w_p*partition
            end = (
                w_max + w_max*strobe_num if partition + 1 == n_partition
                else w_min + w_max*strobe_num + w_p + w_p*partition
            )

            strobe_window = deque(hash_list[start: end])
            min_index, min_w = argmin(strobe_window)
            min_index = min_index + w_min + w_max*strobe_num + w_p*partition
            tmp_index.append(
                (
                    strobe_window,
                    min_w,
                    min_index,
                    start, end
                )
            )
        tmp_index_dict[strobe_num] = tmp_index

    for i in range(len(hash_list) - w_max*order):  # temporary iteration
        index_hash = hash_list[i]
        positions = [i, ]

        for strobe_num in range(0, order-1):
            tmp_index = []
            for window_numer, window in enumerate(tmp_index_dict[strobe_num]):
                # updating windows
                strobe_window, min_w, min_index, start, end = window
                new_w = hash_list[i + end]
                min_index, min_w = update_queue(
                    strobe_window, min_w, min_index, new_w, i, start, end
                )

                # update tmp_index_dict
                tmp_index_dict[strobe_num][window_numer] = (
                    strobe_window, min_w, min_index, start, end
                )
                tmp_index.append((min_index, min_w))

            next_i, next_index_hash = tmp_index[index_hash % n_partition]
            positions.append(next_i)
            index_hash = index_hash + (strobe_num+1) * (-1)**(strobe_num+1) * next_index_hash

        # decide whether kmers should be sampled instead of mixedstrobes
        if hash_list[i] % denominator >= numerator:
            positions = [i + strobe * k_size for strobe in range(order)]
            index_hash = hash(seq[i:i+k_size * order])

        yield positions, index_hash


def mixedhybridstrobes(seq: str, k_size: int, strobe_w_min_offset: int,
                       strobe_w_max_offset: int, w: int,
                       order: int = 2, strobe_fraction: float = 0.5) -> dict:
    """
    Mixed protocol to produce specified hybridstrobes and k-mers content fractions

    :param seq: a string with a nucleotide sequence
    :param k_size: length of all strobes (len(strobe_1) +  ... + len(strobe_n))
    :param strobe_w_min_offset: minimum window offset to the previous window (wMin > 0)
    :param strobe_w_max_offset: maximum window offset to the previous window (wMin <= wMax)
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :param strobe_fraction: a fraction of sampled strobemers, rest kmers (0 <= strobe_fraction <= 1)
    :returns: a dictionary with positions along the string as keys and the mixed hybridstrobes/kmers as value
    """
    mixed_output = dict()
    fraction = Fraction(str(strobe_fraction))
    denominator = fraction.denominator
    numerator = fraction.numerator
    assert strobe_w_min_offset > 0, "Minimum strobemer offset has to be greater than 0 in this implementation"

    if k_size % order != 0:
        print("WARNING: kmer size is not evenly divisible with 2, will use {0} as kmer size: ".format(k_size - k_size % 2))
        k_size = k_size - k_size % order
    m_size = k_size//order

    if w == 1:
        mixed_output = {
            tuple(positions): index_hash
            for positions, index_hash in seq_to_mixedhybridstrobes_iter(
                seq, m_size, strobe_w_min_offset, strobe_w_max_offset, w, order,
                denominator, numerator
            )
        }

    else:  # thin out mixedhybridstrobes
        mixedhybridstrobes_tmp = [
            (positions, index_hash)
            for positions, index_hash in seq_to_mixedhybridstrobes_iter(
                seq, m_size, strobe_w_min_offset, strobe_w_max_offset, w, order,
                denominator, numerator
                )
            ]

        thinned_mixedhybridstrobes = thinner(
            [index_hash for (positions, index_hash) in mixedhybridstrobes_tmp],
            w
        )

        hybridstrobes = {}
        for p1, index_hash in thinned_mixedhybridstrobes:
            if p1 < len(mixedhybridstrobes_tmp):
                (positions, index_hash) = mixedhybridstrobes_tmp[p1]
                mixed_output[tuple(positions)] = index_hash

    return mixed_output


def mixedhybridstrobes_iter(seq: str, k_size: int, strobe_w_min_offset: int,
                         strobe_w_max_offset: int, w: int, order: int = 2,
                         strobe_fraction: float = 1,
                         buffer_size: int = 10000000) -> Iterator[tuple]:
    """
    Generator for creating mixedhybridstrobes (less memory requiring)

    :param seq: a string with a nucleotide sequence
    :param k_size: length of all strobes (len(strobe_1) +  ... + len(strobe_n))
    :param strobe_w_min_offset: minimum window offset to the previous window (wMin > 0)
    :param strobe_w_max_offset: maximum window offset to the previous window (wMin <= wMax)
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :param buffer_size: size of buffer at the end of the sequence
    :returns: an iterator for creating mixedhybridstrobes
    """
    for i in range(0, len(seq), buffer_size):
        substring = seq[i: i+buffer_size]
        for p, m in mixedhybridstrobes(substring, k_size, strobe_w_min_offset, strobe_w_max_offset, w, order, strobe_fraction).items():
            yield p, m
