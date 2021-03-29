import sys
import operator
import random 
import copy
from collections import defaultdict, deque

MAX = sys.maxsize

# def argmin(values):
#     min_index, min_value = min(enumerate(values), key=operator.itemgetter(1))
#     return min_index, min_value

def argmin(array):
    min_index = array.index(min(array))
    min_val =  array[min_index]
    return min_index, min_val

def thinner(hash_list, w):
    """
        Input: a list with hash values 
        Output: A list with tuples: (pos in original list, minimim hash value) for each window of w hashes
    """
    window_hashes = deque(hash_list[:w])
    min_index, curr_min_hash = argmin(window_hashes)
    thinned_hash_list = [ (min_index, curr_min_hash) ]

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
            thinned_hash_list.append( (min_index + i + 1 - w, curr_min_hash) )

        # Previous minimizer still in window, we only need to compare with the recently added kmer 
        elif new_hash < curr_min_hash:
            curr_min_hash = new_hash
            thinned_hash_list.append( (i, curr_min_hash) )


    return thinned_hash_list


def minimizers(seq, k_size, w_size):
    # kmers = [seq[i:i+k_size] for i in range(len(seq)-k_size) ]
    w = w_size - k_size
    window_kmers = deque([seq[i:i+k_size] for i in range(w +1)])
    curr_min = min(window_kmers)
    minimizers = [ (curr_min, list(window_kmers).index(curr_min)) ]

    for i in range(w+1,len(seq) - k_size):
        new_kmer = seq[i:i+k_size]
        # updateing window
        discarded_kmer = window_kmers.popleft()
        window_kmers.append(new_kmer)

        # we have discarded previous windows minimizer, look for new minimizer brute force
        if curr_min == discarded_kmer: 
            curr_min = min(window_kmers)
            minimizers.append( (curr_min, list(window_kmers).index(curr_min) + i - w ) )

        # Previous minimizer still in window, we only need to compare with the recently added kmer 
        elif new_kmer < curr_min:
            curr_min = new_kmer
            minimizers.append( (curr_min, i) )

    return minimizers


def spaced_kmers(seq, k_size, span_size, positions, w):
    '''
        Parameters:
        Positions : A set of positions to consider for the spaced k-mer

        Returns: a dictionary with positions along the string as keys and the spaced kmer as value 
    '''
    assert len(positions) == k_size

    if w > 1:
        hash_seq_list = [(i, hash( "".join([seq[i + j] for j in range(span_size) if j in positions ]))) for i in range(len(seq) - span_size +1)]
        hash_seq_list_thinned = thinner([h for i,h in hash_seq_list], w) # produce a subset of positions, still with samme index as in full sequence
        spaced_kmers_pos = {i : h for i, h in hash_seq_list_thinned }

    else:
        spaced_kmers_pos = {i :  hash("".join([seq[i + j] for j in range(span_size) if j in positions ])) for i in range(len(seq) - span_size +1)}
    # print(positions, len(positions), span_size)
    # well, this is not the most time efficient way to sample spaced kmers but works for simulations...
    return spaced_kmers_pos


def spaced_kmers_iter(seq, k_size, span_size, positions):
    '''
        Parameters:
        Positions : A set of positions to consider for the spaced k-mer

        Returns: a dictionary with positions along the string as keys and the spaced kmer as value 
    '''
    assert len(positions) == k_size
    # print(positions, len(positions), span_size)
    # well, this is not the most time efficient way to sample spaced kmers but works for simulations...
    for i in range(len(seq) - span_size +1):
        yield  hash("".join([seq[i + j] for j in range(span_size) if j in positions ]))



def kmers(seq, k_size, w):
    if w > 1:
        hash_seq_list = [(i, hash(seq[i:i+k_size])) for i in range(len(seq) - k_size +1)]
        hash_seq_list_thinned = thinner([h for i,h in hash_seq_list], w) # produce a subset of positions, still with samme index as in full sequence
        kmers_pos = {i : h for i, h in hash_seq_list_thinned }

    else:
        kmers_pos = {i : hash(seq[i:i+k_size]) for i in range(len(seq) - k_size +1)}
    
    return kmers_pos


def kmer_iter(seq, k_size, w):
    if w > 1:
        hash_seq_list = [(i, hash(seq[i:i+k_size])) for i in range(len(seq) - k_size +1)]
        hash_seq_list_thinned = thinner([h for i,h in hash_seq_list], w) # produce a subset of positions, still with samme index as in full sequence
        for p, h in hash_seq_list_thinned:
            yield p, h
    else:
        hash_seq_list = [(i, hash(seq[i:i+k_size])) for i in range(len(seq) - k_size +1)]
        for p, h in hash_seq_list:   
            yield p, h


def randstrobe_order2(hash_seq_list, start, stop, hash_m1, prime):
    min_index, min_value = argmin([ (hash_m1 + hash_seq_list[i][1]) % prime for i in range(start, stop)])
    min_hash_val = hash_m1 - hash_seq_list[start + min_index][1]
    return min_index, min_hash_val


def seq_to_randstrobes2_iter(seq, k_size, strobe_w_min_offset, strobe_w_max_offset, prime, w):
    hash_seq_list = [(i, hash(seq[i:i+k_size])) for i in range(len(seq) - k_size +1)]
    if w > 1:
        hash_seq_list_thinned = thinner([h for i,h in hash_seq_list], w) # produce a subset of positions, still with samme index as in full sequence
    else:
        hash_seq_list_thinned = hash_seq_list
    
    # assert len(hash_seq_list[:-k_size]) == len(hash_seq_list) - k_size

    for (p, hash_m1) in hash_seq_list_thinned: #[:-k_size]:
        if p >= len(hash_seq_list) - k_size:
            break
        # hash_m1 = hash_seq_list[p]
        window_p_start = p + strobe_w_min_offset if p + strobe_w_max_offset <= len(hash_seq_list) else max( (p + strobe_w_min_offset) -  (p + strobe_w_max_offset - len(hash_seq_list)), p )
        window_p_end = min(p + strobe_w_max_offset, len(hash_seq_list))
        # print(window_p_start, window_p_end)
        min_index, hash_value = randstrobe_order2(hash_seq_list, window_p_start, window_p_end, hash_m1, prime)
        p2 = window_p_start + min_index
        yield p, p2, hash_value



def randstrobe_order3(hash_seq_list, start1, stop1, start2, stop2, hash_m1, prime):
    min_index1, min_value = argmin([ (hash_m1 + hash_seq_list[i][1]) % prime for i in range(start1, stop1)])
    min_hash_val = hash_m1 - hash_seq_list[start1 + min_index1][1]

    min_index2, min_value = argmin([ (min_hash_val + hash_seq_list[i][1]) % prime for i in range(start2, stop2)])
    min_hash_val = min_hash_val + 2*hash_seq_list[start2 + min_index2][1]

    return min_index1, min_index2, min_hash_val

def seq_to_randstrobes3_iter(seq, k_size, strobe_w_min_offset, strobe_w_max_offset, prime, w):
    hash_seq_list = [(i, hash(seq[i:i+k_size])) for i in range(len(seq) - k_size +1)]
    if w > 1:
        hash_seq_list_thinned = thinner([h for i,h in hash_seq_list], w) # produce a subset of positions, still with samme index as in full sequence
    else:
        hash_seq_list_thinned = hash_seq_list
    
    # assert len(hash_seq_list[:-k_size]) == len(hash_seq_list) - k_size

    for (p, hash_m1) in hash_seq_list_thinned: #[:-k_size]:
        if p >= len(hash_seq_list) - 2*k_size:
            break
        # hash_m1 = hash_seq_list[p]
        window_p2_start = p + strobe_w_max_offset + strobe_w_min_offset if p + 2*strobe_w_max_offset <= len(hash_seq_list) else max( (p + strobe_w_max_offset + strobe_w_min_offset) -  (p+2*strobe_w_max_offset - len(hash_seq_list)), p )
        window_p2_end = min(p + 2*strobe_w_max_offset, len(hash_seq_list))

        window_p1_start = p + strobe_w_min_offset if p + 2*strobe_w_max_offset <= len(hash_seq_list) else max(p,  len(hash_seq_list)  + 2*(strobe_w_min_offset - strobe_w_max_offset))
        window_p1_end = min(p + strobe_w_max_offset, len(hash_seq_list))
        # print(window_p1_start, window_p1_end,  window_p2_start, window_p2_end, len(seq))
        # assert window_p1_start < window_p1_end
        # print(window_p1_start, window_p1_end)
        min_index_s1, min_index_s2, hash_value = randstrobe_order3(hash_seq_list, window_p1_start, window_p1_end, window_p2_start, window_p2_end, hash_m1, prime)
        p2 = window_p1_start + min_index_s1
        p3 = window_p2_start + min_index_s2
        yield p, p2, p3, hash_value



# def randstrobe_order3(subseq, m_size, w_1, w_2):
#     # print(len(subseq),m_size, w_1, w_2, m_size, m_size+ w_1 - m_size + 1, [i for i in range(m_size + w_1, m_size + w_1 + w_2 - m_size + 1)])
#     k1 = subseq[0:m_size]
#     min_index, min_value = argmin([ hash(k1+subseq[i:i+m_size]) for i in range(m_size, m_size+ w_1 - m_size + 1)])
#     min_k2 = subseq[m_size + min_index: m_size+ min_index+m_size]

#     min_index, min_value = argmin([ hash(k1 + min_k2 + subseq[i:i+m_size]) for i in range(m_size + w_1, m_size + w_1 + w_2 - m_size + 1)])
#     min_k3 = subseq[m_size + w_1 + min_index: m_size + w_1+ min_index+m_size]

#     return k1 + min_k2 + min_k3


# def randstrobe_order4(subseq, m_size, w_1, w_2, w_3):
#     k1 = subseq[0:m_size]
#     min_index, min_value = argmin([ hash(k1+subseq[i:i+m_size]) for i in range(m_size, m_size+ w_1 - m_size + 1)])
#     min_k2 = subseq[m_size + min_index: m_size+ min_index+m_size]

#     min_index, min_value = argmin([ hash(k1 + min_k2 + subseq[i:i+m_size]) for i in range(m_size + w_1, m_size + w_1 + w_2 - m_size + 1)])
#     min_k3 = subseq[m_size + w_1 + min_index: m_size + w_1+ min_index+m_size]

#     min_index, min_value = argmin([ hash(k1 + min_k2 + min_k3 + subseq[i:i+m_size]) for i in range(m_size + w_1 + w_2, m_size + w_1 + w_2 + w_3 - m_size + 1)])
#     min_k4 = subseq[m_size + w_1 + min_index: m_size + w_1+ min_index+m_size]

#     return k1 + min_k2 + min_k3 + min_k4


def randstrobes(seq, k_size, strobe_w_min_offset, strobe_w_max_offset, w, order = 2):
    prime = 997
    assert strobe_w_min_offset > 0, "Minimum strobemer offset has to be greater than 0 in this implementation"
    if order == 2:
        if k_size % 2 != 0:
            print("WARNING: kmer size is not evenly divisible with 2, will use {0} as kmer size: ".format(k_size - k_size % 2))
            k_size = k_size - k_size % 2
        m_size = k_size//2
        randstrobes = {(p1,p2): h for p1,p2,h in seq_to_randstrobes2_iter(seq, m_size, strobe_w_min_offset, strobe_w_max_offset, prime, w)}
        return randstrobes

    elif order == 3:
        if k_size % 3 != 0:
            print("WARNING: kmer size is not evenly divisible with 3, will use {0} as kmer size: ".format(k_size - k_size % 3))
            k_size = k_size - k_size % 3
        m_size = k_size//3
        randstrobes = {(p1,p2,p3): h for p1,p2,p3,h in seq_to_randstrobes3_iter(seq, m_size, strobe_w_min_offset, strobe_w_max_offset, prime, w)}
        return randstrobes

    elif order == 4:
        raise NotImplementedError
        # w_1 = kwargs["w_1"]
        # w_2 = kwargs["w_2"]  
        # w_3 = kwargs["w_3"]  
        # if k_size % 4 != 0:
        #     print("WARNING: kmer size is not evenly divisible with 4, will use {0} as kmer size: ".format(k_size - k_size % 4))
        #     k_size = k_size - k_size % 4
        # m_size = k_size//4
        # randstrobes = {p: randstrobe_order4(seq[p:min(p+m_size+w_1+w_2+w_3, len(seq))], m_size, 
        #                                         w_1 if w_1 + w_2 + w_3 < len(seq) - p - m_size else (w_1 - ( w_1 + w_2 + w_3 - len(seq)-p - m_size)//3), 
        #                                         w_2 if w_1 + w_2 + w_3 < len(seq) - p - m_size else (w_2 - ( w_1 + w_2 + w_3 - len(seq)-p - m_size)//3),
        #                                         w_3 if w_1 + w_2 + w_3 < len(seq) - p - m_size else (w_3 - ( w_1 + w_2 + w_3 - len(seq)-p - m_size)//3))
        #                                         for p in range(len(seq) - k_size +1)}
        # return randstrobes


# def randstrobes_iter(seq, k_size, order = 2, **kwargs):
#     """
#         Low memory consumption due to not precalculating hash values, but more time consuming than randstrobes function.
#         TODO: Implement buffer window that precalculates hash values of a subsecuence window that we use to calculate randstrobes
#         from. This will be as fast as randstrobes without much memory overhead.
#     """

#     if order == 2:
#         w_1 = kwargs["w_1"]
#         if k_size % 2 != 0:
#             print("WARNING: kmer size is not evenly divisible with 2, will use {0} as kmer size: ".format(k_size - k_size % 2))
#             k_size = k_size - k_size % 2
#         m_size = k_size//2
#         for p in range(len(seq) - k_size +1):
#             yield randstrobe_order2(seq[p:min(p+m_size+w_1, len(seq))], m_size)

#     elif order == 3:
#         w_1 = kwargs["w_1"]
#         w_2 = kwargs["w_2"]  
#         if k_size % 3 != 0:
#             print("WARNING: kmer size is not evenly divisible with 3, will use {0} as kmer size: ".format(k_size - k_size % 3))
#             k_size = k_size - k_size % 3
#         m_size = k_size//3
#         for p in range(len(seq) - k_size +1):
#             yield randstrobe_order3(seq[p:min(p+m_size+w_1+w_2, len(seq))], m_size, 
#                                     w_1 if w_1 + w_2 < len(seq) - p - m_size else (w_1 - ( w_1 + w_2 - len(seq)-p - m_size)//2), 
#                                     w_2 if w_1 + w_2 < len(seq) - p - m_size else (w_2 - ( w_1 + w_2 - len(seq)-p - m_size)//2))

#     elif order == 4:
#         w_1 = kwargs["w_1"]
#         w_2 = kwargs["w_2"]  
#         w_3 = kwargs["w_3"]  
#         if k_size % 4 != 0:
#             print("WARNING: kmer size is not evenly divisible with 4, will use {0} as kmer size: ".format(k_size - k_size % 4))
#             k_size = k_size - k_size % 4
#         m_size = k_size//4
#         for p in range(len(seq) - k_size +1):
#             yield randstrobe_order4(seq[p:min(p+m_size+w_1+w_2+w_3, len(seq))], m_size, 
#                                     w_1 if w_1 + w_2 + w_3 < len(seq) - p - m_size else (w_1 - ( w_1 + w_2 + w_3 - len(seq)-p - m_size)//3), 
#                                     w_2 if w_1 + w_2 + w_3 < len(seq) - p - m_size else (w_2 - ( w_1 + w_2 + w_3 - len(seq)-p - m_size)//3),
#                                     w_3 if w_1 + w_2 + w_3 < len(seq) - p - m_size else (w_3 - ( w_1 + w_2 + w_3 - len(seq)-p - m_size)//3))



# def minstrobe_order2(subseq, m_size):
#     k1 = subseq[0:m_size]
#     min_index, min_value = argmin([ hash(k1) - hash(subseq[i:i+m_size]) for i in range(m_size, len(subseq) - m_size + 1)])
#     min_k2 = subseq[m_size+ min_index:m_size+ min_index+m_size]
#     # print(len(k1 + min_k2))
#     return k1 + min_k2



# def minstrobe_order3(subseq, m_size, w_1, w_2):
#     # print(len(subseq),m_size, w_1, w_2, m_size, m_size+ w_1 - m_size + 1, [i for i in range(m_size + w_1, m_size + w_1 + w_2 - m_size + 1)])
#     k1 = subseq[0:m_size]
#     min_index, min_value = argmin([ hash(k1) - hash(subseq[i:i+m_size]) for i in range(m_size, m_size+ w_1 - m_size + 1)])
#     min_k2 = subseq[m_size + min_index: m_size+ min_index+m_size]

#     min_index, min_value = argmin([ (hash(k1) - hash(min_k2))/ hash(subseq[i:i+m_size]) for i in range(m_size + w_1, m_size + w_1 + w_2 - m_size + 1)])
#     min_k3 = subseq[m_size + w_1 + min_index: m_size + w_1+ min_index+m_size]

#     return k1 + min_k2 + min_k3


# def minstrobes(seq, k_size, order = 2, **kwargs):

#     if order == 2:
#         w_1 = kwargs["w_1"]
#         if k_size % 2 != 0:
#             print("WARNING: kmer size is not evenly divisible with 2, will use {0} as kmer size: ".format(k_size - k_size % 2))
#             k_size = k_size - k_size % 2
#         m_size = k_size//2
#         minstrobes = {p : minstrobe_order2(seq[p:min(p+w_1, len(seq))], m_size) for p in range(len(seq) - k_size +1)}
#         return minstrobes

#     elif order == 3:
#         w_1 = kwargs["w_1"]
#         w_2 = kwargs["w_2"]  
#         if k_size % 3 != 0:
#             print("WARNING: kmer size is not evenly divisible with 3, will use {0} as kmer size: ".format(k_size - k_size % 3))
#             k_size = k_size - k_size % 3
#         m_size = k_size//3
#         minstrobes = {p: minstrobe_order3(seq[p:min(p+m_size+w_1+w_2, len(seq))], m_size, min(w_1 + w_2, len(seq)-p - m_size)//2, min(w_1 + w_2, len(seq)-p - m_size)//2) for p in range(len(seq) - k_size +1)}
#         return minstrobes


def seq_to_minstrobes2_iter(seq, k_size, strobe_w_min_offset, strobe_w_max_offset, prime, w):
    hash_seq_list = [(i, hash(seq[i:i+k_size])) for i in range(len(seq) - k_size +1)]
    strobes = deque(thinner([h for i,h in hash_seq_list], strobe_w_max_offset - strobe_w_min_offset)) # produce a subset of positions, still with samme index as in full sequence

    if w > 1:
        hash_seq_list_thinned = thinner([h for i,h in hash_seq_list], w) # produce a subset of positions, still with same index as in full sequence
    else:
        hash_seq_list_thinned = hash_seq_list
    
    # assert len(hash_seq_list[:-k_size]) == len(hash_seq_list) - k_size

    for (p, hash_m1) in hash_seq_list_thinned: #[:-k_size]:
        if p >= len(hash_seq_list) - k_size:
            break
        # print(p,len(hash_seq_list) - k_size, len(hash_seq_list), len(seq), len(strobes), strobes)
        if p + k_size + strobe_w_min_offset < len(seq):
            while strobes[0][0] < min(p + k_size + strobe_w_min_offset, len(hash_seq_list)-1):
                l = strobes.popleft()

        # print(p, len(hash_seq_list) - k_size, len(hash_seq_list), len(seq), len(strobes), strobes[0])
        p2, hash_val = strobes[0]
        hash_value = hash_m1 - hash_val
        yield p, p2, hash_value



def seq_to_minstrobes3_iter(seq, k_size, strobe_w_min_offset, strobe_w_max_offset, prime, w):
    hash_seq_list = [(i, hash(seq[i:i+k_size])) for i in range(len(seq) - k_size +1)]
    strobes = deque(thinner([h for i,h in hash_seq_list], strobe_w_max_offset - strobe_w_min_offset)) # produce a subset of positions, still with samme index as in full sequence
    strobes2 = copy.deepcopy(strobes)
    if w > 1:
        hash_seq_list_thinned = thinner([h for i,h in hash_seq_list], w) # produce a subset of positions, still with same index as in full sequence
    else:
        hash_seq_list_thinned = hash_seq_list
    
    # assert len(hash_seq_list[:-k_size]) == len(hash_seq_list) - k_size

    for (p, hash_m1) in hash_seq_list_thinned: #[:-k_size]:
        if p >= len(hash_seq_list) - 2*k_size:
            break
        # print(p,len(hash_seq_list) - k_size, len(hash_seq_list), len(seq), len(strobes), strobes)

        if p + strobe_w_max_offset + strobe_w_min_offset < len(seq):
            while strobes2[0][0] <  min(p + strobe_w_min_offset + strobe_w_max_offset, len(hash_seq_list)-1):
                l = strobes2.popleft()

        if p + k_size + strobe_w_min_offset < len(seq):
            while strobes[0][0] <  min(p + k_size + strobe_w_min_offset, len(hash_seq_list)-1):
                l = strobes.popleft()
        # print(p, len(hash_seq_list) - k_size, len(hash_seq_list), len(seq), len(strobes), strobes[0], len(strobes2))

        p2, hash_val2 = strobes[0]
        p3, hash_val3 = strobes2[0]
        hash_value = hash_m1 - hash_val2 + 2*hash_val3
        yield p, p2, p3, hash_value




def minstrobes(seq, k_size, strobe_w_min_offset, strobe_w_max_offset, w, order = 2):
    prime = 997
    assert strobe_w_min_offset > 0, "Minimum strobemer offset has to be greater than 0 in this implementation"
    if order == 2:
        if k_size % 2 != 0:
            print("WARNING: kmer size is not evenly divisible with 2, will use {0} as kmer size: ".format(k_size - k_size % 2))
            k_size = k_size - k_size % 2
        m_size = k_size//2
        minstrobes = {(p1,p2): h for p1,p2,h in seq_to_minstrobes2_iter(seq, m_size, strobe_w_min_offset, strobe_w_max_offset, prime, w)}
        return minstrobes

    elif order == 3:
        if k_size % 3 != 0:
            print("WARNING: kmer size is not evenly divisible with 3, will use {0} as kmer size: ".format(k_size - k_size % 3))
            k_size = k_size - k_size % 3
        m_size = k_size//3
        minstrobes = {(p1,p2,p3): h for p1,p2,p3,h in seq_to_minstrobes3_iter(seq, m_size, strobe_w_min_offset, strobe_w_max_offset, prime, w)}
        return minstrobes

    else:
        raise NotImplementedError



def minstrobes_iter(seq, k_size, strobe_w_min_offset, strobe_w_max_offset, w, order = 2, buffer_size = 10000000):
    
    for i in range(0, len(seq), buffer_size):
        substring = seq[i:i+buffer_size] 
        # print(substring, len(substring))
        for p, m in minstrobes(substring, k_size, strobe_w_min_offset, strobe_w_max_offset, w, order = order).items():
            yield m


def randstrobes_iter(seq, k_size, strobe_w_min_offset, strobe_w_max_offset, w, order = 2, buffer_size = 10000000):
    
    for i in range(0, len(seq), buffer_size):
        substring = seq[i:i+buffer_size] 
        for p, m in randstrobes(substring, k_size, strobe_w_min_offset, strobe_w_max_offset, w, order = order).items():
            yield m


def update_queue(q, curr_min, min_index, new_hash, i, start_offset, end_offset):
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



def seq_to_hybridstrobes2_iter(seq, k_size, w_min, w_max, w):
    hash_list = [ hash(seq[i:i+k_size]) for i in range(len(seq) - k_size +1)]
    n_partition = 2
    w_p = (w_max - w_min ) // n_partition

    win1 = deque(hash_list[w_min : w_min + w_p])
    min_index1, min_w1 = argmin(win1)
    min_index1 = min_index1 + w_min

    win2 = deque(hash_list[w_min+w_p : w_min + 2*w_p])
    min_index2, min_w2 = argmin(win2)
    min_index2 = min_index2 + w_min + w_p 

    # win3 = deque(hash_list[w_min+2*w_p : w_min + 3*w_p])
    # min_index3, min_w3 = argmin(win3)
    # min_index3 = min_index3 + w_min+2*w_p

    # win4 = deque(hash_list[w_min+3*w_p : w_min + 4*w_p])
    # min_index4, min_w4 = argmin(win4)
    # min_index4 = min_index4 + w_min+2*w_p

    for i in range(len(hash_list) - w_min - n_partition*w_p): # temporary iteration
        m1 = hash_list[i]

        # updating windows
        new_w1 = hash_list[i + w_min + w_p]
        min_index1, min_w1 = update_queue(win1, min_w1, min_index1, new_w1, i, w_min, w_min + w_p)
        # print(len(win1), win1)
        new_w2 = hash_list[i + w_min + 2*w_p]
        min_index2, min_w2 = update_queue(win2, min_w2, min_index2, new_w2, i, w_min + w_p,  w_min + 2*w_p)

        # new_w3 = hash_list[i+ w_min + 3*w_p]
        # min_index3, min_w3 = update_queue(win3, min_w3, min_index3, new_w3, i, w_min + 2*w_p,  w_min + 3*w_p)

        # new_w4 = hash_list[i+ w_min + 4*w_p]
        # min_index4, min_w4 = update_queue(win4, min_w4, min_index4, new_w4, i, w_min + 3*w_p,  w_min + 4*w_p)

        # print(i, min_index1, min_w1, min_w2, min_w3)
        r =  m1 % n_partition
        if r == 0:
            # print(i, 1,m1 - min_w1)
            yield i, min_index1, m1 - min_w1
        elif r == 1:
            # print(i, 2,m1 - min_w2)
            yield i, min_index2, m1 - min_w2
        # elif r == 2:
        #     # print(i, 3, m1 - min_w3)
        #     yield i, min_index3, m1 - min_w3
        # else:
        #     # print(i, 3, m1 - min_w3)
        #     yield i, min_index4, m1 - min_w4


def seq_to_hybridstrobes3_iter(seq, k_size, w_min, w_max, w):
    hash_list = [ hash(seq[i:i+k_size]) for i in range(len(seq) - k_size +1)]
    n_partition = 2
    w_p = (w_max - w_min ) // n_partition

    s1_win1 = deque(hash_list[w_min : w_min + w_p])
    min_index1, min_w1 = argmin(s1_win1)
    min_index1 = min_index1 + w_min

    s1_win2 = deque(hash_list[w_min+w_p : w_max])
    min_index2, min_w2 = argmin(s1_win2)
    min_index2 = min_index2 + w_min + w_p 

    s2_win1 = deque(hash_list[w_min + w_max : w_min + w_max + w_p])
    min_index3, min_w3 = argmin(s2_win1)
    min_index3 = min_index3 + w_min + w_max

    s2_win4 = deque(hash_list[w_min + w_max + w_p: 2*w_max])
    min_index4, min_w4 = argmin(s2_win4)
    min_index4 = min_index4 + w_min + w_max + w_p

    for i in range(len(hash_list) - n_partition*w_max): # temporary iteration
        m1 = hash_list[i]

        # updating windows
        new_w1 = hash_list[i + w_min + w_p]
        min_index1, min_w1 = update_queue(s1_win1, min_w1, min_index1, new_w1, i, w_min, w_min + w_p)
        # print(len(win1), win1)
        new_w2 = hash_list[i + w_max]
        min_index2, min_w2 = update_queue(s1_win2, min_w2, min_index2, new_w2, i, w_min + w_p,  w_max)

        new_w3 = hash_list[i+ w_min + w_max + w_p]
        min_index3, min_w3 = update_queue(s2_win1, min_w3, min_index3, new_w3, i, w_min + w_max, w_min + w_max + w_p )

        new_w4 = hash_list[i+ 2*w_max]
        min_index4, min_w4 = update_queue(s2_win4, min_w4, min_index4, new_w4, i, w_min + w_max + w_p, 2*w_max)

        # print(i, min_index1, min_w1, min_w2, min_w3)
        r =  m1 % n_partition
        if r == 0:
            # print(i, 1,m1 - min_w1)
            i2 = min_index1
            m2 = min_w1
            # yield i, min_index1, m1 - min_w1
        elif r == 1:
            # print(i, 2,m1 - min_w2)
            i2 = min_index2
            m2 = min_w2
            # yield i, min_index2, m1 - min_w2
        
        r2 =  (m1 - m2) % n_partition
        if r2 == 0:
            # print(i, 3, m1 - min_w3)
            # print(i, m1, m2, min_w3, m1 - m2 + 2*min_w3)
            yield i, i2, min_index3, m1 - m2 + 2*min_w3
        elif r2 == 1:
            # print(i, m1, m2, min_w4, m1 - m2 + 2*min_w4)
            yield i, i2, min_index4, m1 - m2 + 2*min_w4


def hybridstrobes(seq, k_size, strobe_w_min_offset, strobe_w_max_offset, w, order = 2):
    assert strobe_w_min_offset > 0, "Minimum strobemer offset has to be greater than 0 in this implementation"
    if order == 2:
        if k_size % 2 != 0:
            print("WARNING: kmer size is not evenly divisible with 2, will use {0} as kmer size: ".format(k_size - k_size % 2))
            k_size = k_size - k_size % 2
        m_size = k_size//2
        if w == 1:
            hybridstrobes = {(p1,p2): h for p1,p2,h in seq_to_hybridstrobes2_iter(seq, m_size, strobe_w_min_offset, strobe_w_max_offset, w)}
            # print(len(hybridstrobes))

        else: # thin out hybridstrobes
            hybridstrobes_tmp = [ (p1, p2, h) for p1,p2,h in seq_to_hybridstrobes2_iter(seq, m_size, strobe_w_min_offset, strobe_w_max_offset, w)]
            thinned_hybridstrobes = thinner([h for (p1, p2,h) in hybridstrobes_tmp], w)
            hybridstrobes = {}
            for p1, h in thinned_hybridstrobes:
                if p1 < len(hybridstrobes_tmp):
                    (p1, p2, h) = hybridstrobes_tmp[p1]
                    hybridstrobes[(p1,p2)] = h
        return hybridstrobes

    elif order == 3:
        if k_size % 3 != 0:
            print("WARNING: kmer size is not evenly divisible with 3, will use {0} as kmer size: ".format(k_size - k_size % 3))
            k_size = k_size - k_size % 3
        m_size = k_size//3
        if w == 1:
            hybridstrobes = {(p1,p2, p3): h for p1,p2,p3,h in seq_to_hybridstrobes3_iter(seq, m_size, strobe_w_min_offset, strobe_w_max_offset, w)}
            # print(len(hybridstrobes))
        else: # thin out hybridstrobes
            hybridstrobes_tmp = [ (p1, p2, p3, h) for p1,p2,p3,h in seq_to_hybridstrobes3_iter(seq, m_size, strobe_w_min_offset, strobe_w_max_offset, w)]
            thinned_hybridstrobes = thinner([h for (p1, p2, p3, h) in hybridstrobes_tmp], w)
            hybridstrobes = {}
            for p1, h in thinned_hybridstrobes:
                if p1 < len(hybridstrobes_tmp):
                    (p1, p2, p3, h) = hybridstrobes_tmp[p1]
                    hybridstrobes[(p1,p2, p3)] = h
        return hybridstrobes

    else:
        raise NotImplementedError

