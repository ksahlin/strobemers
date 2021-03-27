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

    elif order == 4:
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


# def minstrobes_iter(seq, k_size, order = 2, **kwargs):
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
#             yield minstrobe_order2(seq[p:min(p+w_1, len(seq))], m_size)

#     elif order == 3:
#         w_1 = kwargs["w_1"]
#         w_2 = kwargs["w_2"]  
#         if k_size % 3 != 0:
#             print("WARNING: kmer size is not evenly divisible with 3, will use {0} as kmer size: ".format(k_size - k_size % 3))
#             k_size = k_size - k_size % 3
#         m_size = k_size//3
#         for p in range(len(seq) - k_size +1):
#             yield minstrobe_order3(seq[p:min(p+m_size+w_1+w_2, len(seq))], m_size, min(w_1 + w_2, len(seq)-p - m_size)//2, min(w_1 + w_2, len(seq)-p - m_size)//2)

