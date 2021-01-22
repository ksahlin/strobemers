import operator

def argmin(values):
    min_index, min_value = min(enumerate(values), key=operator.itemgetter(1))
    return min_index, min_value


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


def kmers(seq, k_size):
    return set([seq[i:i+k_size] for i in range(len(seq) - k_size +1)])



def randomer_order2(subseq, m_size):
    k1 = subseq[0:m_size]
    f = lambda x: x
    mod = 2**26
    min_index, min_value = argmin([ hash(k1+ subseq[i:i+m_size]) for i in range(m_size, len(subseq) - m_size + 1)])
    min_k2 = subseq[m_size+ min_index:m_size+ min_index+m_size]
    # print(len(k1 + min_k2))
    return k1 + min_k2



def randomer_order3(subseq, m_size, N_1, N_2):
    # print(len(subseq),m_size, N_1, N_2, m_size, m_size+ N_1 - m_size + 1, [i for i in range(m_size + N_1, m_size + N_1 + N_2 - m_size + 1)])
    k1 = subseq[0:m_size]
    f = lambda x: x
    mod = 2**26
    min_index, min_value = argmin([ hash(k1+subseq[i:i+m_size]) for i in range(m_size, m_size+ N_1 - m_size + 1)])
    min_k2 = subseq[m_size + min_index: m_size+ min_index+m_size]

    min_index, min_value = argmin([ hash(k1 + min_k2 + subseq[i:i+m_size]) for i in range(m_size + N_1, m_size + N_1 + N_2 - m_size + 1)])
    min_k3 = subseq[m_size + N_1 + min_index: m_size + N_1+ min_index+m_size]

    return k1 + min_k2 + min_k3


def randomers(seq, k_size, order = 2, **kwargs):

    if order == 2:
        N_1 = kwargs["N_1"]
        if k_size % 2 != 0:
            print("WARNING: kmer size is not evenly divisible with 2, will use {0} as kmer size: ".format(k_size - k_size % 2))
            k_size = k_size - k_size % 2
        m_size = k_size//2
        randomers = {p : randomer_order2(seq[p:min(p+N_1, len(seq))], m_size) for p in range(len(seq) - k_size +1)}
        return randomers

    elif order == 3:
        N_1 = kwargs["N_1"]
        N_2 = kwargs["N_2"]  
        if k_size % 3 != 0:
            print("WARNING: kmer size is not evenly divisible with 3, will use {0} as kmer size: ".format(k_size - k_size % 3))
            k_size = k_size - k_size % 3
        m_size = k_size//3
        randomers = {p: randomer_order3(seq[p:min(p+m_size+N_1+N_2, len(seq))], m_size, min(N_1 + N_2, len(seq)-p - m_size)//2, min(N_1 + N_2, len(seq)-p - m_size)//2) for p in range(len(seq) - k_size +1)}
        return randomers





