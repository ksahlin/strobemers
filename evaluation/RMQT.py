
from time import time
import queue 
import random
import copy




class Node:
    __slots__ = ['d', 'j', 'Cj', "j_max"]
    def __init__(self, d, j, Cj, j_max):
        self.d = d
        self.j = j
        self.Cj = Cj
        self.j_max = j_max


def construct_tree(tree, leafs, n):   
    # assign values to leaves of  
    # the segment tree  
    for i in range(n):  
        tree[n + i] = leafs[i] 
      
    # assign values to remaining nodes  
    # to compute maximum in a given range  
    for i in range(n - 1, 0, -1):  
        # print(i)
        # print(tree)
        max_coord_node = max([tree[2 * i ], tree[2 * i + 1]], key = lambda x: x.d ) # should always be left one though
        tree_node = Node(max_coord_node.d, max_coord_node.j, max_coord_node.Cj, max_coord_node.j_max)
        tree[i] = tree_node


def range_query(tree, l, r, n):
    # print(l, r)
    assert l <= r
    # search for both at the same time 
    l_pos = 1 # root position
    r_pos = 1 # root position
    V_biss = set()
    V_prime = set()
    v = 1
    while True: # not yet reached a leaf

        l_pos_l_child = 2*l_pos
        l_pos_r_child = 2*l_pos + 1
        
        r_pos_l_child = 2*r_pos
        r_pos_r_child = 2*r_pos + 1

        ## traversing downwards for right path
        if r >= tree[r_pos_l_child].d:
            r_pos = r_pos_r_child
        else: # tree[left_child_index].d <= c:
            r_pos = r_pos_l_child         

        ## traversing downwards for left path
        if  l > tree[l_pos_l_child].d:
            l_pos = l_pos_r_child 
        else:
            l_pos = l_pos_l_child 

        if l_pos == r_pos:
            v = l_pos

        ## Adding nodes to sets V' and V'' if the two paths separates into two (See page 20-21 in GSAD book, Makinen et al)
        if l_pos != r_pos and v != l_pos >> 1: # paths have split and v is not the immediate predecendant
            if r_pos == r_pos_r_child: # we took a right step for right path
                V_biss.add(r_pos_l_child)

            if l_pos == l_pos_l_child: # we took a left step for left path
                V_prime.add(l_pos_r_child)
        
        # exit if we are at leafs
        if l_pos >= n and r_pos >= n:
            break

    R = tree[r_pos].d
    if R == r:
        V_biss.add(r_pos)
    L = tree[l_pos].d
    if L >= l:
        V_prime.add(l_pos)

    # #####
    # tmp = sorted(V_prime, key = lambda x: tree[x].Cj)
    # if len(tmp) > 1:
    #     if tree[ tmp[-1] ].Cj  == tree[ tmp[-2] ].Cj:
    #         chosen = max(V_prime, key = lambda x: tree[x].Cj)
    #         all_ = [zz for zz in V_prime if tree[ chosen ].Cj == tree[ zz ].Cj]
    #         print("LOOOOOL", tree[ chosen ].j_max, [tree[ zz ].j_max for zz in all_])
    
    # tmp = sorted(V_biss, key = lambda x: tree[x].Cj)
    # if len(tmp) > 1:
    #     if tree[ tmp[-1] ].Cj  == tree[ tmp[-2] ].Cj:
    #         chosen = max(V_biss, key = lambda x: tree[x].Cj)
    #         all_ = [zz for zz in V_biss if tree[ chosen ].Cj == tree[ zz ].Cj]
    #         print("LOOOOOL", tree[ chosen ].j_max, [tree[ zz ].j_max for zz in all_])
    # ########


    if V_prime:
        # vl = max(V_prime, key = lambda x: tree[x].Cj)
        # The sorted function guarantees largest j_max is returned to be consistent (i.e., closest mem)
        vl = max(sorted(V_prime, key = lambda x: tree[x].j_max, reverse = True), key = lambda x: tree[x].Cj)
    else:
        vl = None

    if V_biss:
        # vr = max(V_biss, key = lambda x: tree[x].Cj)
        # The sorted function guarantees largest j_max is returned to be consistent (i.e., closest mem)
        vr = max(sorted(V_biss, key = lambda x: tree[x].j_max, reverse = True), key = lambda x: tree[x].Cj)

    else:
        vr = None
       
    if vl is not None and vr is not None:
        v_max_pos = max(sorted([vl,vr], key = lambda x: tree[x].j_max, reverse = True), key = lambda x: tree[x].Cj)
        # v_max_pos = max([vl,vr], key = lambda x: tree[x].Cj)
        C_max = tree[v_max_pos].Cj
    elif vl is not None:
        v_max_pos = vl
        C_max = tree[vl].Cj
    elif vr is not None:
        v_max_pos = vr
        C_max = tree[vr].Cj
    else:
        print("BUG", l,r)
        print(tree)
        sys.exit()
    
    # print(l,r , "v", v, "RQ2",V_prime, V_biss, "RQ2",C_max)
    return C_max, tree[v_max_pos].j_max, v_max_pos



def update(tree, leaf_pos, value, n): 
    # change the index to leaf node first  
    pos = leaf_pos + n
    # update the value at the leaf node  
    # at the exact index 
    tree[pos].Cj = value  
    while (pos > 1): 
          
        # move up one level at a time in the tree  
        pos >>= 1  
        # print('updating: ', pos)
        # update the values in the nodes  
        # in the next higher level 
        # cur_best = max([tree[2 * pos], tree[2 * pos + 1]], key=lambda x: x.Cj)
        cur_best = max( sorted([tree[2 * pos], tree[2 * pos + 1]], key = lambda x: x.j_max, reverse = True), key=lambda x: x.Cj)
        # if tree[2 * pos].Cj == tree[2 * pos + 1].Cj:
        #     print("OMGGGG",cur_best.j_max, tree[2 * pos].j_max, tree[2 * pos + 1].j_max)
        # tree[pos].Cj = max(tree[2 * pos].Cj,  
        #                    tree[2 * pos + 1].Cj)  
        tree[pos].Cj = cur_best.Cj
        tree[pos].j_max = cur_best.j_max




def argmax(iterable):
    return max(enumerate(iterable), key=lambda x: x[1])[0]

def max_both(iterable):
    return max(enumerate(iterable), key=lambda x: x[1])


def reconstruct_solution(mems, C, trace_vector):
    solution_index = argmax(C)
    value = C[solution_index]
    # print()
    solution = []
    while solution_index > 0:
        solution.append(mems[solution_index - 1])  # trace vector is shifted on to the right so need to remove 1 from vectore to get j_index 
        solution_index = trace_vector[solution_index]
        # print(solution_index)
        # if solution_index is None:
        #     break
    return value, solution[::-1]

def reconstruct_all_solutions(mems, all_C_max_indicies, trace_vector):
    # solution_index = argmax(C)
    solutions = []
    for solution_index in all_C_max_indicies:
        value = C[solution_index]
        solution = []
        while solution_index > 0:
            solution.append(mems[solution_index - 1])  # trace vector is shifted on to the right so need to remove 1 from vectore to get j_index 
            solution_index = trace_vector[solution_index]
        solutions.append( solution[::-1] )
    return value, solutions

def make_leafs_power_of_2(mems):
    nodes = []
    nodes.append( Node(-1, -1, -2**32, -1) ) # add an start node in case a search is smaller than any d coord in tree
    for i, mem in enumerate(mems):
        # if i > 3: break
        m = Node(mem.d, mem.j, -2**32, mem.j)
        nodes.append(m)

    for i in range(20):
        if len(nodes) == 2**i or len(nodes) == 2**(i+1):
            break
        elif 2**i < len(nodes) < 2**(i+1):
            remainder = 2**(i+1) - len(nodes) 
            # print(remainder,  2**(i+1), len(nodes))
            for j in range(remainder):
                nodes.append( Node(-j-1, -j - 2, -2**32, -j - 2) ) # fill up nodes to have leaves a power of 2
            break

    leafs = sorted(nodes, key= lambda x: x.d)
    # n = len(leafs)
    return leafs


def all_solutions_c_max_indicies(C, C_max):
    return [i for i, c in enumerate(C) if c == C_max] 

## driver code for test if called as external script 

if __name__ == '__main__':
    from collections import namedtuple
    import colinear_solver 

    # construct sorted leafs
    
    def generate_mems(nr_mems):
        mem = namedtuple('Mem', ['x', 'y', 'c', 'd', 'val','j'])
        mem_choords = []
        for i in range(nr_mems):
            ref_start = random.randint(1,10000000)
            mem_length = random.randint(10,10)
            ref_stop = ref_start + mem_length
            read_start = random.randint(1,10000000)
            m = (ref_start, ref_stop,  read_start, read_start + mem_length, mem_length)
            mem_choords.append(m)
        mem = namedtuple('Mem', ['x', 'y', 'c', 'd', 'val','j'])
        mems = []
        for i, m in enumerate(sorted(mem_choords, key=lambda x: x[1])):
            m = mem(m[0], m[1],  m[2], m[3], m[4], i)
            mems.append(m)
        return mems

    mems = generate_mems(100)
    
    mem = namedtuple('Mem', ['x', 'y', 'c', 'd', 'val','j'])
    mems = [ mem(1, 10,  1, 10, 10, 0), mem(11, 20,  11, 20, 10, 1),
            mem(21, 30,  1, 10, 10, 2), mem(31, 40,  11, 20, 10, 3)]

    
    st = time()
    solution_quadratic, mem_solution_value = colinear_solver.read_coverage(mems)
    for sol in solution_quadratic:
        print(mem_solution_value, [mem.j for mem in sol])
    print("Total time quadratic method:", time()- st)
    print()


    st = time()
    T_leafs = make_leafs_power_of_2(mems)
    I_leafs = make_leafs_power_of_2(mems)
    n = len(T_leafs)
    mem_to_leaf_index = {l.j : i for i,l in enumerate(T_leafs)}
    T = [0 for i in range(2 * n) ]  
    I = [0 for i in range(2 * n) ]  
    construct_tree(T, T_leafs, n)
    construct_tree(I, I_leafs, n)
    time_construct = time()- st
    print("Time constructing tree:", time_construct)


    ############  T and I ###################
    #########################################
    st = time()
    # leafs = sorted(copy.deepcopy(nodes), key= lambda x: x.d)
    # I_leafs = copy.deepcopy(leafs)
    # for leaf in I_leafs:
    #     leaf.Cj -= leaf.d 
    # print([l.Cj for l in I_leafs])

    C = [0]* (len(mems) + 1) #(len(leafs))
    trace_vector = [None]*(len(mems) + 1)
    # print("leaf index (first mem):",mem_to_leaf_index[0],mem_to_leaf_index[-1],mem_to_leaf_index[1], "total leafs",len(mem_to_leaf_index), "len mems:", len(mems))
    update(T, mem_to_leaf_index[-1], 0, n) # point update 
    update(I, mem_to_leaf_index[-1], 0, n) # point update 

    for j, mem in enumerate(mems):
        # print(mem)
        # print("vals T:", [l.Cj for l in leafs])
        # print("vals I:", [l.Cj for l in I_leafs])
        leaf_to_update = mem_to_leaf_index[j]
        # print()
        c = mem.c
        T_max, j_prime_a, node_pos  = range_query(T, -1, c-1, len(T_leafs)) 
        print("C_a:",  T_max +  mem.d - mem.c + 1, j_prime_a, node_pos, leaf_to_update )
        # print("T TREE:", [(s, zz.j, zz.d, zz.Cj, zz.j_max) for s, zz in enumerate(T) if type(zz) != int])
        C_a =  T_max +  mem.d - mem.c + 1  # add the mem_length to T since disjoint

        if T_max < 0:
            print("BUG", T_max)
            sys.exit()

        
        d = mem.d
        I_max, j_prime_b, node_pos  = range_query(I, c, d, len(I_leafs))         
        print("C_b:", I_max +  mem.d, I_max, j_prime_b, node_pos, leaf_to_update )
        # print( I_max, mem.d, mems[j_prime_b].d, mems[j_prime_b])
        # print("I TREE:", [(s, zz.j, zz.d, zz.Cj, zz.j_max) for s, zz in enumerate(I) if type(zz) != int])
        C_b =  I_max +  mem.d #- mems[j_prime_b].d   # add the part of the mem that is not overlapping

        # if C_b < 0:
        #     print("BUG")
        #     sys.exit()
        # if C_a == C_b:
        #     print(j, C_b, "LOOOOL solve this by taking the j that is largest!",j_prime_a, j_prime_b)
        #     # j_prime = max(j_prime_a, j_prime_b)

        # else:
        index, value = max_both([C_a, C_b])
        if index == 0: # Updating with C_a
            j_prime = j_prime_a
        else: # Updating with C_b
            j_prime = j_prime_b

        C[j+1] = value

        if j_prime < 0: # any of the additional leaf nodes (with negative index) we add to make number of leafs 2^n
            trace_vector[j+1] = 0
        elif value == 0: # first j (i.e. j=0) 
            trace_vector[j+1]= 0
        else:
            trace_vector[j+1] = j_prime +1

        update(T, leaf_to_update, value, n) # point update 
        update(I, leaf_to_update, value - mem.d, n) # point update 

    # print(trace_vector)
    # print(C)
    # print(trace_vector)
    C_max, solution = reconstruct_solution(mems, C, trace_vector)

    all_C_max_indicies = all_solutions_c_max_indicies(C,C_max)
    # print(C, C_max)
    print("number solutions with the same score:", all_solutions_c_max_indicies(C, C_max))
    C_max2, solutions = reconstruct_all_solutions(mems, all_C_max_indicies, trace_vector)
    for sol in solutions:
        print(C_max2, [mem.j for mem in sol])

    time_find = time()- st
    # print(C)
    # print(trace_vector)
    # print(C_max, [mem.j for mem in solution])

    # print([ m2.d >= m1.d for m1, m2 in zip(solution[:-1], solution[1:]) ])
    # print([ m2.y >= m1.y for m1, m2 in zip(solution[:-1], solution[1:]) ])
    # assert all( m2.d >= m1.d for m1, m2 in zip(solution[:-1], solution[1:]) )
    # assert all( m2.y >= m1.y for m1, m2 in zip(solution[:-1], solution[1:]) )
    print("Time find nlogn solution:", time_find)
    print("total nlog n", time_construct + time_find)
    for i in range(len(solutions)):
        assert [mem.j for mem in solutions[i]] == [mem.j for mem in solution_quadratic[i]]
    # print("time querying RQ method 2:", time()- st)  




