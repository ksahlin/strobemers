//
//  index.hpp
//  cpp_course
//
//  Created by Kristoffer Sahlin on 4/21/21.
//

#ifndef index_hpp
#define index_hpp

#include <stdio.h>
#include <string>
#include <vector>
#include <deque>
#include "robin_hood.h"

uint64_t hash(std::string kmer);
static inline uint64_t hash64(uint64_t key, uint64_t mask);
typedef robin_hood::unordered_map< uint64_t , std::vector<unsigned int> > seq_index;
void generate_kmer_index(seq_index &h, int k, std::string &seq, unsigned int ref_index);

void generate_minstrobe2_index(seq_index &h, int n, int k, int w_min, int w_max, std::string &seq, unsigned int ref_index);

struct strobemer2 {
    uint64_t hashval;
    unsigned short int p1;
    unsigned short int p2;
};

struct kmer {
    uint64_t hashval;
    unsigned short int p1;
};

#endif /* index_hpp */



