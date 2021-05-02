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
#include <tuple>
#include "robin_hood.h"


uint64_t hash(std::string kmer);
static inline uint64_t hash64(uint64_t key, uint64_t mask);
typedef robin_hood::unordered_map< uint64_t , std::vector< std::tuple<unsigned int, unsigned int>> > seq_index1;
typedef robin_hood::unordered_map< uint64_t , std::vector< std::tuple<unsigned int, unsigned int, unsigned int>> > seq_index2;
typedef robin_hood::unordered_map< uint64_t , std::vector< std::tuple<unsigned int, unsigned int, unsigned int, unsigned int>> > seq_index3;
//typedef robin_hood::unordered_map< uint64_t , std::vector<unsigned int> > seq_index;

void generate_kmer_index(seq_index1 &h, int k, std::string &seq, unsigned int ref_index);

void generate_minstrobe2_index(seq_index2 &h, int n, int k, int w_min, int w_max, std::string &seq, unsigned int ref_index);

static inline void make_string_to_hashvalues(std::string &seq, std::vector<uint64_t> &string_hashes, int k, uint64_t kmask);
void generate_randstrobe2_index(seq_index2 &h, int n, int k, int w_min, int w_max, std::string &seq, unsigned int ref_index);
static inline void get_next_strobe(std::vector<uint64_t> &string_hashes, uint64_t strobe_hashval, unsigned int &strobe_pos_next, uint64_t &strobe_hashval_next,  unsigned int w_start, unsigned int w_end, uint64_t q);

void generate_hybridstrobe2_index(seq_index2 &h, int n, int k, int w_min, int w_max, std::string &seq, unsigned int ref_index);

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



