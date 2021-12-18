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
#include "xxhash.h"

uint64_t hash(std::string kmer);
static inline uint64_t hash64(uint64_t key, uint64_t mask);

typedef std::vector< std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int>> mers_vector;

void string_to_hash_wang(std::string &seq, std::vector<uint64_t> &string_hashes, std::vector<unsigned int> &pos_to_seq_choord, int k);
void string_to_hash_xxhash(std::string &seq, std::vector<uint64_t> &string_hashes, std::vector<unsigned int> &pos_to_seq_choord, int k);
void string_to_hash_nohash(std::string &seq, std::vector<uint64_t> &string_hashes, std::vector<unsigned int> &pos_to_seq_choord, int k);

static inline void get_next_strobe(std::vector<uint64_t> &string_hashes, uint64_t strobe_hashval, unsigned int &strobe_pos_next, uint64_t &strobe_hashval_next,  unsigned int w_start, unsigned int w_end, uint64_t q);
mers_vector link_2_strobes_method2(int w_min, int w_max, std::vector<uint64_t> &string_hashes, std::vector<unsigned int> &pos_to_seq_choord, unsigned int ref_index);
mers_vector link_3_strobes_method2(int w_min, int w_max, std::vector<uint64_t> &string_hashes, std::vector<unsigned int> &pos_to_seq_choord, unsigned int ref_index);

typedef robin_hood::unordered_map< unsigned int, std::vector< std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int>>> three_pos_index;
//mers_vector construct_flat_vector_three_pos(three_pos_index &tmp_index, uint64_t &unique_elements);
typedef robin_hood::unordered_map< uint64_t, std::tuple<uint64_t, unsigned int >> kmer_lookup;
void process_flat_vector(mers_vector &flat_vector, uint64_t &unique_elements);
unsigned int index_vector(mers_vector  &mers_vector, kmer_lookup &mers_index, float f);


struct hit {
    unsigned int query_s;
    unsigned int query_e;
    unsigned int ref_s;
    unsigned int ref_e;
};

struct nam {
    unsigned int ref_id;
    unsigned int query_s;
    unsigned int query_e;
    unsigned int ref_s;
    unsigned int ref_e;
    unsigned int n_hits = 0;
    unsigned int previous_query_start;
    unsigned int previous_ref_start;
//    uint64_t copy_id; // If many hits, keep track of which it in order of left to right on the reference
};
#endif /* index_hpp */



