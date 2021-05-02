//
//  index.cpp
//  cpp_course
//
//  Created by Kristoffer Sahlin on 4/21/21.
//

#include "index.hpp"
#include <iostream>
#include <math.h>       /* pow */

/**********************************************************
 *
 * hash kmer into uint64
 *
 * *******************************************************/
// copy from http://www.cse.yorku.ca/~oz/hash.html:

uint64_t hash(std::string kmer)
{
    unsigned long hash = 5381;
    int c;
    for (std::string::size_type i=0; i< kmer.length(); i++) {
        c = kmer[i];
        hash = ((hash << 5) + hash) + c; /* hash * 33 + c */
    }
    return hash;
}

/**********************************************************
 *
 * hash kmer into uint64
 *
 * *******************************************************/
// copy from minimap2:sketch.c :
static inline uint64_t hash64(uint64_t key, uint64_t mask)
{
    key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
    key = key ^ key >> 24;
    key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
    key = key ^ key >> 14;
    key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
    key = key ^ key >> 28;
    key = (key + (key << 31)) & mask;
    return key;
}//hash64


static unsigned char seq_nt4_table[256] = {
        0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
}; //seq_nt4_table

static inline uint64_t kmer_to_uint64(std::string &kmer, uint64_t kmask)
{
    uint64_t bkmer = 0;

    for (char i : kmer) {
        int c = seq_nt4_table[(uint8_t)i];
        bkmer = (bkmer << 2 | c) & kmask;

    }
    return bkmer;
}

// std::vector<strobemer2>() seq_to_hybridstrobes2(int n, int k, int w_min, int w_max, std::string &seq, unsigned int ref_index)
// std::vector<strobemer3>() seq_to_hybridstrobes3(int n, int k, int w_min, int w_max, std::string &seq, unsigned int ref_index)
// std::vector<kmer>() seq_to_kmers(int n, int k, int w_min, int w_max, std::string &seq, unsigned int ref_index)


void generate_kmer_index(seq_index1 &h, int k, std::string &seq, unsigned int ref_index)
{
    robin_hood::hash<std::string> robin_hash;
    uint64_t kmask=(1ULL<<2*k) - 1;
    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    for (int i = 0; i <= seq.length() - k; i++) {
        auto kmer = seq.substr(i, k);
//        uint64_t bkmer = kmer_to_uint64(kmer, kmask);
//        uint64_t kmer_hashval = hash64(bkmer, kmask);
        uint64_t kmer_hashval;
        kmer_hashval = robin_hash(kmer);
//        uint64_t kmer3;
//        kmer3 = hash(kmer);
        if (h.find(kmer_hashval) == h.end()){ // Not  in  index
            h[kmer_hashval] = std::vector< std::tuple<unsigned int, unsigned int> >();  //std::vector<unsigned int>(); // Initialize key with null vector
            std::tuple<unsigned int, unsigned int> s (ref_index, i);
            h[kmer_hashval].push_back(s);

        }
        else{
            std::tuple<unsigned int, unsigned int> s (ref_index, i);
            h[kmer_hashval].push_back(s);
        }

    }
}


// initialize queue and current minimum and position
static inline void initialize_window(std::string &seq, std::deque <uint64_t> &q, uint64_t &q_min_val, int &q_min_pos, int w_min, int w_max, int k, uint64_t &kmask){
    robin_hood::hash<std::string> robin_hash;
    for (int i = w_min; i < w_max; i++) {
        auto strobe = seq.substr(i, k);
        uint64_t strobe_hashval;
        strobe_hashval = robin_hash(strobe);
//        uint64_t bstrobe = kmer_to_uint64(strobe, kmask);
//        uint64_t strobe_hashval = hash64(bstrobe, kmask);
//        uint64_t strobe_hashval;
//        strobe_hashval = hash(strobe);

//        std::cout << seq << std::endl;
//        std::cout << std::string(i, ' ') << strobe << " " << strobe_hashval << " " << i << std::endl;
//        std::cout << "INIT " << i << " " << strobe << " " << strobe_hashval  << std::endl;
        q.push_back(strobe_hashval);
        if (strobe_hashval < q_min_val) {
            q_min_val = strobe_hashval;
            q_min_pos = i;
        }
    }
}

// update queue and current minimum and position
static inline void update_window(std::deque <uint64_t> &q, uint64_t &q_min_val, int &q_min_pos, uint64_t new_strobe_hashval, int w_min, int w_max, int i ){
    uint64_t popped_val;
    popped_val = q.front();
    q.pop_front();
    q.push_back(new_strobe_hashval);
//    std::cout << "LOL " << popped_val << " " << q_min_val  << std::endl;
    if (popped_val == q_min_val){ // we popped the minimum value, find new brute force
//        std::cout << "OK "  << std::endl;
        q_min_val = UINT64_MAX;
        q_min_pos = -1;
//        std::cout << "OK " << q_min_val  << std::endl;
//        std::cout << "OK " << q_min_pos  << std::endl;
        for (int j = 0; j <= q.size()-1; j++) {
//            std::cout << q[j] << " " << j << " " << i + w_min  << std::endl;
            if (q[j] < q_min_val) {
                q_min_val = q[j];
                q_min_pos = i + w_min + j + 1;
            }
        }
//        std::cout << "Changed: " << q_min_pos << " " << q_min_val << std::endl;
    }
    else if ( new_strobe_hashval < q_min_val ) { // the new value added to queue is the new minimum
        q_min_val = new_strobe_hashval;
        q_min_pos = i + w_max;
    }
}


void generate_minstrobe2_index(seq_index2 &h, int n, int k, int w_min, int w_max, std::string &seq, unsigned int ref_index)
{
    if (seq.length() < w_max) {
        return;
    }

    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    robin_hood::hash<std::string> robin_hash;
    unsigned int seq_length = seq.length();
    // initialize the deque
    std::deque <uint64_t> q;
    uint64_t kmask=(1ULL<<2*k) - 1;
    uint64_t q_min_val = UINT64_MAX;
    int q_min_pos = -1;
    initialize_window(seq, q, q_min_val, q_min_pos, w_min, w_max, k, kmask);

//    std::cout << seq << std::endl;

    // create the minstrobes
    for (unsigned int i = 0; i <= seq_length - k; i++) {
        auto strobe1 = seq.substr(i, k);
//        uint64_t bstrobe = kmer_to_uint64(strobe1, kmask);
//        uint64_t strobe_hashval = hash64(bstrobe, kmask);
        uint64_t strobe_hashval;
        strobe_hashval = robin_hash(strobe1);

        uint64_t minstrobe_h = (strobe_hashval << k) ^ q_min_val;


        if (h.find(minstrobe_h) == h.end()){ // Not  in  index
            h[minstrobe_h] = std::vector< std::tuple<unsigned int, unsigned int, unsigned int> >();  //std::vector<unsigned int>(); // Initialize key with null vector
            std::tuple<unsigned int, unsigned int, unsigned int> s (ref_index, i, q_min_pos);
            h[minstrobe_h].push_back(s);
        }
        else{
            std::tuple<unsigned int, unsigned int, unsigned int> s (ref_index, i, q_min_pos);
            h[minstrobe_h].push_back(s);
        }

        // update queue and current minimum and position
        if (i + w_max <= seq_length - k){
            auto new_strobe = seq.substr(i + w_max, k);
            uint64_t new_bstrobe = kmer_to_uint64(new_strobe, kmask);
            uint64_t new_strobe_hashval = hash64(new_bstrobe, kmask);

            update_window(q, q_min_val, q_min_pos, new_strobe_hashval, w_min, w_max, i );
        }
        else if ((i + w_min + 1 < seq_length - k) && (seq_length - k < i + w_max) ){
            uint64_t new_strobe_hashval =  UINT64_MAX;
            update_window(q, q_min_val, q_min_pos, new_strobe_hashval, w_min, w_max, i );
        }
        else{
            return;
        }

//        std::cout << std::string(i, ' ') << strobe1 << std::string(q_min_pos - (i+k), ' ') << std::string(k, 'X') << std::endl;

    }
}


static inline void get_next_strobe(std::vector<uint64_t> &string_hashes, uint64_t strobe_hashval, unsigned int &strobe_pos_next, uint64_t &strobe_hashval_next,  unsigned int w_start, unsigned int w_end, uint64_t q){
    uint64_t min_val = UINT64_MAX;
    unsigned int min_pos;
    min_pos = -1;
    for (auto i = w_start; i <= w_end; i++) {
        uint64_t res = (strobe_hashval + string_hashes[i]) & q ;
        if (res < min_val){
            min_val = res;
            strobe_pos_next = i;
            strobe_hashval_next = string_hashes[i];
        }
    }
}

static inline void make_string_to_hashvalues(std::string &seq, std::vector<uint64_t> &string_hashes, int k, uint64_t kmask){
    robin_hood::hash<std::string> robin_hash;
    for (int i = 0; i <= seq.length() - k; i++) {
        auto strobe1 = seq.substr(i, k);
//        uint64_t bstrobe = kmer_to_uint64(strobe1, kmask);
//        uint64_t strobe_hashval = hash64(bstrobe, kmask);
        uint64_t strobe_hashval;
        strobe_hashval = robin_hash(strobe1);

        string_hashes.push_back(strobe_hashval);
    }

}


void generate_randstrobe2_index(seq_index2 &h, int n, int k, int w_min, int w_max, std::string &seq, unsigned int ref_index)
{
    if (seq.length() < w_max) {
        return;
    }

    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    uint64_t kmask=(1ULL<<2*k) - 1;
    uint64_t q = pow (2, 10) - 1;
    // make string of strobes into hashvalues all at once to avoid repetitive k-mer to hash value computations
    std::vector<uint64_t> string_hashes;
    make_string_to_hashvalues(seq, string_hashes, k, kmask);
    unsigned int seq_length = string_hashes.size();
    uint64_t strobe_hash = 1223;
    uint64_t randstrobe_h = 12203;
//    std::cout << seq << std::endl;

    // create the randstrobes
    for (unsigned int i = 0; i <= seq_length; i++) {

        if ((i % 1000000) == 0 ){
            std::cout << i << " strobemmers created." << std::endl;
        }
        unsigned int strobe_pos_next;
        uint64_t strobe_hashval_next;

        if (i + w_max <= seq_length){
            unsigned int w_start = i+w_min;
            unsigned int w_end = i+w_max;
            uint64_t strobe_hash;
            strobe_hash = string_hashes[i];
            get_next_strobe(string_hashes, strobe_hash, strobe_pos_next, strobe_hashval_next, w_start, w_end, q);
        }
        else if ((i + w_min + 1 < seq_length) && (seq_length < i + w_max) ){
            unsigned int w_start = i+w_min;
            unsigned int w_end = seq_length -1;
            uint64_t strobe_hash;
            strobe_hash = string_hashes[i];
            get_next_strobe(string_hashes, strobe_hash, strobe_pos_next, strobe_hashval_next, w_start, w_end, q);
        }
        else{
            return;
        }

        uint64_t randstrobe_h = (string_hashes[i] << k) ^ strobe_hashval_next;


        if (h.find(randstrobe_h) == h.end()){ // Not  in  index
            h[randstrobe_h] = std::vector< std::tuple<unsigned int, unsigned int, unsigned int> >();  //std::vector<unsigned int>(); // Initialize key with null vector
            std::tuple<unsigned int, unsigned int, unsigned int> s (ref_index, i, strobe_pos_next);
            h[randstrobe_h].push_back(s);

        }
        else{
            std::tuple<unsigned int, unsigned int, unsigned int> s (ref_index, i, strobe_pos_next);
            h[randstrobe_h].push_back(s);
        }

//        auto strobe1 = seq.substr(i, k);
//        std::cout << std::string(i, ' ') << strobe1 << std::string(strobe_pos_next - (i+k), ' ') << std::string(k, 'X') << std::endl;

    }
}



void generate_hybridstrobe2_index(seq_index2 &h, int n, int k, int w_min, int w_max, std::string &seq, unsigned int ref_index)
{
    if (seq.length() < w_max) {
        return;
    }

    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    unsigned int seq_length = seq.length();
    uint64_t kmask=(1ULL<<2*k) - 1;
    int x = (w_max - w_min) /3;

    // initialize deque1
    std::deque <uint64_t> q1;
    uint64_t q1_min_val = UINT64_MAX;
    int q1_min_pos = -1;
    initialize_window(seq, q1, q1_min_val, q1_min_pos, w_min, w_min+ x, k, kmask);

    // initialize deque2
    std::deque <uint64_t> q2;
    uint64_t q2_min_val = UINT64_MAX;
    int q2_min_pos = -1;
    initialize_window(seq, q2, q2_min_val, q2_min_pos, w_min+x, w_min+2*x, k, kmask);

    // initialize deque3
    std::deque <uint64_t> q3;
    uint64_t q3_min_val = UINT64_MAX;
    int q3_min_pos = -1;
    initialize_window(seq, q3, q3_min_val, q3_min_pos, w_min+2*x, w_max, k, kmask);

//    std::cout << seq << std::endl;

    // create the hybridstrobes
    for (unsigned int i = 0; i <= seq_length - k; i++) {
        auto strobe1 = seq.substr(i, k);
        uint64_t bstrobe = kmer_to_uint64(strobe1, kmask);
        uint64_t strobe_hashval = hash64(bstrobe, kmask);

        int strobe2_pos;
        uint64_t strobe2_val;
        uint r =  strobe_hashval % 3;
        if (r == 0){
            strobe2_val = q1_min_val;
            strobe2_pos = q1_min_pos;
        }
        else if (r == 1){
            strobe2_val = q2_min_val;
            strobe2_pos = q2_min_pos;
        }
        else{
            strobe2_val = q3_min_val;
            strobe2_pos = q3_min_pos;
        }

        uint64_t hybridstrobe_h;
        hybridstrobe_h = (strobe_hashval << k) ^ strobe2_val;


        if (h.find(hybridstrobe_h) == h.end()){ // Not  in  index
            h[hybridstrobe_h] = std::vector< std::tuple<unsigned int, unsigned int, unsigned int> >();  //std::vector<unsigned int>(); // Initialize key with null vector
            std::tuple<unsigned int, unsigned int, unsigned int> s (ref_index, i, strobe2_pos);
            h[hybridstrobe_h].push_back(s);
        }
        else{
            std::tuple<unsigned int, unsigned int, unsigned int> s (ref_index, i, strobe2_pos);
            h[hybridstrobe_h].push_back(s);
        }

        // update queue and current minimum and position
        if (i + w_max <= seq_length - k){
            auto new_strobe1 = seq.substr(i + w_min+x, k);
            uint64_t new_bstrobe1 = kmer_to_uint64(new_strobe1, kmask);
            uint64_t new_strobe_hashval1 = hash64(new_bstrobe1, kmask);
            update_window(q1, q1_min_val, q1_min_pos, new_strobe_hashval1, w_min, w_min+x, i );

            auto new_strobe2 = seq.substr(i + w_min+2*x, k);
            uint64_t new_bstrobe2 = kmer_to_uint64(new_strobe2, kmask);
            uint64_t new_strobe_hashval2 = hash64(new_bstrobe2, kmask);
            update_window(q2, q2_min_val, q2_min_pos, new_strobe_hashval2, w_min+x, w_min+2*x, i );

            auto new_strobe3 = seq.substr(i + w_max, k);
            uint64_t new_bstrobe3 = kmer_to_uint64(new_strobe3, kmask);
            uint64_t new_strobe_hashval3 = hash64(new_bstrobe3, kmask);
            update_window(q3, q3_min_val, q3_min_pos, new_strobe_hashval3, w_min+2*x, w_max, i );
        }
        // TODO: Fix the below code to generate hybridstrobes all the way to the end of the sequence.
        //  It probably creates and empty queue which leads to an infinity loop or something.
//        else if ((i + w_min + 1 < seq_length - k) && (seq_length - k < i + w_min+x) ){
//            uint64_t new_strobe_hashval =  UINT64_MAX;
//            update_window(q1, q1_min_val, q1_min_pos, new_strobe_hashval, w_min, w_min+x, i );
//            update_window(q2, q2_min_val, q2_min_pos, new_strobe_hashval, w_min+x, w_min+2*x, i );
//            update_window(q3, q3_min_val, q3_min_pos, new_strobe_hashval, w_min+2*x, w_max, i );
//        }
//        else if ((i + w_min + 1 < seq_length - k) && (seq_length - k < i + w_min+2*x) ){
//            uint64_t new_strobe_hashval =  UINT64_MAX;
//            update_window(q2, q2_min_val, q2_min_pos, new_strobe_hashval, w_min+x, w_min+2*x, i );
//            update_window(q3, q3_min_val, q3_min_pos, new_strobe_hashval, w_min+2*x, w_max, i );
//        }
//        else if ((i + w_min + 1 < seq_length - k) && (seq_length - k < i + w_max) ){
//            uint64_t new_strobe_hashval =  UINT64_MAX;
//            update_window(q3, q3_min_val, q3_min_pos, new_strobe_hashval, w_min, w_max, i );
//        }
        else{
            return;
        }

//        std::cout << std::string(i, ' ') << strobe1 << std::string(strobe2_pos - (i+k), ' ') << std::string(k, 'X') << std::endl;

    }
}

