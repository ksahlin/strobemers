//
//  index.cpp
//  cpp_course
//
//  Created by Kristoffer Sahlin on 4/21/21.
//

#include "index.hpp"
#include <iostream>
#include <math.h>       /* pow */
#include <bitset>
#include "xxhash.h"
#include "wyhash.h"

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
//
//
//
//
//
//
//
// copy from minimap2:sketch.c : Thomas Wang Hash
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


void string_to_hash_nohash(std::string &seq, std::vector<uint64_t> &string_hashes, std::vector<unsigned int> &pos_to_seq_choord, int k) {
    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    uint64_t kmask=(1ULL<<2*k) - 1;

    unsigned int hash_count = 0;
    int l;
    int i;
    uint64_t x = 0;
    for (int i = l = 0; i < seq.length(); i++) {
        int c = seq_nt4_table[(uint8_t) seq[i]];
        if (c < 4) { // not an "N" base
            x = (x << 2 | c) & kmask;                  // forward strand
            if (++l >= k) { // we find a k-mer
                string_hashes.push_back(x); // no hash
                pos_to_seq_choord.push_back( i - k + 1);
                hash_count ++;
            }
        } else {
            l = 0, x = 0; // if there is an "N", restart
        }
    }
}


void string_to_hash_wang(std::string &seq, std::vector<uint64_t> &string_hashes, std::vector<unsigned int> &pos_to_seq_choord, int k) {
    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    uint64_t kmask=(1ULL<<2*k) - 1;
    unsigned int hash_count = 0;
    int l;
    int i;
    uint64_t x = 0;
    for (int i = l = 0; i < seq.length(); i++) {

        int c = seq_nt4_table[(uint8_t) seq[i]];
        if (c < 4) { // not an "N" base
            x = (x << 2 | c) & kmask;                  // forward strand
            if (++l >= k) { // we find a k-mer
                uint64_t hash_k = hash64(x, kmask); // Thomas Wang hash
                string_hashes.push_back(hash_k);
                pos_to_seq_choord.push_back( i - k + 1);
                hash_count ++;
            }
        } else {
            l = 0, x = 0; // if there is an "N", restart
        }
    }
}

void string_to_hash_xxhash(std::string &seq, std::vector<uint64_t> &string_hashes, std::vector<unsigned int> &pos_to_seq_choord, int k) {
    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    uint64_t kmask=(1ULL<<2*k) - 1;
    unsigned int hash_count = 0;
    int l;
    int i;
    uint64_t x = 0;
    for (int i = l = 0; i < seq.length(); i++) {
        int c = seq_nt4_table[(uint8_t) seq[i]];
        if (c < 4) { // not an "N" base
            x = (x << 2 | c) & kmask;                  // forward strand
            if (++l >= k) { // we find a k-mer
                uint64_t hash_k = XXH3_64bits_withSeed(&x, sizeof(x), 0);
                string_hashes.push_back(hash_k);
                pos_to_seq_choord.push_back( i - k + 1);
                hash_count ++;
            }
        } else {
            l = 0, x = 0; // if there is an "N", restart
        }
    }
}


void string_to_hash_wyhash(std::string &seq, std::vector<uint64_t> &string_hashes, std::vector<unsigned int> &pos_to_seq_choord, int k) {
    uint64_t _wyp[4];
    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    uint64_t kmask=(1ULL<<2*k) - 1;
    unsigned int hash_count = 0;
    int l;
    int i;
    uint64_t x = 0;
    for (int i = l = 0; i < seq.length(); i++) {
        int c = seq_nt4_table[(uint8_t) seq[i]];
        if (c < 4) { // not an "N" base
            x = (x << 2 | c) & kmask;                  // forward strand
            if (++l >= k) { // we find a k-mer
                uint64_t hash_k = wyhash(&x, sizeof(x), 0, _wyp);
//                uint64_t hash_k = XXH3_64bits_withSeed(&x, sizeof(x), 0);
                string_hashes.push_back(hash_k);
                pos_to_seq_choord.push_back( i - k + 1);
                hash_count ++;
            }
        } else {
            l = 0, x = 0; // if there is an "N", restart
        }
    }
}





inline void get_next_strobe_shen(std::vector<uint64_t> &string_hashes, uint64_t strobe_hashval, unsigned int &strobe_pos_next, uint64_t &strobe_hashval_next,  unsigned int w_start, unsigned int w_end, uint64_t q){
    uint64_t min_val = UINT64_MAX;
    for (auto i = w_start; i <= w_end; i++) {
        uint64_t res = (strobe_hashval + string_hashes[i]) & q ;
        if (res < min_val){
            min_val = res;
            strobe_pos_next = i;
            strobe_hashval_next = string_hashes[i];
        }
    }
}


mers_vector link_2_strobes_shen(int w_min, int w_max, std::vector<uint64_t> &string_hashes, std::vector<unsigned int> &pos_to_seq_choord, unsigned int ref_index)
{
    mers_vector randstrobes2;
    uint64_t q = pow (2, 16) - 1;

    if (string_hashes.size() == 0) {
        return randstrobes2;
    }
//    std::cout << seq << std::endl;
    int seq_length = string_hashes.size();
    // create the randstrobes
    for (unsigned int i = 0; i <= seq_length; i++) {

//        if ((i % 1000000) == 0 ){
//            std::cout << i << " strobemers created." << std::endl;
//        }
        unsigned int strobe_pos_next;
        uint64_t strobe_hashval_next;

        if (i + w_max <= seq_length - 1){
            unsigned int w_start = i+w_min;
            unsigned int w_end = i+w_max;

            uint64_t strobe_hash;
            strobe_hash = string_hashes[i];
            get_next_strobe_shen(string_hashes, strobe_hash, strobe_pos_next, strobe_hashval_next, w_start, w_end, q);
        }
        else if ((i + w_min + 1 < seq_length) && (seq_length <= i + w_max) ){
            unsigned int w_start = i+w_min;
            unsigned int w_end = seq_length -1;
            uint64_t strobe_hash;
            strobe_hash = string_hashes[i];
            get_next_strobe_shen(string_hashes, strobe_hash, strobe_pos_next, strobe_hashval_next, w_start, w_end, q);
        }
        else{
            return randstrobes2;
        }

//        uint64_t hash_randstrobe2 = (string_hashes[i] << k) ^ strobe_hashval_next;
        uint64_t hash_randstrobe2 = (string_hashes[i]/2) + (strobe_hashval_next/3);

        unsigned int seq_pos_strobe1 = pos_to_seq_choord[i];
        unsigned int seq_pos_strobe2 = pos_to_seq_choord[strobe_pos_next];
        std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int> s (hash_randstrobe2, ref_index, seq_pos_strobe1, seq_pos_strobe2, seq_pos_strobe2);
        randstrobes2.push_back(s);


//        auto strobe1 = seq.substr(i, k);
//        std::cout << std::string(i, ' ') << strobe1 << std::string(strobe_pos_next - (i+k), ' ') << std::string(k, 'X') << std::endl;

    }
    return randstrobes2;
}



inline void get_next_strobe_sahlin1(std::vector<uint64_t> &string_hashes, uint64_t strobe_hashval, unsigned int &strobe_pos_next, uint64_t &strobe_hashval_next,  unsigned int w_start, unsigned int w_end, uint64_t p){
    uint64_t min_val = UINT64_MAX;
    for (auto i = w_start; i <= w_end; i++) {
        uint64_t res = (strobe_hashval + string_hashes[i]) % p ;
        if (res < min_val){
            min_val = res;
            strobe_pos_next = i;

            strobe_hashval_next = string_hashes[i];
        }
    }
}

mers_vector link_2_strobes_sahlin1(int w_min, int w_max, std::vector<uint64_t> &string_hashes, std::vector<unsigned int> &pos_to_seq_choord, unsigned int ref_index)
{
    mers_vector randstrobes2;
    uint64_t p = 997;

    if (string_hashes.size() == 0) {
        return randstrobes2;
    }
//    std::cout << seq << std::endl;
    int seq_length = string_hashes.size();
    // create the randstrobes
    for (unsigned int i = 0; i <= seq_length; i++) {

//        if ((i % 1000000) == 0 ){
//            std::cout << i << " strobemers created." << std::endl;
//        }
        unsigned int strobe_pos_next;
        uint64_t strobe_hashval_next;

        if (i + w_max <= seq_length - 1){
            unsigned int w_start = i+w_min;
            unsigned int w_end = i+w_max;
            uint64_t strobe_hash;
            strobe_hash = string_hashes[i];
            get_next_strobe_sahlin1(string_hashes, strobe_hash, strobe_pos_next, strobe_hashval_next, w_start, w_end, p);
        }
        else if ((i + w_min + 1 < seq_length) && (seq_length <= i + w_max) ){
            unsigned int w_start = i+w_min;
            unsigned int w_end = seq_length -1;
            uint64_t strobe_hash;
            strobe_hash = string_hashes[i];
            get_next_strobe_sahlin1(string_hashes, strobe_hash, strobe_pos_next, strobe_hashval_next, w_start, w_end, p);
        }
        else{
            return randstrobes2;

        }

//        uint64_t hash_randstrobe2 = (string_hashes[i] << k) ^ strobe_hashval_next;
        uint64_t hash_randstrobe2 = (string_hashes[i]/2) + (strobe_hashval_next/3);

        unsigned int seq_pos_strobe1 = pos_to_seq_choord[i];
        unsigned int seq_pos_strobe2 = pos_to_seq_choord[strobe_pos_next];
        std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int> s (hash_randstrobe2, ref_index, seq_pos_strobe1, seq_pos_strobe2, seq_pos_strobe2);
        randstrobes2.push_back(s);


//        auto strobe1 = seq.substr(i, k);
//        std::cout << std::string(i, ' ') << strobe1 << std::string(strobe_pos_next - (i+k), ' ') << std::string(k, 'X') << std::endl;

    }
    return randstrobes2;
}


inline void get_next_strobe_sahlin2(std::vector<uint64_t> &string_hashes, uint64_t strobe_hashval, unsigned int &strobe_pos_next, uint64_t &strobe_hashval_next,  unsigned int w_start, unsigned int w_end, uint64_t q){
    uint64_t min_val = UINT64_MAX;
    std::bitset<64> b;
    for (auto i = w_start; i <= w_end; i++) {
        b = (strobe_hashval ^ string_hashes[i]) & q;
        uint64_t res = b.count();
        if (res < min_val){
            min_val = res;
            strobe_pos_next = i;
            strobe_hashval_next = string_hashes[i];
        }
    }
}

mers_vector link_2_strobes_sahlin2(int w_min, int w_max, std::vector<uint64_t> &string_hashes, std::vector<unsigned int> &pos_to_seq_choord, unsigned int ref_index)
{
    mers_vector randstrobes2;
    uint64_t q = pow (2, 8) - 1;
    if (string_hashes.size() == 0) {
        return randstrobes2;
    }
    //
//    std::cout << seq << std::endl;
    int seq_length = string_hashes.size();
    // create the randstrobes
    for (unsigned int i = 0; i <= seq_length; i++) {

//        if ((i % 1000000) == 0 ){
//            std::cout << i << " strobemers created." << std::endl;
//        }
        unsigned int strobe_pos_next;
        uint64_t strobe_hashval_next;

        if (i + w_max <= seq_length - 1){
            unsigned int w_start = i+w_min;
            unsigned int w_end = i+w_max;
            uint64_t strobe_hash;
            strobe_hash = string_hashes[i];
            get_next_strobe_sahlin2(string_hashes, strobe_hash, strobe_pos_next, strobe_hashval_next, w_start, w_end,q);
        }
        else if ((i + w_min + 1 < seq_length) && (seq_length <= i + w_max) ){
            unsigned int w_start = i+w_min;
            unsigned int w_end = seq_length -1;
            uint64_t strobe_hash;
            strobe_hash = string_hashes[i];
            get_next_strobe_sahlin2(string_hashes, strobe_hash, strobe_pos_next, strobe_hashval_next, w_start, w_end,q);
        }
        else{
            return randstrobes2;
        }

//        uint64_t hash_randstrobe2 = (string_hashes[i] << k) ^ strobe_hashval_next;
        uint64_t hash_randstrobe2 = (string_hashes[i]/2) + (strobe_hashval_next/3);

        unsigned int seq_pos_strobe1 = pos_to_seq_choord[i];
        unsigned int seq_pos_strobe2 = pos_to_seq_choord[strobe_pos_next];
        std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int> s (hash_randstrobe2, ref_index, seq_pos_strobe1, seq_pos_strobe2, seq_pos_strobe2);
        randstrobes2.push_back(s);


//        auto strobe1 = seq.substr(i, k);
//        std::cout << std::string(i, ' ') << strobe1 << std::string(strobe_pos_next - (i+k), ' ') << std::string(k, 'X') << std::endl;


    }
    return randstrobes2;
}


inline void get_next_strobe_guo_pibri(std::vector<uint64_t> &string_hashes, uint64_t strobe_hashval, unsigned int &strobe_pos_next, uint64_t &strobe_hashval_next,  unsigned int w_start, unsigned int w_end){
    uint64_t min_val = UINT64_MAX;
    std::bitset<64> b;
    for (auto i = w_start; i <= w_end; i++) {
        uint64_t res = (strobe_hashval ^ string_hashes[i]);
        if (res < min_val){
            min_val = res;
            strobe_pos_next = i;
            strobe_hashval_next = string_hashes[i];
        }
    }
}

mers_vector link_2_strobes_guo_pibri(int w_min, int w_max, std::vector<uint64_t> &string_hashes, std::vector<unsigned int> &pos_to_seq_choord, unsigned int ref_index)
{
    mers_vector randstrobes2;

    if (string_hashes.size() == 0) {
        return randstrobes2;
    }
//    std::cout << seq << std::endl;
    int seq_length = string_hashes.size();
    // create the randstrobes
    for (unsigned int i = 0; i <= seq_length; i++) {

//        if ((i % 1000000) == 0 ){
//            std::cout << i << " strobemers created." << std::endl;
//        }
        unsigned int strobe_pos_next;
        uint64_t strobe_hashval_next;

        if (i + w_max <= seq_length - 1){
            unsigned int w_start = i+w_min;
            unsigned int w_end = i+w_max;

            uint64_t strobe_hash;
            strobe_hash = string_hashes[i];
            get_next_strobe_guo_pibri(string_hashes, strobe_hash, strobe_pos_next, strobe_hashval_next, w_start, w_end);
        }
        else if ((i + w_min + 1 < seq_length) && (seq_length <= i + w_max) ){
            unsigned int w_start = i+w_min;
            unsigned int w_end = seq_length -1;
            uint64_t strobe_hash;
            strobe_hash = string_hashes[i];
            get_next_strobe_guo_pibri(string_hashes, strobe_hash, strobe_pos_next, strobe_hashval_next, w_start, w_end);
        }
        else{
            return randstrobes2;
        }

//        uint64_t hash_randstrobe2 = (string_hashes[i] << k) ^ strobe_hashval_next;
        uint64_t hash_randstrobe2 = (string_hashes[i]/2) + (strobe_hashval_next/3);

        unsigned int seq_pos_strobe1 = pos_to_seq_choord[i];
        unsigned int seq_pos_strobe2 = pos_to_seq_choord[strobe_pos_next];
        std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int> s (hash_randstrobe2, ref_index, seq_pos_strobe1, seq_pos_strobe2, seq_pos_strobe2);
        randstrobes2.push_back(s);


//        auto strobe1 = seq.substr(i, k);
//        std::cout << std::string(i, ' ') << strobe1 << std::string(strobe_pos_next - (i+k), ' ') << std::string(k, 'X') << std::endl;

    }
    return randstrobes2;
}


typedef struct { uint64_t high; uint64_t low; } int128;


inline void get_next_strobe_liu_patro_li(std::vector<uint64_t> &string_hashes, uint64_t strobe_hashval, unsigned int &strobe_pos_next, uint64_t &strobe_hashval_next,  unsigned int w_start, unsigned int w_end){
    uint64_t min_val = UINT64_MAX;
    // uint64_t _wyp[4];
    int128 strobeconcat;
    strobeconcat.high = strobe_hashval;
    for (auto i = w_start; i <= w_end; i++) {
        strobeconcat.low = string_hashes[i];
        uint64_t res = XXH3_64bits_withSeed(&strobeconcat, sizeof(int128), 0);


        if (res < min_val){
            min_val = res;
            strobe_pos_next = i;
            strobe_hashval_next = string_hashes[i];
        }
    }
}


inline void get_next_strobe_liu_patro_li_wyhash(std::vector<uint64_t> &string_hashes, uint64_t strobe_hashval, unsigned int &strobe_pos_next, uint64_t &strobe_hashval_next,  unsigned int w_start, unsigned int w_end){
    uint64_t min_val = UINT64_MAX;
    int128 strobeconcat;
    uint64_t _wyp[4];
    srand(strobe_hashval);
    for (int i = 0; i < 4; i++)
        _wyp[i] = rand();

    strobeconcat.high = strobe_hashval;
    for (auto i = w_start; i <= w_end; i++) {
        strobeconcat.low = string_hashes[i];
        uint64_t res = wyhash(&strobeconcat, sizeof(int128), 0, _wyp);
        // std::cout << res << ' ' << _wyp[3] << ' ' << strobeconcat.high << std::endl;

        if (res < min_val){
            min_val = res;
            strobe_pos_next = i;
            strobe_hashval_next = string_hashes[i];
        }
    }
}


mers_vector link_2_strobes_liu_patro_li(int w_min, int w_max, std::vector<uint64_t> &string_hashes, std::vector<unsigned int> &pos_to_seq_choord, unsigned int ref_index)
{
    mers_vector randstrobes2;

    if (string_hashes.size() == 0) {
        return randstrobes2;
    }
//    std::cout << seq << std::endl;
    int seq_length = string_hashes.size();
    // create the randstrobes
    for (unsigned int i = 0; i <= seq_length; i++) {

//        if ((i % 1000000) == 0 ){
//
//            std::cout << i << " strobemers created." << std::endl;
//        }
        unsigned int strobe_pos_next;
        uint64_t strobe_hashval_next;

        if (i + w_max <= seq_length - 1){
            unsigned int w_start = i+w_min;
            unsigned int w_end = i+w_max;
            uint64_t strobe_hash;
            strobe_hash = string_hashes[i];
            get_next_strobe_liu_patro_li(string_hashes, strobe_hash, strobe_pos_next, strobe_hashval_next, w_start, w_end);
        }
        else if ((i + w_min + 1 < seq_length) && (seq_length <= i + w_max) ){
            unsigned int w_start = i+w_min;
            unsigned int w_end = seq_length -1;
            uint64_t strobe_hash;
            strobe_hash = string_hashes[i];
            get_next_strobe_liu_patro_li(string_hashes, strobe_hash, strobe_pos_next, strobe_hashval_next, w_start, w_end);
        }
        else{
            return randstrobes2;
        }

//        uint64_t hash_randstrobe2 = (string_hashes[i] << k) ^ strobe_hashval_next;
        uint64_t hash_randstrobe2 = (string_hashes[i]/2) + (strobe_hashval_next/3);

        unsigned int seq_pos_strobe1 = pos_to_seq_choord[i];
        unsigned int seq_pos_strobe2 = pos_to_seq_choord[strobe_pos_next];
        std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int> s (hash_randstrobe2, ref_index, seq_pos_strobe1, seq_pos_strobe2, seq_pos_strobe2);
        randstrobes2.push_back(s);


//        auto strobe1 = seq.substr(i, k);
//        std::cout << std::string(i, ' ') << strobe1 << std::string(strobe_pos_next - (i+k), ' ') << std::string(k, 'X') << std::endl;

    }
    return randstrobes2;
}



mers_vector link_2_strobes_liu_patro_li_wyhash(int w_min, int w_max, std::vector<uint64_t> &string_hashes, std::vector<unsigned int> &pos_to_seq_choord, unsigned int ref_index)
{
    mers_vector randstrobes2;

    if (string_hashes.size() == 0) {
        return randstrobes2;
    }
//    std::cout << seq << std::endl;
    int seq_length = string_hashes.size();
    // create the randstrobes
    for (unsigned int i = 0; i <= seq_length; i++) {

//        if ((i % 1000000) == 0 ){
//
//            std::cout << i << " strobemers created." << std::endl;
//        }
        unsigned int strobe_pos_next;
        uint64_t strobe_hashval_next;

        if (i + w_max <= seq_length - 1){
            unsigned int w_start = i+w_min;
            unsigned int w_end = i+w_max;
            uint64_t strobe_hash;
            strobe_hash = string_hashes[i];
            get_next_strobe_liu_patro_li_wyhash(string_hashes, strobe_hash, strobe_pos_next, strobe_hashval_next, w_start, w_end);
        }
        else if ((i + w_min + 1 < seq_length) && (seq_length <= i + w_max) ){
            unsigned int w_start = i+w_min;
            unsigned int w_end = seq_length -1;
            uint64_t strobe_hash;
            strobe_hash = string_hashes[i];
            get_next_strobe_liu_patro_li_wyhash(string_hashes, strobe_hash, strobe_pos_next, strobe_hashval_next, w_start, w_end);
        }
        else{
            return randstrobes2;
        }

//        uint64_t hash_randstrobe2 = (string_hashes[i] << k) ^ strobe_hashval_next;
        uint64_t hash_randstrobe2 = (string_hashes[i]/2) + (strobe_hashval_next/3);

        unsigned int seq_pos_strobe1 = pos_to_seq_choord[i];
        unsigned int seq_pos_strobe2 = pos_to_seq_choord[strobe_pos_next];
        std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int> s (hash_randstrobe2, ref_index, seq_pos_strobe1, seq_pos_strobe2, seq_pos_strobe2);
        randstrobes2.push_back(s);


//        auto strobe1 = seq.substr(i, k);
//        std::cout << std::string(i, ' ') << strobe1 << std::string(strobe_pos_next - (i+k), ' ') << std::string(k, 'X') << std::endl;

    }
    return randstrobes2;
}


//
//mers_vector link_3_strobes_method2(int w_min, int w_max, std::vector<uint64_t> &string_hashes, std::vector<unsigned int> &pos_to_seq_choord, unsigned int ref_index)
//{
//    mers_vector randstrobes3;
//    uint64_t q = pow (2, 16) - 1;
//    if (string_hashes.size() == 0) {
//        return randstrobes3;
//    }
////    std::cout << seq << std::endl;
//    int seq_length = string_hashes.size();
//    // create the randstrobes
//    for (unsigned int i = 0; i <= seq_length; i++) {
//
////        if ((i % 1000000) == 0 ){
////            std::cout << i << " randstrobes created." << std::endl;
////        }
//        uint64_t strobe_hash;
//        strobe_hash = string_hashes[i];
//
//        unsigned int strobe_pos_next1;
//        uint64_t strobe_hashval_next1;
//        unsigned int strobe_pos_next2;
//        uint64_t strobe_hashval_next2;
//
//        if (i + 2*w_max <= seq_length - 1){
//            unsigned int w1_start = i+w_min;
//            unsigned int w1_end = i+w_max;
//            get_next_strobe(string_hashes, strobe_hash, strobe_pos_next1, strobe_hashval_next1, w1_start, w1_end, q);
//
//            unsigned int w2_start = i+w_max + w_min;
//            unsigned int w2_end = i+2*w_max;
////            uint64_t conditional_next = strobe_hash ^ strobe_hashval_next1;
//            get_next_strobe(string_hashes, strobe_hashval_next1, strobe_pos_next2, strobe_hashval_next2, w2_start, w2_end, q);
//        }
//        else if ((i + 2*w_min + 1 < seq_length) && (seq_length <= i + 2*w_max) ){
////            unsigned int w_start = i+w_min;
////            unsigned int w_end = seq_length -1;
////            uint64_t strobe_hash;
////            strobe_hash = string_hashes[i];
////            get_next_strobe(string_hashes, strobe_hash, strobe_pos_next1, strobe_hashval_next1, w_start, w_end, q);
//
//            int overshot;
//            overshot = i + 2*w_max - seq_length;
//            unsigned int w1_start = i+w_min;
//            unsigned int w1_end = i+w_max - overshot/2;
//            get_next_strobe(string_hashes, strobe_hash, strobe_pos_next1, strobe_hashval_next1, w1_start, w1_end, q);
//
//            unsigned int w2_start = i+w_max - overshot/2 + w_min;
//            unsigned int w2_end = i+2*w_max - overshot;
////            uint64_t conditional_next = strobe_hash ^ strobe_hashval_next1;
//            get_next_strobe(string_hashes, strobe_hashval_next1, strobe_pos_next2, strobe_hashval_next2, w2_start, w2_end, q);
//        }
//        else{
//            return randstrobes3;
//        }
//
////        uint64_t hash_randstrobe3 = (((strobe_hash << k) ^ strobe_hashval_next1) << k) ^ strobe_hashval_next2;
//        uint64_t hash_randstrobe3 = (strobe_hash/3) + (strobe_hashval_next1/4) + (strobe_hashval_next2/5);
//
//        unsigned int seq_pos_strobe1 = pos_to_seq_choord[i];
//        unsigned int seq_pos_strobe2 =  pos_to_seq_choord[strobe_pos_next1]; //seq_pos_strobe1 + (strobe_pos_next1 - i); //
//        unsigned int seq_pos_strobe3 =  pos_to_seq_choord[strobe_pos_next2]; //seq_pos_strobe1 + (strobe_pos_next2 - i); //
////        std::cout << i << " " << strobe_pos_next1 << " " << strobe_pos_next2 << " " << seq_pos_strobe1 << " " << seq_pos_strobe2 << " " << seq_pos_strobe3 << " " << pos_to_seq_choord.size() << std::endl;
//
//        // TODO: Take care of corner case (tmep if statement below. Some values in end of string produce a cororidnate of 0 for the last strobe. Probably an off-by-one error in the calculation of the strobe coord in the last strobe window
//        if (strobe_pos_next2 ==  seq_length){
////            std::cout << "OMGGGGGGG " << i << " " << seq_pos_strobe1 << " " << seq_pos_strobe2 << " " << seq_pos_strobe3 << std::endl;
//            seq_pos_strobe3 = seq_length-1;
//        }
//        std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int> s (hash_randstrobe3, ref_index, seq_pos_strobe1, seq_pos_strobe2, seq_pos_strobe3);
//        randstrobes3.push_back(s);
//
//
////        auto strobe1 = seq.substr(i, k);
////        std::cout << std::string(i, ' ') << strobe1 << std::string(strobe_pos_next1 - (i+k), ' ') << std::string(k, 'X') << std::string(strobe_pos_next2 - strobe_pos_next1 - k, ' ') << std::string(k, 'X') << std::endl;
////        std::cout << i << " " << strobe_pos_next1 << " " << strobe_pos_next2 << " " << seq_length << std::endl;
//    }
//    return randstrobes3;
//}


void process_flat_vector(mers_vector &flat_vector, uint64_t &unique_elements){

    //    flat_array sort
    std::sort(flat_vector.begin(), flat_vector.end());

    uint64_t prev_k;
    std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int> t = flat_vector[0];
    prev_k = std::get<0>(t);
    uint64_t curr_k;
    unique_elements = 1;
    for ( auto &t : flat_vector ) {
//        std::cout << t << std::endl;
        curr_k = std::get<0>(t);
        if (curr_k != prev_k){
            unique_elements ++;
        }
        prev_k = curr_k;
    }

//    return flat_vector;
}


unsigned int index_vector(mers_vector &flat_vector, kmer_lookup &mers_index, float f){

    std::cout << "Flat vector size: " << flat_vector.size() << std::endl;
//    kmer_lookup mers_index;
    uint64_t offset = 0;
    uint64_t prev_offset = 0;
    unsigned int count = 0;

    unsigned int tot_occur_once = 0;
    unsigned int tot_high_ab = 0;
    unsigned int tot_mid_ab = 0;
    std::vector<unsigned int> strobemer_counts;

    uint64_t prev_k;
    std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int> t = flat_vector[0];
    prev_k = std::get<0>(t);
    uint64_t curr_k;

    for ( auto &t : flat_vector ) {
//        std::cout << t << std::endl;
        curr_k = std::get<0>(t);
        if (curr_k == prev_k){
            count ++;
        }
        else {
            if (count == 1){
                tot_occur_once ++;
            }
            else if (count > 100){
                tot_high_ab ++;
//                std::cout << count << std::endl;
            }
            else{
                tot_mid_ab ++;
            }
            strobemer_counts.push_back(count);

            std::tuple<unsigned int, unsigned int> s(prev_offset, count);
            mers_index[prev_k] = s;
            count = 1;
            prev_k = curr_k;
            prev_offset = offset;
        }
        offset ++;
    }

    // last k-mer
    std::tuple<unsigned int, unsigned int> s(prev_offset, count);
    mers_index[curr_k] = s;

    std::cout << "Total strobemers count: " << offset << std::endl;
    std::cout << "Total strobemers occur once: " << tot_occur_once << std::endl;
    std::cout << "Total strobemers highly abundant > 100: " << tot_high_ab << std::endl;
    std::cout << "Total strobemers mid abundance (between 2-100): " << tot_mid_ab << std::endl;
    std::cout << "Total distinct strobemers stored: " << mers_index.size() << std::endl;
    if (tot_high_ab >= 1) {
        std::cout << "Ratio distinct to highly abundant: " << mers_index.size() / tot_high_ab << std::endl;
    }
    if (tot_mid_ab >= 1) {
        std::cout << "Ratio distinct to non distinct: " << mers_index.size() / (tot_high_ab + tot_mid_ab) << std::endl;
    }
    // get count for top -f fraction of strobemer count to filter them out
    std::sort(strobemer_counts.begin(), strobemer_counts.end(), std::greater<int>());

    unsigned int index_cutoff = strobemer_counts.size()*f;
    std::cout << "Filtered cutoff index: " << index_cutoff << std::endl;
    unsigned int filter_cutoff =  strobemer_counts[index_cutoff];
    filter_cutoff = filter_cutoff > 30 ? filter_cutoff : 30; // cutoff is around 30-50 on hg38. No reason to have a lower cutoff than this if aligning to a smaller genome or contigs.
    std::cout << "Filtered cutoff count: " << filter_cutoff << std::endl;
    std::cout << "" << std::endl;
    std::cout << "" << std::endl;
    return filter_cutoff;
}





