#include <iostream>
#include <fstream>
#include <unordered_map>
#include "robin_hood.h"
#include <vector>
#include <string>
#include <chrono>  // for high_resolution_clock

#include "index.hpp"
//#include "khashl.h" // hash table

typedef robin_hood::unordered_map< unsigned int , std::string > sequences;
typedef robin_hood::unordered_map< unsigned int, std::string > idx_to_acc;

//KHASHL_MAP_INIT(, kc_c1_t, index1, uint64_t, std::vector< std::tuple<unsigned int, unsigned int>> > , kh_hash_uint64, kh_eq_generic)


static void read_fasta(sequences &seqs, idx_to_acc &acc_map, std::string fn)
{
    std::ifstream file(fn);
    std::string line, seq;
    int ref_index = 0;
    while (getline(file, line)) {
        if (line[0] == '>') {
//            std::cout << ref_index << " " << line << std::endl;
            if (seq.length() > 0){
                seqs[ref_index -1] = seq;
//                std::cout << ref_index - 1 << " here " << seq << " " << seq.length() << " " << seq.size() << std::endl;
//                generate_kmers(h, k, seq, ref_index);
            }
            acc_map[ref_index] = line;
            ref_index++;
            seq = "";
        }
        else {
            seq += line;
        }
    }
    if (seq.length() > 0){
        seqs[ref_index -1] = seq;
//        std::cout << ref_index -1 << " here2 " << seq << std::endl;
//        generate_kmers(h, k, seq, ref_index);
    }
    file.close();
}

static void print_diagnostics_new(std::vector< std::tuple<uint64_t, unsigned int, unsigned int>>  kmers, idx_to_acc &acc_map) {
    uint64_t tot_index_size = 0;
    for (int i = 0; i < kmers.size(); ++i)
    {
        // access using []
        auto t = kmers[i];
//        std::cout << "(" << std::get<0>(t) << ", " << std::get<1>(t) << ", " << std::get<2>(t) << ")";
        tot_index_size += sizeof(t);
//        std::cout << sizeof(t) << std::endl;
    }
    std::cout << "Total size of index vector : " << tot_index_size/1000000  << " Mb." << std::endl;

}

static void print_diagnostics_new2(std::vector< std::tuple<uint64_t, unsigned int, unsigned int>> &mers_vector, robin_hood::unordered_map< uint64_t, std::tuple<uint64_t, unsigned int >> mers_index ) {
    uint64_t tot_flat_vector_size = 0;
    for (int i = 0; i < mers_vector.size(); ++i)
    {
        // access using []
        auto t = mers_vector[i];
//        std::cout << "(" << std::get<0>(t) << ", " << std::get<1>(t) << ", " << std::get<2>(t) << "), ";
        tot_flat_vector_size += sizeof(t);
    }
    std::cout << "Total size of flat kmer vector : " << tot_flat_vector_size/1000000  << " Mb." << std::endl;

    uint64_t tot_hashtable_index_size = 0;
    for (auto &it : mers_index)
    {
//        std::cout << it.first << ": (" << std::get<0>(it.second) << ", " << std::get<1>(it.second) << "), " ;
        tot_hashtable_index_size += sizeof(it.first);
        tot_hashtable_index_size += sizeof(it.second);
    }
    std::cout << "Total size of hash table index : " << tot_hashtable_index_size/1000000  << " Mb." << std::endl;
}

static void print_diagnostics_new3(std::vector< std::tuple<uint64_t, unsigned int, unsigned int, unsigned int>> &mers_vector, robin_hood::unordered_map< uint64_t, std::tuple<uint64_t, unsigned int >> mers_index ) {
    uint64_t tot_flat_vector_size = 0;
    for (int i = 0; i < mers_vector.size(); ++i)
    {
        // access using []
        auto t = mers_vector[i];
//        std::cout << "(" << std::get<0>(t) << ", " << std::get<1>(t) << ", " << std::get<2>(t) << ", " << std::get<3>(t) << "), ";
        tot_flat_vector_size += sizeof(t);
    }
    std::cout << "Total size of flat kmer vector : " << tot_flat_vector_size/1000000  << " Mb." << std::endl;

    uint64_t tot_hashtable_index_size = 0;
    for (auto &it : mers_index)
    {
//        std::cout << it.first << ": (" << std::get<0>(it.second) << ", " << std::get<1>(it.second) << "), " ;
        tot_hashtable_index_size += sizeof(it.first);
        tot_hashtable_index_size += sizeof(it.second);
    }
    std::cout << "Total size of hash table index : " << tot_hashtable_index_size/1000000  << " Mb." << std::endl;
}

static void print_diagnostics(seq_index1 &h, idx_to_acc &acc_map) {
    uint64_t tot_index_size = 0;
    for (auto &it : h)
    {
//        std::cout << "\nKey: " << it.first << ",  values: ";
        tot_index_size += sizeof(it.first);
        for (auto &t : it.second) // it.second is the vector, i is a tuple
        {
//            std::cout << "tuple: " << std::get<0>(t) << " " << std::get<1>(t)  << std::endl;
            tot_index_size += sizeof(t);
        }
//        std::cout << "Size of key : " << sizeof(it.first)  << " byte" << "Size of vector : " << sizeof(it.second)  << " byte" << std::endl;
        tot_index_size += sizeof(it.second);

    }

//    // Traversing an unordered map
//    for (auto x : acc_map)
//        std::cout << x.first << " " << x.second << std::endl;

    std::cout << "Total size of index : " << tot_index_size/1000000  << " Mb." << std::endl;

}

static void print_diagnostics2(seq_index2 &h, idx_to_acc &acc_map) {
    uint64_t tot_index_size = 0;

    for (auto &it : h)
    {
//        std::cout << "\nKey: " << it.first << ",  values: ";
        tot_index_size += sizeof(it.first);
        for (auto &t : it.second) // it.second is the vector, i is a tuple
        {
//            std::cout << "(" << std::get<0>(t) << " " << std::get<1>(t)  << " " << std::get<2>(t) << ") ";
            tot_index_size += sizeof(t);
        }
//        std::cout << ". Size of key : " << sizeof(it.first)  << " byte" << "Size of vector : " << sizeof(it.second)  << " byte" << std::endl;
        tot_index_size += sizeof(it.second);

    }

//    // Traversing an unordered map
//    for (auto x : acc_map) {
//        std::cout << x.first << " " << x.second << std::endl;
//    }

    std::cout << "Total size of index : " << tot_index_size/1000000  << " Mb." << std::endl;

}


int main (int argc, char *argv[])
{
//    std::string filename  = "example.txt";
    std::string filename  = "ecoli.fa";
//    std::string filename  = "hg18.fa";
//    std::string choice = "kmer_index";
    std::string choice = "minstrobe_index";
//   std::string choice = "hybridstrobe_index";
//   std::string choice = "randstrobe_index";
    int n = 2;
    int k = 15;
    int w_min = 20;
    int w_max = 40;
    assert(k <= w_min && "k have to be smaller than w_min");
    std::string* file_p;
    file_p = &filename;
    sequences ref_seqs;
    idx_to_acc acc_map;
    read_fasta(ref_seqs, acc_map, filename);

    // Record index creation start time
    auto start = std::chrono::high_resolution_clock::now();

    // CREATE INDEX OF REF SEQUENCES
    if (choice == "kmer_index" ){
        one_pos_index tmp_index; // hash table holding all reference k-mers
        seq_index1 h;
        for (auto x : ref_seqs){
//            generate_kmer_index(h, k, x.second, x.first);
            std::vector<std::tuple<uint64_t, unsigned int, unsigned int>> kmers; // pos, chr_id, kmer hash value
            kmers = seq_to_kmers(k, x.second, x.first);
            tmp_index[x.first] = kmers;
//            print_diagnostics_new(kmers, acc_map);
        }

        std::vector< std::tuple<uint64_t, unsigned int, unsigned int>> mers_vector;
        mers_vector = construct_flat_vector_one_pos(tmp_index); // construct flat array or all kmers
        tmp_index.clear();
        robin_hood::unordered_map< uint64_t, std::tuple<uint64_t, unsigned int >> mers_index; // k-mer -> (offset in flat_vector, occurence count )
        mers_index = index_vector_one_pos(mers_vector); // construct index over flat array
        print_diagnostics_new2(mers_vector, mers_index);
//        print_diagnostics(h, acc_map);
    }

    else if (choice == "hybridstrobe_index" ){
        seq_index2 h;
        if (n == 2 ){
            for (auto x : ref_seqs){
                ;
//                generate_hybridstrobe2_index(h, n, k, w_min, w_max, x.second, x.first);
            }
            print_diagnostics2(h, acc_map);
        }
//        else if (n == 3){ for (auto x : ref_seqs){generate_hybridstrobe3_index(h, k, x.second, x.first);}}
    }
    else if (choice == "minstrobe_index" ){
        if (n == 2 ){
            two_pos_index tmp_index; // hash table holding all reference randstrobes
//            seq_index2 h;
            for (auto x : ref_seqs){
//                generate_minstrobe2_index(h, n, k, w_min, w_max, x.second, x.first);
                std::vector<std::tuple<uint64_t, unsigned int, unsigned int, unsigned int>> minstrobes2; // pos, chr_id, kmer hash value
                minstrobes2 = seq_to_minstrobes2(n, k, w_min, w_max, x.second, x.first);
                tmp_index[x.first] = minstrobes2;
            }
            std::vector< std::tuple<uint64_t, unsigned int, unsigned int, unsigned int>> mers_vector;
            mers_vector = construct_flat_vector_two_pos(tmp_index); // construct flat array or all strobemers
            tmp_index.clear();
            robin_hood::unordered_map< uint64_t, std::tuple<uint64_t, unsigned int >> mers_index; // k-mer -> (offset in flat_vector, occurence count )
            mers_index = index_vector_two_pos(mers_vector); // construct index over flat array
            print_diagnostics_new3(mers_vector, mers_index);
//            print_diagnostics2(h, acc_map);
        }
//        else if (n == 3){ for (auto x : ref_seqs){generate_minstrobe3_index(h, k, x.second, x.first);}}
    }
    else if (choice == "randstrobe_index" ){
        if (n == 2 ){
            two_pos_index tmp_index; // hash table holding all reference randstrobes
//            seq_index2 h;
            for (auto x : ref_seqs){
//                generate_randstrobe2_index(h, n, k, w_min, w_max, x.second, x.first);
                std::vector<std::tuple<uint64_t, unsigned int, unsigned int, unsigned int>> randstrobes2; // pos, chr_id, kmer hash value
                randstrobes2 = seq_to_randstrobes2(n, k, w_min, w_max, x.second, x.first);
                tmp_index[x.first] = randstrobes2;
            }
            std::vector< std::tuple<uint64_t, unsigned int, unsigned int, unsigned int>> mers_vector;
            mers_vector = construct_flat_vector_two_pos(tmp_index); // construct flat array or all strobemers
            tmp_index.clear();
            robin_hood::unordered_map< uint64_t, std::tuple<uint64_t, unsigned int >> mers_index; // k-mer -> (offset in flat_vector, occurence count )
            mers_index = index_vector_two_pos(mers_vector); // construct index over flat array
            print_diagnostics_new3(mers_vector, mers_index);
//            print_diagnostics2(h, acc_map);
        }
//        else if (n == 3){ for (auto x : ref_seqs){generate_randstrobe3_index(h, k, x.second, x.first);}}
    }
    else {
        std::cout << choice << "not implemented : " << std::endl;
    }

    // Record index creation end time
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Total time generating index: " << elapsed.count() << " s\n" <<  std::endl;



}

