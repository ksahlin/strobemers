#include <iostream>
#include <fstream>
#include <unordered_map>
#include "robin_hood.h"
#include <vector>
#include <string>
#include <chrono>  // for high_resolution_clock

#include "index.hpp"




typedef robin_hood::unordered_map< std::string , std::string > queries;
typedef robin_hood::unordered_map< unsigned int , std::string > references;
typedef robin_hood::unordered_map< unsigned int, std::string > idx_to_acc;

typedef robin_hood::unordered_map< uint64_t, std::tuple<uint64_t, unsigned int >> vector_index;

typedef std::vector< std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int>> mers_vector;



static void read_references(references &seqs, idx_to_acc &acc_map, std::string fn)
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

//static void print_diagnostics_new(mers_vector_order1  kmers, idx_to_acc &acc_map) {
//    uint64_t tot_index_size = 0;
//    for (int i = 0; i < kmers.size(); ++i)
//    {
//        // access using []
//        auto t = kmers[i];
////        std::cout << "(" << std::get<0>(t) << ", " << std::get<1>(t) << ", " << std::get<2>(t) << ")";
//        tot_index_size += sizeof(t);
////        std::cout << sizeof(t) << std::endl;
//    }
//    std::cout << "Total size of index vector : " << tot_index_size/1000000  << " Mb." << std::endl;
//
//}

//static void print_diagnostics_new2(mers_vector_order1 &mers_vector, vector_index &mers_index ) {
//    uint64_t tot_flat_vector_size = 0;
//    for (int i = 0; i < mers_vector.size(); ++i)
//    {
//        // access using []
//        auto t = mers_vector[i];
////        std::cout << "(" << std::get<0>(t) << ", " << std::get<1>(t) << ", " << std::get<2>(t) << "), ";
//        tot_flat_vector_size += sizeof(t);
//    }
//    std::cout << "Total size of flat kmer vector : " << tot_flat_vector_size/1000000  << " Mb." << std::endl;
//
//    uint64_t tot_hashtable_index_size = 0;
//    for (auto &it : mers_index)
//    {
////        std::cout << it.first << ": (" << std::get<0>(it.second) << ", " << std::get<1>(it.second) << "), " ;
//        tot_hashtable_index_size += sizeof(it.first);
//        tot_hashtable_index_size += sizeof(it.second);
//    }
//    std::cout << "Total size of hash table index : " << tot_hashtable_index_size/1000000  << " Mb." << std::endl;
//}

//static void print_diagnostics_new3(mers_vector_order2 &mers_vector, vector_index &mers_index ) {
//    uint64_t tot_flat_vector_size = 0;
//    for (int i = 0; i < mers_vector.size(); ++i)
//    {
//        // access using []
//        auto t = mers_vector[i];
////        std::cout << "(" << std::get<0>(t) << ", " << std::get<1>(t) << ", " << std::get<2>(t) << ", " << std::get<3>(t) << "), ";
//        tot_flat_vector_size += sizeof(t);
//    }
//    std::cout << "Total size of flat kmer vector : " << tot_flat_vector_size/1000000  << " Mb." << std::endl;
//
//    uint64_t tot_hashtable_index_size = 0;
//    for (auto &it : mers_index)
//    {
////        std::cout << it.first << ": (" << std::get<0>(it.second) << ", " << std::get<1>(it.second) << "), " ;
//        tot_hashtable_index_size += sizeof(it.first);
//        tot_hashtable_index_size += sizeof(it.second);
//    }
//    std::cout << "Total size of hash table index : " << tot_hashtable_index_size/1000000  << " Mb." << std::endl;
//}

//static void print_diagnostics(seq_index1 &h, idx_to_acc &acc_map) {
//    uint64_t tot_index_size = 0;
//    for (auto &it : h)
//    {
////        std::cout << "\nKey: " << it.first << ",  values: ";
//        tot_index_size += sizeof(it.first);
//        for (auto &t : it.second) // it.second is the vector, i is a tuple
//        {
////            std::cout << "tuple: " << std::get<0>(t) << " " << std::get<1>(t)  << std::endl;
//            tot_index_size += sizeof(t);
//        }
////        std::cout << "Size of key : " << sizeof(it.first)  << " byte" << "Size of vector : " << sizeof(it.second)  << " byte" << std::endl;
//        tot_index_size += sizeof(it.second);
//
//    }
//
////    // Traversing an unordered map
////    for (auto x : acc_map)
////        std::cout << x.first << " " << x.second << std::endl;
//
//    std::cout << "Total size of index : " << tot_index_size/1000000  << " Mb." << std::endl;
//
//}


//static void print_diagnostics2(seq_index2 &h, idx_to_acc &acc_map) {
//    uint64_t tot_index_size = 0;
//
//    for (auto &it : h)
//    {
////        std::cout << "\nKey: " << it.first << ",  values: ";
//        tot_index_size += sizeof(it.first);
//        for (auto &t : it.second) // it.second is the vector, i is a tuple
//        {
////            std::cout << "(" << std::get<0>(t) << " " << std::get<1>(t)  << " " << std::get<2>(t) << ") ";
//            tot_index_size += sizeof(t);
//        }
////        std::cout << ". Size of key : " << sizeof(it.first)  << " byte" << "Size of vector : " << sizeof(it.second)  << " byte" << std::endl;
//        tot_index_size += sizeof(it.second);
//
//    }
//
////    // Traversing an unordered map
////    for (auto x : acc_map) {
////        std::cout << x.first << " " << x.second << std::endl;
////    }
//
//    std::cout << "Total size of index : " << tot_index_size/1000000  << " Mb." << std::endl;
//
//}


static void print_diagnostics_new4(mers_vector &mers_vector, vector_index mers_index ) {
    uint64_t tot_flat_vector_size = 0;
    for (int i = 0; i < mers_vector.size(); ++i)
    {
        // access using []
        auto t = mers_vector[i];
//        std::cout << "(" << std::get<0>(t) << ", " << std::get<1>(t) << ", " << std::get<2>(t) << ", " << std::get<3>(t) << ", " << std::get<4>(t) << "), ";
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

//
//mers_vector find_nams(mers_vector &kmers, mers_vector_order2 &mers_vector, vector_index &mers_index){
//
//}

int main (int argc, char *argv[])
{

    ///////////////////// INPUT /////////////////////////

//    std::string filename  = "example2.txt";
    std::string reads_filename  = "example2_reads.txt";
    std::string filename  = "ecoli.fa";
//    std::string filename  = "chr21.fa";
//    std::string choice = "kmers";
//    std::string choice = "minstrobes";
//   std::string choice = "hybridstrobes";
   std::string choice = "randstrobes";
    int n = 3;
    int k = 20;
    int w_min = 21;
    int w_max = 100;
    assert(k <= w_min && "k have to be smaller than w_min");
    std::string* file_p;
    file_p = &filename;
    references ref_seqs;
    idx_to_acc acc_map;
    read_references(ref_seqs, acc_map, filename);

    //////////////////////////////////////////////////////


    //////////// CREATE INDEX OF REF SEQUENCES /////////////////

    // Record index creation start time
    auto start = std::chrono::high_resolution_clock::now();

    three_pos_index tmp_index; // hash table holding all reference mers

    if (choice == "kmers" ){
        for (auto x : ref_seqs){
            mers_vector kmers; //  kmer hash value, pos, chr_id
            kmers = seq_to_kmers(k, x.second, x.first);
            tmp_index[x.first] = kmers;
        }
    }
    else if (choice == "hybridstrobes" ){
        if (n == 2 ){
            for (auto x : ref_seqs){
                mers_vector hybridstrobes2; // pos, chr_id, kmer hash value
                hybridstrobes2 = seq_to_hybridstrobes2(n, k, w_min, w_max, x.second, x.first);
                tmp_index[x.first] = hybridstrobes2;
            }
        }
//        else if (n == 3){ for (auto x : ref_seqs){generate_hybridstrobe3_index(h, k, x.second, x.first);}}
    }
    else if (choice == "minstrobes" ){
        if (n == 2 ){
            for (auto x : ref_seqs){
//                generate_minstrobe2_index(h, n, k, w_min, w_max, x.second, x.first);
                mers_vector minstrobes2; // pos, chr_id, kmer hash value
                minstrobes2 = seq_to_minstrobes2(n, k, w_min, w_max, x.second, x.first);
                tmp_index[x.first] = minstrobes2;
            }
        }
//        else if (n == 3){ for (auto x : ref_seqs){generate_minstrobe3_index(h, k, x.second, x.first);}}
    }
    else if (choice == "randstrobes" ){
        if (n == 2 ){
            for (auto x : ref_seqs){
//                generate_randstrobe2_index(h, n, k, w_min, w_max, x.second, x.first);
                mers_vector randstrobes2; // pos, chr_id, kmer hash value
                randstrobes2 = seq_to_randstrobes2(n, k, w_min, w_max, x.second, x.first);
                tmp_index[x.first] = randstrobes2;
            }
        }
        else if (n == 3){
            for (auto x : ref_seqs){
                mers_vector randstrobes3; // pos, chr_id, kmer hash value
                randstrobes3 = seq_to_randstrobes3(n, k, w_min, w_max, x.second, x.first);
                tmp_index[x.first] = randstrobes3;
            }
        }
    }
    else {
        std::cout << choice << "not implemented : " << std::endl;
    }


    mers_vector all_mers_vector;
    all_mers_vector = construct_flat_vector_three_pos(tmp_index);
    robin_hood::unordered_map< uint64_t, std::tuple<uint64_t, unsigned int >> mers_index; // k-mer -> (offset in flat_vector, occurence count )
    mers_index = index_vector_three_pos(all_mers_vector); // construct index over flat array
    tmp_index.clear();
    print_diagnostics_new4(all_mers_vector, mers_index);

    //////////////////////////////////////////////////////////////////////////

    // Record index creation end time
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Total time generating index: " << elapsed.count() << " s\n" <<  std::endl;



    ///////////////////////////// MAP ///////////////////////////////////////

    // Record matching time
    auto start_map = std::chrono::high_resolution_clock::now();

    std::ifstream query_file(reads_filename);
    std::ofstream output_file;
    output_file.open ("output.tsv");

    std::string line, seq, prev_acc;
    unsigned int q_id = 0;
    while (getline(query_file, line)) {
        if (line[0] == '>') {
            if (seq.length() > 0){
                // generate mers here
                if (choice == "kmers" ){
                    mers_vector kmers; // pos, chr_id, kmer hash value
                    kmers = seq_to_kmers(k, seq, q_id);

                    // Find NAMs
//                    std::vector nams;
//                    nams = find_nams(kmers, all_mers_vector, mers_index);
                    // Output results
//                    output_file << "> " <<  prev_acc << "\n";
//                    output_file << "  " << ref_acc << " " << ref_p << " " << q_pos << " " << "\n";
//                    outfile.write("  {0} {1} {2} {3}\n".format(ref_acc, ref_p, q_pos, k))
                }

            }
            prev_acc = line.substr(1, line.length() -1);
            seq = "";
            q_id ++;
        }
        else {
            seq += line;
        }
    }
    if (seq.length() > 0){
        if (choice == "kmers" ){
            mers_vector kmers; // pos, chr_id, kmer hash value
            kmers = seq_to_kmers(k, seq, q_id);
        }
    }

    query_file.close();
    output_file.close();


    // Record mapping end time
    auto finish_map = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_map = finish_map - start_map;
    std::cout << "Total time generating matches: " << elapsed_map.count() << " s\n" <<  std::endl;


    //////////////////////////////////////////////////////////////////////////


    /////////////////////// FIND AND OUTPUT NAMs ///////////////////////////////





}

