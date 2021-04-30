#include <iostream>
#include <fstream>
#include <unordered_map>
#include "robin_hood.h"
#include <vector>
#include <string>


#include "index.hpp"


typedef robin_hood::unordered_map< unsigned int , std::string > sequences;
typedef robin_hood::unordered_map< unsigned int, std::string > idx_to_acc;


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

static void print_diagnostics(seq_index1 &h, idx_to_acc &acc_map) {

    for (auto &it : h)
    {
        std::cout << "\nKey: " << it.first << ",  values: ";
        for (auto &t : it.second) // it.second is the vector, i is a tuple
        {
            std::cout << "tuple: " << std::get<0>(t) << " " << std::get<1>(t)  << std::endl;
        }
        std::cout << "Size of key : " << sizeof(it.first)  << " byte" << "Size of vector : " << sizeof(it.second)  << " byte" << std::endl;

    }

    // Traversing an unordered map
    for (auto x : acc_map)
        std::cout << x.first << " " << x.second << std::endl;

}

static void print_diagnostics2(seq_index2 &h, idx_to_acc &acc_map) {

    for (auto &it : h)
    {
        std::cout << "\nKey: " << it.first << ",  values: ";
        for (auto &t : it.second) // it.second is the vector, i is a tuple
        {
            std::cout << "(" << std::get<0>(t) << " " << std::get<1>(t)  << " " << std::get<2>(t) << ") ";
        }
        std::cout << ". Size of key : " << sizeof(it.first)  << " byte" << "Size of vector : " << sizeof(it.second)  << " byte" << std::endl;

    }

    // Traversing an unordered map
    for (auto x : acc_map)
        std::cout << x.first << " " << x.second << std::endl;

}

//void generate_index(seq_index &h, sequences &ref_seqs, int k)
//{
//    // Traversing an unordered map
//
//}

int main (int argc, char *argv[])
{
    std::string filename  = "example.txt";
//    std::string choice = "kmer_index";
    std::string choice = "minstrobe_index";
//   std::string choice = "hybridstrobe_index";
//   std::string choice = "randstrobe_index";
    int n = 2;
    int k = 10;
    int w_min = 20;
    int w_max = 40;
    assert(k <= w_min && "k have to be smaller than w_min");
    std::string* file_p;
    file_p = &filename;
    sequences ref_seqs;
    idx_to_acc acc_map;
    read_fasta(ref_seqs, acc_map, filename);
//   TODO: create mapping of reference accession to uint index

    // Traversing an unordered map
    for (auto x : ref_seqs) {
        std::cout << x.first << " " << x.second << std::endl;
    }

    // CREATE INDEX OF REF SEQUENCES
    if (choice == "kmer_index" ){
        seq_index1 h;
        for (auto x : ref_seqs){
            generate_kmer_index(h, k, x.second, x.first);
        }
        print_diagnostics(h, acc_map);
    }
    else if (choice == "hybridstrobe_index" ){
        ;
//        if (n == 2 ){ for (auto x : ref_seqs){generate_hybridstrobe2_index(h, k, x.second, x.first);}}
//        else if (n == 3){ for (auto x : ref_seqs){generate_hybridstrobe3_index(h, k, x.second, x.first);}}
    }
    else if (choice == "minstrobe_index" ){
        if (n == 2 ){
            seq_index2 h;
            for (auto x : ref_seqs){
                generate_minstrobe2_index(h, n, k, w_min, w_max, x.second, x.first);
            }
            print_diagnostics2(h, acc_map);
        }
//        else if (n == 3){ for (auto x : ref_seqs){generate_minstrobe3_index(h, k, x.second, x.first);}}
    }
    else if (choice == "randstrobe_index" ){
        if (n == 2 ){
            seq_index2 h;
            for (auto x : ref_seqs){
                generate_randstrobe2_index(h, n, k, w_min, w_max, x.second, x.first);
            }
            print_diagnostics2(h, acc_map);
        }
//        else if (n == 3){ for (auto x : ref_seqs){generate_randstrobe3_index(h, k, x.second, x.first);}}
    }
    else {
        std::cout << choice << "not implemented : " << std::endl;
    }



}

