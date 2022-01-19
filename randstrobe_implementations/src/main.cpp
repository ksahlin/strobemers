#include <iostream>
#include <fstream>
#include <unordered_map>
#include "robin_hood.h"
#include <vector>
#include <string>
#include <chrono>  // for high_resolution_clock
#include <cassert>
#include <math.h>
#include <omp.h>

#include <zlib.h>
#include "kseq++.hpp"
#include "index.hpp"

//#define THREAD_NUM 4

using namespace klibpp;

//typedef robin_hood::unordered_map< std::string , std::string > queries;
typedef robin_hood::unordered_map< unsigned int , std::string > references;
typedef robin_hood::unordered_map< unsigned int, std::string > idx_to_acc;
typedef robin_hood::unordered_map< std::string, unsigned int > acc_to_idx;

typedef robin_hood::unordered_map< uint64_t, std::tuple<uint64_t, unsigned int >> vector_index;

typedef std::vector< std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int>> mers_vector;


static inline std::string split_string(std::string str, std::string delimiter = " ")
{
    int start = 0;
    int end = str.find(delimiter);
//    std::cout << str.substr(start, end - start) << std::endl;
    return str.substr(start, end - start);
//    while (end != -1) {
//        std::cout << str.substr(start, end - start) << std::endl;
//        start = end + delimiter.size();
//        end = str.find(delimiter, start);
//    }
//    std::cout << str.substr(start, end - start);

}

static uint64_t read_references(std::vector<std::string> &seqs, std::vector<unsigned int> &lengths, idx_to_acc &acc_map, std::string fn)
{
    uint64_t total_ref_seq_size = 0;
    std::ifstream file(fn);
    std::string line, seq;
    int ref_index = 0;
    while (getline(file, line)) {
        if (line[0] == '>') {
//            std::cout << ref_index << " " << line << std::endl;
            if (seq.length() > 0){
//                seqs[ref_index -1] = seq;
                seqs.push_back(seq);
                lengths.push_back(seq.length());
                total_ref_seq_size += seq.length();
//                std::cout << ref_index - 1 << " here " << seq << " " << seq.length() << " " << seq.size() << std::endl;
//                generate_kmers(h, k, seq, ref_index);
            }
//            acc_map[ref_index] = line.substr(1, line.length() -1); //line;
            acc_map[ref_index] = line.substr(1, line.find(' ') -1 ); //line;
            ref_index++;
            seq = "";
        }
        else {
            seq += line;
        }
    }
    if (seq.length() > 0){
//        seqs[ref_index -1] = seq;
        seqs.push_back(seq);
        lengths.push_back(seq.length());
        total_ref_seq_size += seq.length();
//        std::cout << ref_index -1 << " here2 " << seq << std::endl;
//        generate_kmers(h, k, seq, ref_index);
    }
    file.close();
    return total_ref_seq_size;
}


static inline void print_positions(mers_vector &flat_vector, idx_to_acc &acc_map, std::ofstream &output_file, int hash_func, int link_func ) {
    std::string h_method;
    std::string l_method;

    if (hash_func == 4){
        h_method = "wyhash";
    } else if (hash_func == 3){
        h_method = "xxh64";
    } else if (hash_func == 2){
        h_method = "TW";
    } else if (hash_func == 1){
        h_method = "nohash";
    }

    if (link_func == 1)
    {
        l_method = "Sahlin1";
    } else if (link_func == 2){
        l_method = "Shen";
    } else if (link_func == 3){
        l_method = "Sahlin2";
    } else if (link_func == 4){
        l_method = "Guo-Pibri";
    } else if (link_func == 5){
        l_method = "Liu-Patro-Li";
    } else if (link_func == 6){
        l_method = "Liu-Patro-Li_wyhash";
    }


    for (size_t i = 0; i < flat_vector.size(); ++i)
    {
        // access using []
        auto t = flat_vector[i];
//        < std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int>
//        (hash_randstrobe3, ref_index, seq_pos_strobe1, seq_pos_strobe2, seq_pos_strobe3)
//        auto hash_v = std::get<0>(t);
        auto ref_index = std::get<1>(t);
        auto seq_pos_strobe1 = std::get<2>(t);
        auto seq_pos_strobe2 = std::get<3>(t);
        auto seq_pos_strobe3 = std::get<4>(t);
        output_file << h_method << "," << l_method << "," <<  acc_map[ref_index]  << "," << seq_pos_strobe1 << "," << seq_pos_strobe3 << "\n";
    }

}



static inline std::vector<nam> find_nams_unique(mers_vector &query_mers, mers_vector &mers_vector, vector_index &mers_index, int k){
    std::cout << "ENTER FIND NAMS UNIQUE " <<  std::endl;
    robin_hood::unordered_map< unsigned int, std::vector<hit>> hits_per_ref; // [ref_id] -> vector( struct hit)
    uint64_t hit_count_reduced = 0;
    uint64_t hit_count_all = 0;
    uint64_t total_mers = 0;
    for (auto &q : query_mers)
    {
        hit h;
        h.query_s = std::get<2>(q);
        h.query_e = std::get<4>(q) + k;
        total_mers ++;
//        std::cout << h.query_s << " " << h.query_e <<  std::endl;

        uint64_t mer_hashv = std::get<0>(q);
        if (mers_index.find(mer_hashv) != mers_index.end()){ //  In  index
            std::tuple<uint64_t, unsigned int> mer;
            mer = mers_index[mer_hashv];
            uint64_t offset = std::get<0>(mer);
            unsigned int count = std::get<1>(mer);
            if (count == 1){
                auto r = mers_vector[offset];
                unsigned int ref_s = std::get<2>(r);
                unsigned int ref_e = std::get<4>(r) + k;
                unsigned int ref_id = std::get<1>(r);

                h.ref_s = ref_s;
                h.ref_e = ref_e;
                hits_per_ref[ref_id].push_back(h);

                hit_count_all ++;
            }
        }
    }

//    std::cout << "NUMBER OF HITS GENERATED: " << hit_count_all << std::endl;
//    std::cout << "TOTAL STROBEMERS GENERATED: " << total_mers << std::endl;

//    std::cout << "NUMBER OF REDUCED HITS GENERATED: " << hit_count_reduced << std::endl;
    std::vector<nam> open_nams;
    std::vector<nam> final_nams; // [ref_id] -> vector(struct nam)

    for (auto &it : hits_per_ref)
    {
        unsigned int ref_id = it.first;
        std::vector<hit> hits = it.second;
        open_nams = std::vector<nam> (); // Initialize vector
        uint64_t prev_q_start = 0;
        for (auto &h : hits){
            bool is_added = false;
            for (auto & o : open_nams) {

                // Extend NAM
                if ( ( o.previous_query_start < h.query_s) && (h.query_s <= o.query_e ) && ( o.previous_ref_start <= h.ref_s) && (h.ref_s <= o.ref_e) ){ // && (o.previous_ref_start <= h.ref_s)  && (o.previous_query_start <= h.query_s) && (hit_copy_id <= o.copy_id)
                    if (h.query_e > o.query_e) {
                        o.query_e = h.query_e;
                    }
                    if (h.ref_e > o.ref_e) {
                        o.ref_e = h.ref_e;
                    }
                    o.previous_query_start = h.query_s;
                    o.previous_ref_start = h.ref_s; // keeping track so that we don't . Can be caused by interleaved repeats.
                    is_added = true;
                    break;
                }

            }

            // Add the hit to open matches
            if (not is_added){
                nam n;
                n.query_s = h.query_s;
                n.query_e = h.query_e;
                n.ref_s = h.ref_s;
                n.ref_e = h.ref_e;
                n.ref_id = ref_id;
                n.previous_query_start = h.query_s;
                n.previous_ref_start = h.ref_s;
//                n.copy_id = hit_copy_id;
                open_nams.push_back(n);
            }


            // Only filter if we have advanced at least k nucleotides
            if (h.query_s > prev_q_start + k) {

                // Output all NAMs from open_matches to final_nams that the current hit have passed
                for (auto &n : open_nams) {
                    if (n.query_e < h.query_s) {
                        final_nams.push_back(n);
                    }
                }

                // Remove all NAMs from open_matches that the current hit have passed
                unsigned int c = h.query_s;
                auto predicate = [c](decltype(open_nams)::value_type const &nam) { return nam.query_e < c; };
                open_nams.erase(std::remove_if(open_nams.begin(), open_nams.end(), predicate), open_nams.end());
                prev_q_start = h.query_s;
            }


        }

        // Add all current open_matches to final NAMs
        for (auto &n : open_nams){
//        for (size_t i = 0; i < open_nams.size(); ++i){
            final_nams.push_back(n);
        }
    }

//    for (auto &n : final_nams){
//        std::cout << n.ref_id << ": (" << n.query_s << ", " << n.query_e << ", " << n.ref_s << ", " << n.ref_e << ")" << std::endl;
//    }


    return final_nams;
}


static inline std::vector<nam> find_nams(mers_vector &query_mers, mers_vector &mers_vector, vector_index &mers_index, int k, acc_to_idx &acc_map_to_idx,  std::string query_acc, unsigned int filter_cutoff){
//    std::cout << "ENTER FIND NAMS " <<  std::endl;
    robin_hood::unordered_map< unsigned int, std::vector<hit>> hits_per_ref; // [ref_id] -> vector( struct hit)
    uint64_t hit_count_reduced = 0;
    uint64_t hit_count_all = 0;
    uint64_t total_mers = 0;
    for (auto &q : query_mers)
//    for (size_t i = 0; i < query_mers.size(); ++i)
    {
        hit h;
        h.query_s = std::get<2>(q);
        h.query_e = std::get<4>(q) + k;
        total_mers ++;
//        std::cout << h.query_s << " " << h.query_e <<  std::endl;

        uint64_t mer_hashv = std::get<0>(q);
        if (mers_index.find(mer_hashv) != mers_index.end()){ //  In  index
//            std::cout << "Found: " << h.query_s <<  std::endl;
            std::tuple<uint64_t, unsigned int> mer;
            mer = mers_index[mer_hashv];
            uint64_t offset = std::get<0>(mer);
            unsigned int count = std::get<1>(mer);
            if (count <= filter_cutoff){
                for(size_t j = offset; j < offset+count; ++j)
                {

                    auto r = mers_vector[j];
                    unsigned int ref_s = std::get<2>(r);
                    unsigned int ref_e = std::get<4>(r) + k;
                    unsigned int ref_id = std::get<1>(r);
                    if ( (acc_map_to_idx.find(query_acc) != acc_map_to_idx.end()) && ref_id ==  acc_map_to_idx[query_acc]) { // No self mappings
                        continue;
                    }
                    h.ref_s = ref_s;
                    h.ref_e = ref_e;
                    hits_per_ref[ref_id].push_back(h);

                    hit_count_all ++;

                }
            }
//            else {
//                std::cout << "Strobemer had count: " <<  count << std::endl;
//            }

        }
    }

//    std::cout << "NUMBER OF HITS GENERATED: " << hit_count_all << " , unique refs: " << hits_per_ref.size() << std::endl;
//    std::cout << "TOTAL STROBEMERS GENERATED: " << total_mers << std::endl;

//    std::cout << "NUMBER OF REDUCED HITS GENERATED: " << hit_count_reduced << std::endl;
    std::vector<nam> open_nams;
    std::vector<nam> final_nams; // [ref_id] -> vector(struct nam)

    for (auto &it : hits_per_ref)
    {
        unsigned int ref_id = it.first;
        std::vector<hit> hits = it.second;
        open_nams = std::vector<nam> (); // Initialize vector
        uint64_t prev_q_start = 0;
        for (auto &h : hits){
            bool is_added = false;
                for (auto & o : open_nams) {

                    // Extend NAM
                    if ( ( o.previous_query_start < h.query_s) && (h.query_s <= o.query_e ) && ( o.previous_ref_start <= h.ref_s) && (h.ref_s <= o.ref_e) ){ // && (o.previous_ref_start <= h.ref_s)  && (o.previous_query_start <= h.query_s) && (hit_copy_id <= o.copy_id)

                        if (h.query_e > o.query_e) {
                            o.query_e = h.query_e;
                        }
                        if (h.ref_e > o.ref_e) {
                            o.ref_e = h.ref_e;
                        }
                        o.previous_query_start = h.query_s;
                        o.previous_ref_start = h.ref_s; // keeping track so that we don't . Can be caused by interleaved repeats.
                        o.n_hits ++;
                        is_added = true;
                        break;
                    }

                }

            // Add the hit to open matches
            if (not is_added){
                nam n;
                n.query_s = h.query_s;
                n.query_e = h.query_e;
                n.ref_s = h.ref_s;
                n.ref_e = h.ref_e;
                n.ref_id = ref_id;
                n.previous_query_start = h.query_s;
                n.previous_ref_start = h.ref_s;
                n.n_hits = 1;
//                n.copy_id = hit_copy_id;
                open_nams.push_back(n);
            }

            // Only filter if we have advanced at least k nucleotides
            if (h.query_s > prev_q_start + k) {

                // Output all NAMs from open_matches to final_nams that the current hit have passed
                for (auto &n : open_nams) {
                    if (n.query_e < h.query_s) {
                        final_nams.push_back(n);
                    }
                }

                // Remove all NAMs from open_matches that the current hit have passed
                unsigned int c = h.query_s;
                auto predicate = [c](decltype(open_nams)::value_type const &nam) { return nam.query_e < c; };
                open_nams.erase(std::remove_if(open_nams.begin(), open_nams.end(), predicate), open_nams.end());
                prev_q_start = h.query_s;
            }


        }

        // Add all current open_matches to final NAMs
        for (auto &n : open_nams){
//        for (size_t i = 0; i < open_nams.size(); ++i){
            final_nams.push_back(n);
        }
    }

//    for (auto &n : final_nams){
//        std::cout << n.ref_id << ": (" << n.query_s << ", " << n.query_e << ", " << n.ref_s << ", " << n.ref_e << ")" << std::endl;
//    }


    return final_nams;
}

static inline std::string reverse_complement(std::string &read) {
    auto read_rev = read;
    std::reverse(read_rev.begin(), read_rev.end()); // reverse
//    std::cout << read_rev << std::endl;
    for (size_t j = 0; j < read_rev.length(); ++j) { // complement
        if (read_rev[j] == 'A') read_rev[j] = 'T';
        else if (read_rev[j] == 'T') read_rev[j] = 'A';
        else if (read_rev[j] == 'C') read_rev[j] = 'G';
        else if (read_rev[j] == 'G') read_rev[j] = 'C';
    }
    return read_rev;
}

static inline bool compareByQueryCoord(const nam &a, const nam &b)
{
    // first sort on ref ID, then on query, then on reference
    return (a.ref_id < b.ref_id) ||
           ( (a.ref_id == b.ref_id) && (a.query_s < b.query_s) ) ||
           ((a.ref_id == b.ref_id) && (a.query_s == b.query_s ) && (a.ref_s < b.ref_s)) ;
}

static inline bool compareByQueryLength(const nam &a, const nam &b)
{
    return (a.query_e - a.query_s) < ( b.query_e - b.query_s);
}

static inline bool score(const nam &a, const nam &b)
{
    return ( (a.n_hits * (a.query_e - a.query_s)) > (b.n_hits * (b.query_e - b.query_s)) );
}

static inline void output_nams(std::vector<nam> &nams, std::ofstream &output_file, std::string query_acc, idx_to_acc &acc_map, bool is_rc) {
    // Output results
    if (is_rc) {
        output_file << "> " << query_acc << " Reverse\n";
    }
    else{
        output_file << "> " << query_acc << "\n";
    }
    for (auto &n : nams) {
        output_file << "  " << acc_map[n.ref_id]  << " " << n.ref_s + 1 << " " << n.query_s + 1 << " " << n.ref_e - n.ref_s << "\n";
//      python: outfile.write("  {0} {1} {2} {3}\n".format(ref_acc, ref_p, q_pos, k))
    }
}

void print_usage() {
    std::cerr << "\n";
    std::cerr << "Randstrobe evaluation\n";
    std::cerr << "\n";
    std::cerr << "StrobeMap [options] <references.fasta> <queries.fast[a/q]>\n";
    std::cerr << "options:\n";
    std::cerr << "\t-n INT number of strobes [2 or 3]\n";
    std::cerr << "\t-k INT strobe length, limited to 32 [20]\n";
    std::cerr << "\t-v INT strobe w_min offset [21]\n";
    std::cerr << "\t-w INT strobe w_max offset [120]\n";
    std::cerr << "\t-t INT number of threads [3]\n";
    std::cerr << "\t-o name of output tsv-file [output.tsv]\n";
    std::cerr << "\t-x Choice of hash function to use; 1: nohash, 2: wanghash, 3:xxhash [1]. \n";
    std::cerr << "\t-l Choice of link function to use; 1: method1 (sahlin modulo), 2: method2 (shen bitwise AND), 3: method3 (guo_pibri XOR), 4: method4 (sahlin bitcount XOR), 5: method5 (Liu-Patro-Li, concatenation), 6: method6 (Liu-Patro-Li, concatenation using wyhash for linking). \n";
//    std::cerr << "\t-C UINT Mask (do not process) strobemer hits with count larger than C [1000]\n";
//    std::cerr << "\t-L UINT Print at most L NAMs per query [1000]. Will print the NAMs with highest score S = n_strobemer_hits * query_span. \n";
//    std::cerr << "\t-S Sort output NAMs for each query based on score. Default is to sort first by ref ID, then by query coordinate, then by reference coordinate. \n";
//    std::cerr << "\t-s Split output into one file per thread and forward/reverse complement mappings. \n\t   This option is used to generate format compatible with uLTRA long-read RNA aligner and requires \n\t   option -o to be specified as a folder path to uLTRA output directory, e.g., -o /my/path/to/uLTRA_output/ \n";
}


// g++ -std=c++11 main.cpp index.cpp -o StrobeMap -O3 -mavx2

int main (int argc, char *argv[])
{


    if (argc < 3) {
        print_usage();
        return 0;
    }

    // Default parameters
    int hash_func = 1;
    int link_func = 1;
    int n = 2;
    int k = 20;
    float f = 0.0002;
    std::string output_file_hits = "hits.tsv"; // q_position_start, q_position_end, r_position_start, r_position_end
    int w_min = 21;
    int w_max = 100;
    int n_threads = 3;
    bool unique = false;
    bool output_specified = false;
    bool wmin_set = false;
    unsigned int filter_cutoff = 1000;
    unsigned int max_lines = 1000;
    bool sort_on_score_set = false;
    int opn = 1;
    while (opn < argc) {
        bool flag = false;
        if (argv[opn][0] == '-') {
            if (argv[opn][1] == 'n') {
                n = std::stoi(argv[opn + 1]);
                opn += 2;
                flag = true;
            } else if (argv[opn][1] == 'k') {
                k = std::stoi(argv[opn + 1]);
                opn += 2;
                flag = true;
            } else if (argv[opn][1] == 'o') {
                output_file_hits = argv[opn + 1];
                opn += 2;
                flag = true;
                output_specified = true;
            } else if (argv[opn][1] == 'v') {
                w_min = std::stoi(argv[opn + 1]);
                opn += 2;
                flag = true;
                wmin_set = true;
            } else if (argv[opn][1] == 'w') {
                w_max = std::stoi(argv[opn + 1]);
                opn += 2;
                flag = true;
            } else if (argv[opn][1] == 'x') {
                hash_func = std::stoi(argv[opn + 1]);
                opn += 2;
                flag = true;
            } else if (argv[opn][1] == 'l') {
                link_func =std::stoi( argv[opn + 1]);
                opn += 2;
                flag = true;
            } else if (argv[opn][1] == 'u') {
                unique = true;
                opn += 1;
                flag = true;
            } else if (argv[opn][1] == 't') {
                n_threads = std::stoi(argv[opn + 1]);
                opn += 2;
                flag = true;
            } else if (argv[opn][1] == 'C') {
                filter_cutoff = std::stoi(argv[opn + 1]);
                opn += 2;
                flag = true;
            } else if (argv[opn][1] == 'L') {
                max_lines = std::stoi(argv[opn + 1]);
                opn += 2;
                flag = true;
            } else if (argv[opn][1] == 'S') {
                sort_on_score_set = true;
                opn += 1;
                flag = true;
            }
            else {
                print_usage();
            }
        }
        if (!flag)
            break;
    }

    if (!wmin_set){
        w_min = k+1; // Update default w_min to k + 1 if user has specified non-default k and not set w_min parameter
    }
    omp_set_num_threads(n_threads); // set number of threads in "parallel" blocks
    std::cout << "Using" << std::endl;
    std::cout << "n: " << n << std::endl;
    std::cout << "k: " << k << std::endl;
    std::cout << "w_min: " << w_min << std::endl;
    std::cout << "w_max: " << w_max << std::endl;
    std::cout << "t: " << n_threads << std::endl;
    std::cout << "C: " << filter_cutoff << std::endl;
    std::cout << "L: " << max_lines << std::endl;

//    assert(k <= (w/2)*w_min && "k should be smaller than (w/2)*w_min to avoid creating short strobemers");
    assert(k > 7 && "You should really not use too small strobe size!");
    int filter_nams = 0;
//    assert(k <= w_min && "k have to be smaller than w_min");
    assert(k <= 32 && "k have to be smaller than 32!");
    // File name to reference
    std::string filename = argv[opn];
    opn++;
    const char *reads_filename = argv[opn];

    auto start = std::chrono::high_resolution_clock::now();
    std::vector<std::string> ref_seqs;
    std::vector<unsigned int> ref_lengths;
    uint64_t total_ref_seq_size;
    idx_to_acc acc_map;

    total_ref_seq_size = read_references(ref_seqs, ref_lengths, acc_map, filename);

    //////////////////////////////////////////////////////


    //////////// CREATE INDEX OF REF SEQUENCES /////////////////

    // Record index creation start time
    auto start_flat_vector = std::chrono::high_resolution_clock::now();
    mers_vector flat_vector;
    flat_vector.reserve(total_ref_seq_size);
    auto start_generating_randstrobes = std::chrono::high_resolution_clock::now();
    unsigned int mer_cnt = 0;

    std::cout << "total_ref_seq_size " << total_ref_seq_size << std::endl;
    std::cout << "Number of refs: " << ref_seqs.size() << std::endl;
    std::chrono::duration<double> elapsed_hash;
    std::chrono::duration<double> elapsed_link;

    std::vector<uint64_t> string_hashes;
    std::vector<unsigned int> pos_to_seq_choord;

    if (link_func == 5) { // method 5 requires combined link and hash
        hash_func = 3;
    } else if (link_func == 6){
        hash_func = 4;
    }

    if (n == 2 ) {
        mers_vector randstrobes2; // pos, chr_id, kmer hash value
        for (size_t i = 0; i < ref_seqs.size(); ++i) {

            string_hashes.reserve(ref_lengths[i]);
            pos_to_seq_choord.reserve(ref_lengths[i]);

            if (link_func == 5) { // method 5 requires combined link and hash
                auto start_hash = std::chrono::high_resolution_clock::now();
                string_to_hash_nohash(ref_seqs[i], string_hashes, pos_to_seq_choord, k);
                auto end_hash = std::chrono::high_resolution_clock::now();
                elapsed_hash += end_hash - start_hash;

                auto start_link = std::chrono::high_resolution_clock::now();
                randstrobes2 = link_2_strobes_liu_patro_li(w_min, w_max, string_hashes, pos_to_seq_choord, i);
                auto end_link = std::chrono::high_resolution_clock::now();
                elapsed_link += end_link - start_link;
                for (auto &t : randstrobes2) {
                    flat_vector.push_back(t);
                }
                string_hashes.clear();
                randstrobes2.clear();
                pos_to_seq_choord.clear();
            } else if (link_func == 6) { // method 6 requires combined link and hash
                auto start_hash = std::chrono::high_resolution_clock::now();
                string_to_hash_nohash(ref_seqs[i], string_hashes, pos_to_seq_choord, k);
                auto end_hash = std::chrono::high_resolution_clock::now();
                elapsed_hash += end_hash - start_hash;

                auto start_link = std::chrono::high_resolution_clock::now();
                randstrobes2 = link_2_strobes_liu_patro_li_wyhash(w_min, w_max, string_hashes, pos_to_seq_choord, i);
                auto end_link = std::chrono::high_resolution_clock::now();
                elapsed_link += end_link - start_link;
                for (auto &t : randstrobes2) {
                    flat_vector.push_back(t);
                }
                string_hashes.clear();
                randstrobes2.clear();
                pos_to_seq_choord.clear();
            } else { // the other methods can be separated
                auto start_hash = std::chrono::high_resolution_clock::now();
                // first hash
                if (hash_func == 1) {
                    string_to_hash_nohash(ref_seqs[i], string_hashes, pos_to_seq_choord, k);

                } else if (hash_func == 2) {
                    string_to_hash_wang(ref_seqs[i], string_hashes, pos_to_seq_choord, k);

                } else if (hash_func == 3) {
                    string_to_hash_xxhash(ref_seqs[i], string_hashes, pos_to_seq_choord, k);
                }
                auto end_hash = std::chrono::high_resolution_clock::now();
                elapsed_hash += end_hash - start_hash;


                // then link
                auto start_link = std::chrono::high_resolution_clock::now();
                if (link_func == 1) {
                    randstrobes2 = link_2_strobes_sahlin1(w_min, w_max, string_hashes, pos_to_seq_choord, i);
                } else if (link_func == 2) {
                    randstrobes2 = link_2_strobes_shen(w_min, w_max, string_hashes, pos_to_seq_choord, i);
                } else if (link_func == 3) {
                    randstrobes2 = link_2_strobes_sahlin2(w_min, w_max, string_hashes, pos_to_seq_choord, i);
                } else if (link_func == 4) {
                    randstrobes2 = link_2_strobes_guo_pibri(w_min, w_max, string_hashes, pos_to_seq_choord, i);
                }
                auto end_link = std::chrono::high_resolution_clock::now();
                elapsed_link += end_link - start_link;
                for (auto &t : randstrobes2) {
                    flat_vector.push_back(t);
                }
                string_hashes.clear();
                randstrobes2.clear();
                pos_to_seq_choord.clear();
            }

//            std::cout << "Done with ref: " << i << std::endl;
        }
    }
//    else if (n == 3){
//        for (size_t i = 0; i < ref_seqs.size(); ++i) {
//            // first hash
//            string_hashes.reserve(ref_lengths[i]);
//            pos_to_seq_choord.reserve(ref_lengths[i]);
//            string_to_hash_wang(ref_seqs[i], string_hashes, pos_to_seq_choord, k);
//
//            // then link
//            mers_vector randstrobes3; // pos, chr_id, kmer hash value
//            randstrobes3 = link_3_strobes_method2(w_min, w_max, string_hashes, pos_to_seq_choord, i);
//            for (auto &t : randstrobes3) {
//                flat_vector.push_back(t);
//            }
//        }
//    }


//    float elapsed_hash_rounded = truncf(elapsed_hash.count() * 10) / 10;


    std::cout << "Total time hashing: " << elapsed_hash.count() << " s\n" <<  std::endl;
    std::cout << "Total time linking: " << elapsed_link.count() << " s\n" <<  std::endl;


    auto finish_generating_randstrobes = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_mers = finish_generating_randstrobes - start_generating_randstrobes;
//    float rounded = truncf(elapsed_mers.count() * 10) / 10;
//    std::cout << rounded << "s";
    std::cout << "Total time generating randstrobes (hashing + linking): " << elapsed_mers.count() << " s\n" <<  std::endl;

    std::cout << "Ref vector actual size: " << flat_vector.size() << std::endl;
    flat_vector.shrink_to_fit();

    std::string output_filename_sampled = "positions.tsv"; // q_position_end (last position) if three strobes
    std::ofstream output_file_sampled; // q_position_end (last position) if three strobes;
    output_file_sampled.open(output_filename_sampled);
    print_positions(flat_vector, acc_map, output_file_sampled, hash_func, link_func ); // TODO
    return 0;




    uint64_t unique_mers = 0;
    process_flat_vector(flat_vector, unique_mers);
    std::cout << "Unique strobemers: " << unique_mers  <<  std::endl;
    auto finish_flat_vector = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_flat_vector = finish_flat_vector - start_flat_vector;
    std::cout << "Total time generating flat vector: " << elapsed_flat_vector.count() << " s\n" <<  std::endl;

    auto start_hash_index = std::chrono::high_resolution_clock::now();
    kmer_lookup mers_index; // k-mer -> (offset in flat_vector, occurence count )
    mers_index.reserve(unique_mers);
    filter_cutoff = index_vector(flat_vector, mers_index, f); // construct index over flat array
    auto finish_hash_index = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_hash_index = finish_hash_index - start_hash_index;
    std::cout << "Total time generating hash table index: " << elapsed_hash_index.count() << " s\n" <<  std::endl;

    //////////////////////////////////////////////////////////////////////////

    // Record index creation end time
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Total time indexing: " << elapsed.count() << " s\n" <<  std::endl;



//    ///////////////////////////// MAP ///////////////////////////////////////
//
//    // Record matching time
//    auto start_map = std::chrono::high_resolution_clock::now();
//
//    //    std::ifstream query_file(reads_filename);
////    KSeq record;
//    gzFile fp = gzopen(reads_filename, "r");
////    auto ks = make_kstream(fp, gzread, mode::in);
//    auto ks = make_ikstream(fp, gzread);
//    int n_q_chunk_size = 500000;
//
////    std::string line, seq, prev_acc;
//    std::string seq_rc;
//    std::string acc = "";
//    unsigned int q_id = 0;
//    unsigned int read_cnt = 0;
//    unsigned int cut_nam_vec_at;
//    unsigned int cut_nam_rc_vec_at;
//    mers_vector query_mers; // pos, chr_id, kmer hash value
//    mers_vector query_mers_rc; // pos, chr_id, kmer hash value
//
//    std::ofstream output_file;
//    output_file.open(output_file_name);
//    while (ks ) {
////        ks >> record;
//         auto records = ks.read(n_q_chunk_size);  // read a chunk of 500000 records
//        int n_it =  records.size();
//         std::cout << "Mapping chunk of " << n_it << " query sequences... " << std::endl;
//        #pragma omp parallel for num_threads(n_threads) shared(read_cnt, output_file, q_id) private(acc,seq_rc, query_mers,query_mers_rc)
//        for(int i = 0; i < n_it; ++i){
//            auto record =records[i];
////            for (auto & record : records){
//            read_cnt ++;
//            acc = split_string(record.name);
//    //        std::cout << acc << std::endl;
//    //        if (!record.comment.empty()) std::cout << record.comment << std::endl;
//    //        std::cout << record.seq << std::endl;
//    //        if (!record.qual.empty()) std::cout << record.qual << std::endl;
//
//            if (n == 2 ){
//                query_mers = seq_to_randstrobes2(n, k, w_min, w_max, record.seq, q_id);
//                seq_rc = reverse_complement(record.seq);
//                query_mers_rc = seq_to_randstrobes2(n, k, w_min, w_max, seq_rc, q_id);
//            }
//
//            else if (n == 3){
//                query_mers = seq_to_randstrobes3(n, k, w_min, w_max, record.seq, q_id);
//                seq_rc = reverse_complement(record.seq);
//                query_mers_rc = seq_to_randstrobes3(n, k, w_min, w_max, seq_rc, q_id);
//            }
//
//            // Find NAMs
//            std::vector<nam> nams; // (r_id, r_pos_start, r_pos_end, q_pos_start, q_pos_end)
//            std::vector<nam> nams_rc; // (r_id, r_pos_start, r_pos_end, q_pos_start, q_pos_end)
//
//            if (unique) {
//                nams = find_nams_unique(query_mers, all_mers_vector, mers_index, k);
//                nams_rc = find_nams_unique(query_mers_rc, all_mers_vector, mers_index, k);
//            }
//            else {
//                nams = find_nams(query_mers, all_mers_vector, mers_index, k, acc_map, acc, filter_cutoff);
//                nams_rc = find_nams(query_mers_rc, all_mers_vector, mers_index, k, acc_map, acc, filter_cutoff);
//            }
//
//            // Sort on score
//            std::sort(nams.begin(), nams.end(), score);
//            std::sort(nams_rc.begin(), nams_rc.end(), score);
//
//            // Take first L NAMs for output
//            cut_nam_vec_at = (max_lines < nams.size()) ? max_lines : nams.size();
//            std::vector<nam> nams_cut(nams.begin(), nams.begin() + cut_nam_vec_at);
//            cut_nam_rc_vec_at = (max_lines < nams_rc.size()) ? max_lines : nams_rc.size();
//            std::vector<nam> nams_rc_cut(nams_rc.begin(), nams_rc.begin() + cut_nam_rc_vec_at);
//
//            //Sort hits based on start choordinate on query sequence
//            if (!sort_on_score_set) {
//                std::sort(nams_cut.begin(), nams_cut.end(), compareByQueryCoord);
//                std::sort(nams_rc_cut.begin(), nams_rc_cut.end(), compareByQueryCoord);
//            }
//
//            // Output results
//            #pragma omp critical (datawrite)
//            {
//                output_nams(nams, output_file, acc, acc_map, false);
//                output_nams(nams_rc, output_file, acc, acc_map, true);
//            };
//    //        std::cout << "Processed " << read_cnt << "reads. " << std::endl;
//
////            if (read_cnt % 10000 == 0){
////                std::cout << "Processed " << read_cnt << "reads. " << std::endl;
////            }
//            q_id ++;
//        }
//    }
//    output_file.close();
//
//    gzclose(fp);
//
//
//
//    // Record mapping end time
//    auto finish_map = std::chrono::high_resolution_clock::now();
//    std::chrono::duration<double> elapsed_map = finish_map - start_map;
//    std::cout << "Total time mapping: " << elapsed_map.count() << " s\n" <<  std::endl;

}

