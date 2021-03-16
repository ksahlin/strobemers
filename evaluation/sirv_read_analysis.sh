#!/bin/bash

# RUN scripts e.g. as:   ./sirv_match_lengths.sh /Users/kxs624/Documents/workspace/strobemers/  /Users/kxs624/tmp/STROBEMERS/sirv_match_lengths/  /Users/kxs624/Documents/data/ont/lc19_pcs109_subsample_full_length_pychopper2_phmmer.fq  100

inbase=$1
outbase=$2
original_reads=$3
nr_reads=$4

experiment_dir=$inbase"/evaluation/"

mkdir -p $outbase


IFS=$'\n'       # make newlines the only separator
# set -f          # disable globbing


# ########### READS TO REFERENCES #######################
# #######################################################

# normalized_hit_length_file=$outbase/"normalized_hit_length.csv"
# nr_hits_file=$outbase/"nr_hits.csv"
# coverage_file=$outbase/"coverage.csv"


# # results_file has format: method,read_acc,ref_id,ref_len,match_length,normalized_match_length,
# echo -n  "method","read_id","ref_id","ref_length","match_length","normalized_match_length"$'\n' > $normalized_hit_length_file
# echo -n  "method","read_id","ref_id","ref_length","nr_hits"$'\n' > $nr_hits_file
# echo -n  "method","read_id","ref_id","ref_length","coverage"$'\n' > $coverage_file

# # align original reads with minimap2 
# # original_reads_mapped=$outbase/original_reads.sam
# # minimap2 -a --eqx -k 10 -w 1 $inbase/test_data/sirv_transcripts.fasta  $original_reads  > $original_reads_mapped
# original_reads_mapped=/Users/kxs624/Documents/data/ont/sirv/cDNA/lc19_pcs109_subsample_full_length_pychopper2_phmmer.sam  #premapped 

# # Subsample nr_reads reads from each transcript
# echo python $experiment_dir/sample_reads.py $original_reads $inbase/data/sirv_transcripts.fasta $original_reads_mapped $outbase/fastq $nr_reads
# # python $experiment_dir/sample_reads.py $original_reads $inbase/data/sirv_transcripts.fasta $original_reads_mapped $outbase/fastq $nr_reads
# ###############

# mkdir -p $outbase/results/strobmers/
# mkdir -p $outbase/refs/

# FILES=$outbase/fastq/*
# for f in $FILES
# do
#   f_base="$(basename -- $f .fastq)"
#   echo "Processing $f_base file..."

#   # run strobmers
#   python strobe_match.py --queries $f --references $outbase/refs/$f_base"_ref.fastq" --outfolder $outbase/results/strobmers/ --prefix $f_base --k 15 --w 10
#   python $experiment_dir/print_hit_statistics.py $outbase/results/strobmers/$f_base.tsv --refs $outbase/refs/$f_base"_ref.fastq" --outfolder $outbase/ --method strobemers  

#   #run kmers
#   python strobe_match.py --queries $f --references $outbase/refs/$f_base"_ref.fastq" --outfolder $outbase/results/kmers/ --prefix $f_base --k 30 --kmer_index --w 10
#   python $experiment_dir/print_hit_statistics.py $outbase/results/kmers/$f_base.tsv --refs $outbase/refs/$f_base"_ref.fastq" --outfolder $outbase/ --method kmers 
# done


# # Plot1: Sort plot sirvs by length on x-axis! normalized_hit_length.csv
# # Plot2: Plot total coverage of hits  covarage.csv
# # Plot3 number of hits per read nr_hits.csv

# python $experiment_dir/plots.py $normalized_hit_length_file $coverage_file $nr_hits_file $plot_file".pdf"

# ########### READS TO REFERENCES #######################
# #######################################################



########### READS TO READS ############################
#######################################################

normalized_hit_length_file=$outbase/"normalized_hit_length_read_vs_read.csv"
nr_hits_file=$outbase/"nr_hits_read_vs_read.csv"
coverage_file=$outbase/"coverage_read_vs_read.csv"


# results_file has format: method,read_acc,ref_id,ref_len,match_length,normalized_match_length,
echo -n  "method","read_id","ref_id","ref_length","match_length","normalized_match_length"$'\n' > $normalized_hit_length_file
echo -n  "method","read_id","ref_id","ref_length","nr_hits"$'\n' > $nr_hits_file
echo -n  "method","read_id","ref_id","ref_length","coverage"$'\n' > $coverage_file

# align original reads with minimap2 
# original_reads_mapped=$outbase/original_reads.sam
# minimap2 -a --eqx -k 10 -w 1 $inbase/test_data/sirv_transcripts.fasta  $original_reads  > $original_reads_mapped
original_reads_mapped=/Users/kxs624/Documents/data/ont/sirv/cDNA/lc19_pcs109_subsample_full_length_pychopper2_phmmer.sam  #premapped 

# Subsample nr_reads reads from each transcript
echo python $experiment_dir/sample_reads.py $original_reads $inbase/data/sirv_transcripts.fasta $original_reads_mapped $outbase/fastq $nr_reads
# python $experiment_dir/sample_reads.py $original_reads $inbase/data/sirv_transcripts.fasta $original_reads_mapped $outbase/fastq $nr_reads
###############

mkdir -p $outbase/results/strobmers/
mkdir -p $outbase/refs/

FILES=$outbase/fastq/*
for f in $FILES
do
  f_base="$(basename -- $f .fastq)"
  echo "Processing $f_base file..."

  # run strobmers
  python strobe_match.py --queries $f --references $f --outfolder $outbase/results/strobmers_r_vs_r/ --prefix $f_base --k 15 --w 10
  python $experiment_dir/print_hit_statistics.py $outbase/results/strobmers_r_vs_r/$f_base.tsv --refs $f --outfolder $outbase/ --method strobemers  --setting r_vs_r

  #run kmers
  python strobe_match.py --queries $f --references $f --outfolder $outbase/results/kmers_r_vs_r/ --prefix $f_base --k 30 --kmer_index --w 10
  python $experiment_dir/print_hit_statistics.py $outbase/results/kmers_r_vs_r/$f_base.tsv --refs $f --outfolder $outbase/ --method kmers --setting r_vs_r
done


# Plot1: Sort plot sirvs by length on x-axis! normalized_hit_length.csv
# Plot2: Plot total coverage of hits  covarage.csv
# Plot3 number of hits per read nr_hits.csv

python $experiment_dir/plots.py $normalized_hit_length_file $coverage_file $nr_hits_file $plot_file".pdf"


