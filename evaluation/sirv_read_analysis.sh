#!/bin/bash

# RUN scripts e.g. as:   ./sirv_match_lengths.sh /Users/kxs624/Documents/workspace/strobemers/  /Users/kxs624/tmp/STROBEMERS/sirv_match_lengths/  /Users/kxs624/Documents/data/ont/lc19_pcs109_subsample_full_length_pychopper2_phmmer.fq  100

inbase=$1
outbase=$2
original_reads=$3
nr_reads=$4

experiment_dir=$inbase"/evaluation/"

mkdir -p $outbase

# alias pypy_run="/Users/kxs624/Downloads/pypy3.7-v7.3.3-osx64/bin/./pypy3"

IFS=$'\n'       # make newlines the only separator
# set -f          # disable globbing


########### READS TO REFERENCES #######################
#######################################################

normalized_hit_length_file=$outbase/"normalized_hit_length.csv"
nr_hits_file=$outbase/"nr_hits.csv"
coverage_file=$outbase/"coverage.csv"


# results_file has format: method,read_acc,ref_id,ref_len,match_length,normalized_match_length,
echo -ne  "method\tread_id\tref_id\tref_length\tmatch_length\tnormalized_match_length"$'\n' > $normalized_hit_length_file
echo -ne  "method\tread_id\tref_id\tref_length\tnr_hits"$'\n' > $nr_hits_file
echo -ne  "method\tread_id\tref_id\tref_length\tcoverage"$'\n' > $coverage_file

# align original reads with minimap2 
# original_reads_mapped=$outbase/original_reads.sam
# minimap2 -a --eqx -k 10 -w 1 $inbase/test_data/sirv_transcripts.fasta  $original_reads  > $original_reads_mapped
original_reads_mapped=/Users/kxs624/Documents/data/ont/sirv/cDNA/lc19_pcs109_subsample_full_length_pychopper2_phmmer.sam  #premapped 

# Subsample nr_reads reads from each transcript
echo python $experiment_dir/sample_reads.py $original_reads $inbase/data/sirv_transcripts.fasta $original_reads_mapped $outbase/fastq $nr_reads
# python $experiment_dir/sample_reads.py $original_reads $inbase/data/sirv_transcripts.fasta $original_reads_mapped $outbase/fastq $nr_reads
###############

mkdir -p $outbase/results/randstrobes2/
mkdir -p $outbase/results/randstrobes3/
mkdir -p $outbase/refs/

FILES=$outbase/fastq/*
for f in $FILES
do
  f_base="$(basename -- $f .fastq)"
  echo "Processing $f_base file..."

  # run randstrobes
  /Users/kxs624/Downloads/pypy3.7-v7.3.3-osx64/bin/./pypy3 strobe_match.py --queries $f --references $outbase/refs/$f_base"_ref.fastq" --outfolder $outbase/results/randstrobes2/ --prefix $f_base --k 15 --n 2 --w 10 --n 2
  python $experiment_dir/print_hit_statistics.py $outbase/results/randstrobes2/$f_base.tsv --refs $outbase/refs/$f_base"_ref.fastq" --outfolder $outbase/ --method 'randstrobes-(2,15,20,70)'

  /Users/kxs624/Downloads/pypy3.7-v7.3.3-osx64/bin/./pypy3 strobe_match.py --queries $f --references $outbase/refs/$f_base"_ref.fastq" --outfolder $outbase/results/randstrobes3/ --prefix $f_base --k 10 --n 3 --w 10 --n 3
  python $experiment_dir/print_hit_statistics.py $outbase/results/randstrobes3/$f_base.tsv --refs $outbase/refs/$f_base"_ref.fastq" --outfolder $outbase/ --method 'randstrobes-(3,10,20,70)'

  # run minstrobes
  /Users/kxs624/Downloads/pypy3.7-v7.3.3-osx64/bin/./pypy3 strobe_match.py --queries $f --references $outbase/refs/$f_base"_ref.fastq" --outfolder $outbase/results/minstrobes2/ --prefix $f_base --k 15 --n 2 --w 10 --n 2 --minstrobe_index
  python $experiment_dir/print_hit_statistics.py $outbase/results/minstrobes2/$f_base.tsv --refs $outbase/refs/$f_base"_ref.fastq" --outfolder $outbase/ --method 'minstrobes-(2,15,20,70)'

  /Users/kxs624/Downloads/pypy3.7-v7.3.3-osx64/bin/./pypy3 strobe_match.py --queries $f --references $outbase/refs/$f_base"_ref.fastq" --outfolder $outbase/results/minstrobes3/ --prefix $f_base --k 10 --n 3 --w 10 --n 3 --minstrobe_index
  python $experiment_dir/print_hit_statistics.py $outbase/results/minstrobes3/$f_base.tsv --refs $outbase/refs/$f_base"_ref.fastq" --outfolder $outbase/ --method 'minstrobes-(3,10,20,70)'

  #run kmers
  /Users/kxs624/Downloads/pypy3.7-v7.3.3-osx64/bin/./pypy3 strobe_match.py --queries $f --references $outbase/refs/$f_base"_ref.fastq" --outfolder $outbase/results/kmers/ --prefix $f_base --k 30 --kmer_index --w 10
  python $experiment_dir/print_hit_statistics.py $outbase/results/kmers/$f_base.tsv --refs $outbase/refs/$f_base"_ref.fastq" --outfolder $outbase/ --method kmers 
done


# # Plot1: Sort plot sirvs by length on x-axis! normalized_hit_length.csv
# # Plot2: Plot total coverage of hits  covarage.csv
# # Plot3 number of hits per read nr_hits.csv

python $experiment_dir/plots.py $coverage_file $nr_hits_file $normalized_hit_length_file  $outbase

#######################################################
#######################################################
#######################################################





# ########### READS TO READS ############################
# #######################################################

# normalized_hit_length_file=$outbase/"normalized_hit_length_read_vs_read.csv"
# nr_hits_file=$outbase/"nr_hits_read_vs_read.csv"
# coverage_file=$outbase/"coverage_read_vs_read.csv"


# # results_file has format: method,read_acc,ref_id,ref_len,match_length,normalized_match_length,
# echo -ne  "method\tread_id\tref_id\tref_length\tmatch_length\tnormalized_match_length"$'\n' > $normalized_hit_length_file
# echo -ne  "method\tread_id\tref_id\tref_length\tnr_hits"$'\n' > $nr_hits_file
# echo -ne  "method\tread_id\tref_id\tref_length\tcoverage"$'\n' > $coverage_file

# # align original reads with minimap2 
# # original_reads_mapped=$outbase/original_reads.sam
# # minimap2 -a --eqx -k 10 -w 1 $inbase/test_data/sirv_transcripts.fasta  $original_reads  > $original_reads_mapped
# original_reads_mapped=/Users/kxs624/Documents/data/ont/sirv/cDNA/lc19_pcs109_subsample_full_length_pychopper2_phmmer.sam  #premapped 

# # Subsample nr_reads reads from each transcript
# echo python $experiment_dir/sample_reads.py $original_reads $inbase/data/sirv_transcripts.fasta $original_reads_mapped $outbase/fastq $nr_reads
# # python $experiment_dir/sample_reads.py $original_reads $inbase/data/sirv_transcripts.fasta $original_reads_mapped $outbase/fastq $nr_reads
# ###############

# mkdir -p $outbase/results/randstrobes2/
# mkdir -p $outbase/results/randstrobes3/
# mkdir -p $outbase/refs/

# FILES=$outbase/fastq/*
# for f in $FILES
# do
#   f_base="$(basename -- $f .fastq)"
#   echo "Processing $f_base file..."

#   # run randstrobes
#   /Users/kxs624/Downloads/pypy3.7-v7.3.3-osx64/bin/./pypy3 strobe_match.py --queries $f --references $f --outfolder $outbase/results/randstrobes2_r_vs_r/ --prefix $f_base --k 15 --n 2 --w 10
#   python $experiment_dir/print_hit_statistics.py $outbase/results/randstrobes2_r_vs_r/$f_base.tsv --refs $f --outfolder $outbase/ --method 'randstrobes-(2,15,20,70)'  --setting r_vs_r

#   /Users/kxs624/Downloads/pypy3.7-v7.3.3-osx64/bin/./pypy3 strobe_match.py --queries $f --references $f --outfolder $outbase/results/randstrobes3_r_vs_r/ --prefix $f_base --k 10 --n 3 --w 10
#   python $experiment_dir/print_hit_statistics.py $outbase/results/randstrobes3_r_vs_r/$f_base.tsv --refs $f --outfolder $outbase/ --method 'randstrobes-(3,10,20,70)'  --setting r_vs_r

#   # run minstrobes
#   /Users/kxs624/Downloads/pypy3.7-v7.3.3-osx64/bin/./pypy3 strobe_match.py --queries $f --references $f --outfolder $outbase/results/minstrobes2_r_vs_r/ --prefix $f_base --k 15 --n 2 --w 10 --minstrobe_index
#   python $experiment_dir/print_hit_statistics.py $outbase/results/minstrobes2_r_vs_r/$f_base.tsv --refs $f --outfolder $outbase/ --method 'minstrobes-(2,15,20,70)'  --setting r_vs_r

#   /Users/kxs624/Downloads/pypy3.7-v7.3.3-osx64/bin/./pypy3 strobe_match.py --queries $f --references $f --outfolder $outbase/results/minstrobes3_r_vs_r/ --prefix $f_base --k 10 --n 3 --w 10 --minstrobe_index
#   python $experiment_dir/print_hit_statistics.py $outbase/results/minstrobes3_r_vs_r/$f_base.tsv --refs $f --outfolder $outbase/ --method 'minstrobes-(3,10,20,70)'  --setting r_vs_r


#   #run kmers
#   /Users/kxs624/Downloads/pypy3.7-v7.3.3-osx64/bin/./pypy3 strobe_match.py --queries $f --references $f --outfolder $outbase/results/kmers_r_vs_r/ --prefix $f_base --k 30 --kmer_index --w 10
#   python $experiment_dir/print_hit_statistics.py $outbase/results/kmers_r_vs_r/$f_base.tsv --refs $f --outfolder $outbase/ --method kmers --setting r_vs_r
# done


# # Plot1: Sort plot sirvs by length on x-axis! normalized_hit_length.csv
# # Plot2: Plot total coverage of hits  covarage.csv
# # Plot3 number of hits per read nr_hits.csv

# python $experiment_dir/plots.py $coverage_file $nr_hits_file $normalized_hit_length_file $outbase/reads_vs_reads

# #######################################################
# #######################################################
# #######################################################

