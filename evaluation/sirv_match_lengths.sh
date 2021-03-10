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



normalized_hit_length_file=$outbase/"normalized_hit_length.csv"
nr_hits_file=$outbase/"nr_hits.csv"
coverage_file=$outbase/"coverage.csv"

plot_file=$outbase/"summary"

# results_file has format: method,read_acc,ref_id,ref_len,match_length,normalized_match_length,
echo -n  "method","read_id","ref_id","ref_length","match_length","normalized_match_length"$'\n' > $normalized_hit_length_file
echo -n  "method","read_id","ref_id","ref_length","nr_hits"$'\n' > $nr_hits_file
echo -n  "method","read_id","ref_id","ref_length","coverage"$'\n' > $coverage_file

# align original reads with minimap2 
# original_reads_mapped=$outbase/original_reads.sam
# minimap2 -a --eqx -k 10 -w 1 $inbase/test_data/sirv_transcripts.fasta  $original_reads  > $original_reads_mapped
original_reads_mapped=/Users/kxs624/Documents/data/ont/sirv/cDNA/lc19_pcs109_subsample_full_length_pychopper2_phmmer.sam  #premapped 

# Subsample nr_reads reads from each transcript
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
  python strobe_match.py --queries $f --references $outbase/refs/$f_base"_ref.fastq" --outfolder $outbase/results/strobmers/ --prefix $f_base --k 15
  python $experiment_dir/print_hit_statistics.py $outbase/results/strobmers/$f_base.tsv --refs $outbase/refs/$f_base"_ref.fastq" --outfolder $outbase/ --method strobemers  

  #run kmers
  python strobe_match.py --queries $f --references $outbase/refs/$f_base"_ref.fastq" --outfolder $outbase/results/kmers/ --prefix $f_base --k 30 --kmer_index
  python $experiment_dir/print_hit_statistics.py $outbase/results/kmers/$f_base.tsv --refs $outbase/refs/$f_base"_ref.fastq" --outfolder $outbase/ --method kmers 
done


# Plot1: Sort plot sirvs by length on x-axis! normalized_hit_length.csv
# Plot2: Plot total coverage of hits  covarage.csv
# Plot3 number of hits per read nr_hits.csv

# python $experiment_dir/plots.py $normalized_hit_length_file $coverage_file $nr_hits_file $plot_file".pdf"



# ##### EXPERMENT FOR k AND w ################
# ############################################
# ############################################

# # for id in $(seq 1 1 10)    
# # do 
# #     mkdir -p $outbase/$id/fastq
# #     python $experiment_dir/sample_reads.py $original_reads $inbase/test_data/sirv_transcriptome.fasta $original_reads_mapped $outbase/$id/fastq > /dev/null
# #     sampled_transcripts=$outbase/$id/fastq/sampled_transcripts.fasta
# #     for depth in 3 5 10 20 
# #     do
# #         echo
# #         echo $id,$depth
# #         reads=$outbase/$id/fastq/$depth
# #         fastq2fasta $reads.fastq $outbase/$id/fastq/$depth.fasta #&> /dev/null


# #         original_eval_out=$outbase/$id/$depth/original/evaluation
# #         mkdir -p $original_eval_out
# #         minimap2 -a --eqx -k 14 -w 4 $sampled_transcripts $reads.fasta  > $reads.sam 2>/dev/null
# #         python $experiment_dir/get_error_rates.py $sampled_transcripts  $reads.sam $original_eval_out/results.csv 

# #         for k_param in 7 8 9
# #         do
# #             for window in 0 2 4 
# #             do
# #                 echo $k_param,$window

# #                 w_param=$(( $k_param + $window ))
# #                 corrected_approximate=$outbase/$id/isoncorrect/$depth/$k_param/$w_param/approximate/
# #                 eval_out=$outbase/$id/isoncorrect/$depth/$k_param/$w_param/approximate/evaluation
# #                 mkdir -p $eval_out
# #                 python $inbase/isONcorrect3 --fastq $reads.fastq  --outfolder $corrected_approximate --k $k_param --w $w_param  &> /dev/null            
# #                 minimap2 -a --eqx -k 14 -w 4 $sampled_transcripts $corrected_approximate/corrected_reads.fastq  > $corrected_approximate/corrected_reads.sam 2>/dev/null
# #                 python $experiment_dir/get_error_rates.py $sampled_transcripts $corrected_approximate/corrected_reads.sam $eval_out/results.csv 
                
# #                 awk -F "," -v awk_id=$id -v awk_depth=$depth -v awk_k=$k_param -v awk_w=$window  '{if (NR!=1) {print awk_id",approx,"awk_depth","awk_k","awk_w","$0}}'  $eval_out/results.csv >> $results_file                
# #                 corrected_exact=$outbase/$id/isoncorrect/$depth/$k_param/$w_param/exact/
# #                 eval_out=$outbase/$id/isoncorrect/$depth/$k_param/$w_param/exact/evaluation
# #                 mkdir -p $eval_out
# #                 python $inbase/isONcorrect3 --fastq $reads.fastq  --outfolder $corrected_exact --k $k_param --w $w_param --exact  &> /dev/null            
# #                 minimap2 -a --eqx -k 14 -w 4 $sampled_transcripts $corrected_exact/corrected_reads.fastq  > $corrected_exact/corrected_reads.sam 2>/dev/null
# #                 python $experiment_dir/get_error_rates.py $sampled_transcripts $corrected_exact/corrected_reads.sam  $eval_out/results.csv 
                
# #                 awk -F "," -v awk_id=$id -v awk_depth=$depth -v awk_k=$k_param -v awk_w=$window '{if (NR!=1) {print awk_id",exact,"awk_depth","awk_k","awk_w","$0}}'  $eval_out/results.csv >> $results_file
# #                 awk -F "," -v awk_id=$id -v awk_depth=$depth -v awk_k=$k_param -v awk_w=$window '{if (NR!=1) {print awk_id",original,"awk_depth","awk_k","awk_w","$0}}'  $original_eval_out/results.csv >> $results_file

# #             done
# #         done
# #     done
# # done

# ############################################
# ############################################


# ####### EXPERMENT FOR xmax ################
# ############################################
# ############################################

# echo -n  "id","type","Depth","xmax","q_acc","r_acc","total_errors","error_rate","subs","ins","del"$'\n' > $results_file

# for id in $(seq 1 1 10)    
# do 
#     mkdir -p $outbase/$id/fastq
#     python $experiment_dir/sample_reads.py $original_reads $inbase/test_data/sirv_transcriptome.fasta $original_reads_mapped $outbase/$id/fastq > /dev/null
#     sampled_transcripts=$outbase/$id/fastq/sampled_transcripts.fasta
#     for depth in 3 5 10 20 
#     do
#         echo
#         echo $id,$depth
#         reads=$outbase/$id/fastq/$depth
#         fastq2fasta $reads.fastq $outbase/$id/fastq/$depth.fasta #&> /dev/null


#         original_eval_out=$outbase/$id/$depth/original/evaluation
#         mkdir -p $original_eval_out
#         minimap2 -a --eqx -k 14 -w 4 $sampled_transcripts $reads.fasta  > $reads.sam 2>/dev/null
#         python $experiment_dir/get_error_rates.py $sampled_transcripts  $reads.sam $original_eval_out/results.csv 

#         for xmax in 40 60 80 100 
#         do
#             corrected_approximate=$outbase/$id/isoncorrect/$depth/$xmax/approximate/
#             eval_out=$outbase/$id/isoncorrect/$depth/$xmax/approximate/evaluation
#             mkdir -p $eval_out
#             python $inbase/isONcorrect3 --fastq $reads.fastq  --outfolder $corrected_approximate --k 9 --w 9 --xmax $xmax  &> /dev/null            
#             minimap2 -a --eqx -k 14 -w 4 $sampled_transcripts $corrected_approximate/corrected_reads.fastq  > $corrected_approximate/corrected_reads.sam 2>/dev/null
#             python $experiment_dir/get_error_rates.py $sampled_transcripts $corrected_approximate/corrected_reads.sam $eval_out/results.csv 
            
#             awk -F "," -v awk_id=$id -v awk_depth=$depth -v awk_x=$xmax '{if (NR!=1) {print awk_id",approx,"awk_depth","awk_x","$0}}'  $eval_out/results.csv >> $results_file                
#             corrected_exact=$outbase/$id/isoncorrect/$depth/$k_param/$w_param/exact/
#             eval_out=$outbase/$id/isoncorrect/$depth/$k_param/$w_param/exact/evaluation
#             mkdir -p $eval_out
#             python $inbase/isONcorrect3 --fastq $reads.fastq  --outfolder $corrected_exact --k 9 --w 9 --xmax $xmax  --exact  &> /dev/null            
#             minimap2 -a --eqx -k 14 -w 4 $sampled_transcripts $corrected_exact/corrected_reads.fastq  > $corrected_exact/corrected_reads.sam 2>/dev/null
#             python $experiment_dir/get_error_rates.py $sampled_transcripts $corrected_exact/corrected_reads.sam  $eval_out/results.csv 
            
#             awk -F "," -v awk_id=$id -v awk_depth=$depth -v awk_x=$xmax   '{if (NR!=1) {print awk_id",exact,"awk_depth","awk_x","$0}}'  $eval_out/results.csv >> $results_file
#             awk -F "," -v awk_id=$id -v awk_depth=$depth -v awk_x=$xmax  '{if (NR!=1) {print awk_id",original,"awk_depth","awk_x","$0}}'  $original_eval_out/results.csv >> $results_file

#         done
#     done
# done

# python $experiment_dir/plot_error_rates.py $results_file $plot_file"_tot.pdf" error_rate

# python $experiment_dir/plot_error_rates.py $results_file $plot_file"_tot.pdf" error_rate
# python $experiment_dir/plot_error_rates.py $results_file $plot_file"_subs.pdf" subs
# python $experiment_dir/plot_error_rates.py $results_file $plot_file"_ind.pdf" ins
# python $experiment_dir/plot_error_rates.py $results_file $plot_file"_del.pdf" del

