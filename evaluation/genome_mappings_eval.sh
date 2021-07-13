#!/bin/bash

# RUN scripts e.g. as:   ./runtime_cpp.sh /Users/kxs624/Documents/data/genomes/human/HG_38/GRCh38_chr21.fa

genome1=$1
genome2=$2
# alias pypy_run="/Users/kxs624/Downloads/pypy3.7-v7.3.3-osx64/bin/./pypy3"

IFS=$'\n'       # make newlines the only separator
# set -f          # disable globbing


#######################################################

# for k in 100 500 #30
# do
# 	# mummer MEM
# 	/usr/bin/time -l mummer -F -maxmatch -l $k -b  $genome1 $genome2 > tmp.tsv 2> runtime.txt
# 	echo -n "MUMmer & MEM & " $k " & " 
# 	python genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2
# 	# # mummer MUM
# 	/usr/bin/time -l mummer -F -l $k -mum -b $genome1 $genome2 > tmp.tsv 2> runtime.txt
# 	echo -n "MUMmer & MUM & " $k " & " 
# 	python genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2
# done

# # Strobemap all
# /usr/bin/time -l StrobeMap -k 30 -v 31 -c kmers -o tmp.tsv $genome1 $genome2 2> runtime.txt 1> stdout.txt
# echo -n "StrobeMap & all & 30 & " 
# python genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2

# # Strobemap unique
# /usr/bin/time -l StrobeMap -k 30 -v 31 -c kmers -u -o tmp.tsv $genome1 $genome2 2> runtime.txt 1> stdout.txt
# echo -n "StrobeMap & unique & 30 & " 
# python genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2

# /usr/bin/time -l StrobeMap -n 3 -k 10 -v 11 -w 100 -c randstrobes -o tmp.tsv $genome1 $genome2 2> runtime.txt 1> stdout.txt
# echo -n "StrobeMap & all & (3, 10, 11, 100) & " 
# python genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2

# /usr/bin/time -l StrobeMap -u -n 3 -k 10 -v 11 -w 100 -c randstrobes -o tmp.tsv $genome1 $genome2 2> runtime.txt 1> stdout.txt
# echo -n "StrobeMap & unique & (3, 10, 11, 100) & " 
# python genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2

# /usr/bin/time -l StrobeMap -n 2 -k 30 -v 31 -w 100 -c randstrobes -o tmp.tsv $genome1 $genome2 2> runtime.txt  1> stdout.txt
# echo -n "StrobeMap & all & (2, 30, 31, 100) & " 
# python genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2

# /usr/bin/time -l StrobeMap -u -n 2 -k 30 -v 31 -w 100 -c randstrobes -o tmp.tsv $genome1 $genome2 2> runtime.txt 1> stdout.txt
# echo -n "StrobeMap & unique & (2, 30, 31, 100) & " 
# python genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2

# /usr/bin/time -l StrobeMap -n 3 -k 30 -v 31 -w 100 -c randstrobes -o tmp.tsv $genome1 $genome2 2> runtime.txt  1> stdout.txt
# echo -n "StrobeMap & all & (3, 30, 31, 100) & " 
# python genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2

# /usr/bin/time -l StrobeMap -u -n 3 -k 30 -v 31 -w 100 -c randstrobes -o tmp.tsv $genome1 $genome2 2> runtime.txt 1> stdout.txt
# echo -n "StrobeMap & unique & (3, 30, 31, 100) & " 
# python genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2

# /usr/bin/time -l StrobeMap -n 2 -k 30 -v 200 -w 400 -c randstrobes -o tmp.tsv $genome1 $genome2 2> runtime.txt 1> stdout.txt
# echo -n "StrobeMap & all & (2, 30, 200, 400) & " 
# python genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2

# /usr/bin/time -l StrobeMap -n 3 -k 30 -v 200 -w 400 -c randstrobes -o tmp.tsv $genome1 $genome2 2> runtime.txt 1> stdout.txt
# echo -n "StrobeMap & all & (3, 30, 200, 400) & " 
# python genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2

# /usr/bin/time -l StrobeMap -u -n 3 -k 30 -v 200 -w 400 -c randstrobes -o tmp.tsv $genome1 $genome2 2> runtime.txt 1> stdout.txt
# echo -n "StrobeMap & unique & (3, 30, 200, 400) & " 
# python genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2


/usr/bin/time -l StrobeMap -n 2 -k 30 -v 500 -w 600 -c randstrobes -o tmp.tsv $genome1 $genome2 2> runtime.txt 1> stdout.txt
echo -n "StrobeMap & all & (2, 30, 500, 600) & " 
python genome_mapping_metrics.py tmp.tsv runtime.txt --refs $genome2


#######################################################
#######################################################
#######################################################


