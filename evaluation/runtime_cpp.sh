#!/bin/bash

# RUN scripts e.g. as:   ./runtime_cpp.sh /Users/kxs624/Documents/data/genomes/human/HG_38/GRCh38_chr21.fa

genome=$1

# alias pypy_run="/Users/kxs624/Downloads/pypy3.7-v7.3.3-osx64/bin/./pypy3"

IFS=$'\n'       # make newlines the only separator
# set -f          # disable globbing


########### READS TO REFERENCES #######################
#######################################################

let "k= 30"
let "v= k+1"
let "w= $v+$size"
echo kmers,$n,$k,$size,$v,$w,$genome
StrobeMap_indextime -k $k -v $v -w $w -c kmers -o /Users/kxs624/Downloads/strobemap_test.tsv $genome $genome
echo

for n in 2 3 
do

	for size in 20 40 60 80 100
	do
	  let "k= 30/$n"
	  let "v= k+1"
	  let "w= $v+$size"
	  echo minstrobes,$n,$k,$size,$v,$w
	  StrobeMap_indextime -k $k -n $n -v $v -w $w -c minstrobes  -o /Users/kxs624/Downloads/strobemap_test.tsv $genome $genome
	done
	echo

	for size in 20 40 60 80 100
	do
	  let "k= 30/$n"
	  let "v= k+1"
	  let "w= $v+$size"
	  echo hybridstrobes,$n,$k,$size,$v,$w
	  StrobeMap_indextime -k $k -n $n -v $v -w $w -c hybridstrobes  -o /Users/kxs624/Downloads/strobemap_test.tsv $genome $genome
	done
	echo


	for size in 20 40 60 80 100
	do
	  let "k= 30/$n"
	  let "v= k+1"
	  let "w= $v+$size"
	  echo randstrobes,$n,$k,$size,$v,$w
	  StrobeMap_indextime -k $k -n $n -v $v -w $w -c randstrobes  -o /Users/kxs624/Downloads/strobemap_test.tsv $genome $genome
	done
	echo
done


#######################################################
#######################################################
#######################################################


