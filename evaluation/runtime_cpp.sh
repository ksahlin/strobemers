#!/bin/bash


#######################################################
#######################################################
########## INFORMATION ABOUT THE SCRIPT ###############
#######################################################
#######################################################

# RUN scripts e.g. as:   ./runtime_cpp.sh /Users/kxs624/Documents/data/genomes/human/HG_38/GRCh38_chr21.fa

# StrobeMap_indextime called belo is a modified version where
# line 681 and 683 in main.cpp has been uncommented and then compiled. 
# This assures that only the construction time of strobemers/k-mers is 
# counted and the only output produced is the total time in seconds.



genome1=$1
genome2=$2
genome3=$3


IFS=$'\n'       # make newlines the only separator
# set -f          # disable globbing




let "k= 30"
let "v= k+1"
let "w= $v+10"
echo -n kmers "&" $k "& "  #,0,$v,$w,$genome1
# /usr/bin/time -l # for peak mem
StrobeMap_indextime -k $k -v $v -w $w -c kmers -o /Users/kxs624/Downloads/strobemap_test.tsv $genome1 $genome1
echo -n "& "
StrobeMap_indextime -k $k -v $v -w $w -c kmers -o /Users/kxs624/Downloads/strobemap_test.tsv $genome2 $genome2
echo -n "& "
StrobeMap_indextime -k $k -v $v -w $w -c kmers -o /Users/kxs624/Downloads/strobemap_test.tsv $genome3 $genome3
echo

for n in 2 3 
do

	for size in 20 40 60 80 100
	do
	  let "k= 30/$n"
	  let "v= k+1"
	  let "w= $v+$size"
	  echo -n minstrobes "&" "("$n, $k, $v, $w")" "& " #,$v,$w
	  # /usr/bin/time -l # for peak mem
	  StrobeMap_indextime -k $k -n $n -v $v -w $w -c minstrobes  -o /Users/kxs624/Downloads/strobemap_test.tsv $genome1 $genome1
	  echo -n "& "
	  StrobeMap_indextime -k $k -n $n -v $v -w $w -c minstrobes  -o /Users/kxs624/Downloads/strobemap_test.tsv $genome2 $genome2
	  echo -n "& "
	  StrobeMap_indextime -k $k -n $n -v $v -w $w -c minstrobes  -o /Users/kxs624/Downloads/strobemap_test.tsv $genome3 $genome3
	  echo
	done


	for size in 20 40 60 80 100
	do
	  let "k= 30/$n"
	  let "v= k+1"
	  let "w= $v+$size"
	  echo -n hybridstrobes "&" "("$n, $k, $v, $w")" "& " #,$v,$w
	  # /usr/bin/time -l # for peak mem
	  StrobeMap_indextime -k $k -n $n -v $v -w $w -c randstrobes  -o /Users/kxs624/Downloads/strobemap_test.tsv $genome1 $genome1
	  echo -n "& "
	  StrobeMap_indextime -k $k -n $n -v $v -w $w -c randstrobes  -o /Users/kxs624/Downloads/strobemap_test.tsv $genome2 $genome2
	  echo -n "& "
	  StrobeMap_indextime -k $k -n $n -v $v -w $w -c randstrobes  -o /Users/kxs624/Downloads/strobemap_test.tsv $genome3 $genome3
	  echo
	done
	


	for size in 20 40 60 80 100
	do
	  let "k= 30/$n"
	  let "v= k+1"
	  let "w= $v+$size"
	  echo -n randstrobes "&" "("$n, $k, $v, $w")" "& " #,$v,$w
	  # /usr/bin/time -l # for peak mem
	  StrobeMap_indextime -k $k -n $n -v $v -w $w -c hybridstrobes  -o /Users/kxs624/Downloads/strobemap_test.tsv $genome1 $genome1
	  echo -n "& "
	  StrobeMap_indextime -k $k -n $n -v $v -w $w -c hybridstrobes  -o /Users/kxs624/Downloads/strobemap_test.tsv $genome2 $genome2
	  echo -n "& "
	  StrobeMap_indextime -k $k -n $n -v $v -w $w -c hybridstrobes  -o /Users/kxs624/Downloads/strobemap_test.tsv $genome3 $genome3
	  echo
	done
	
done


#######################################################
#######################################################
#######################################################


