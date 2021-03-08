Strobemers
===========

A repository for generating strobemers and evaluation.  


First off, this is a prototype implementation created for the analysis in the paper describing strobemers. As Python is terribly inefficient (and memory requiring) with string manipulation, any real implementation and usage of strobemers should probably happen in a low level language. However, the Python code is sufficient as proof of concept. 

The repository consists of a library and a tool. The library `indexing.py` contains functions and generators for creating the datastructures used in the evaluation of the paper. The tool `strobe_match.py` is a program which roughly has the same interface as `nucmer`. `strobe_match.py` takes a reference and queries file in fasta or fastq format. It produces MAMs (Maximal Approximate non-overlapping Matches) between the queries and references and outputs then in a format simular to nucmer/MUMmer. Details Below.


# Using the library

The `indexing.py` module located in the `modules` folder contains functions for generating k-mers, spaced k-mers, minimizers, and strobemers (both minstrobe and randstrobe) of order 2 and 3. For randstrobes, there are two ways to create them. The first way is with the function `randstrobes`, which takes a string, k-mer size, and window size and returns a dictionary with positions as keys and the randstrobe sequences (strings) as values. The second way is to call `randstrobes_iter` which is a generator. Similarly to `randstrobes`, `randstrobes_iter` takes a string, k-mer size, and window size, but instead yields randstrobes from the sequence. `randstrobes_iter` can be used as follows

```
from modules import indexing
all_mers = defaultdict(int)
for s in indexing.randstrobes_iter(seq, k_size, order = 2, w_1 = 50 ):
    all_mers[s] += 1
```

or with order 3 strobemers:

```
from modules import indexing
all_mers = defaultdict(int)
for s in indexing.randstrobes_iter(seq, k_size, order = 3, w_1 = 25, w_2 = 25 ):
    all_mers[s] += 1
```

minstobes has the same interface as randstrobes with functions `minstrobes` and `minstrobes_iter`

The script `analysis1and1.py` at the top level in this repository contains code for generating experiments described in section "strobemer vs k-mer matching" in the paper. The script `analysis3.py` at the top level in this repository contains code for generating experiments described in section "strobemer vs k-mer uniqueness" in the paper. The  script `analysis4.py` is simply a copy of analysis3, where I had planned to implement som proof of concept analysis of biological data if I can find some reasonable experiment using this prototype library. 


# Using the tool

Currently, the tool only implements order 2 randstrobes and kmers. The tool produces MAMs (Maximal Approximate non-overlapping Matches; see explanation below) for both strobemers and kmers. Test data is found in the folder `data` in this repository.
Here are some example uses:

```
# Generate randstrobe matches (randstrobe parametrization (2,15,50)) 
# between ONT SIRV reads and the true reference sequences

python strobe_match.py --queries data/sirv_transcripts.fasta \
                       --references data/ONT_sirv_cDNA_seqs.fasta \
                       --outfile strobemer_hits.txt --k 15  --strobe_w_size 50


# Generate kmer matches (k=30) 
# between ONT SIRV reads and the true reference sequences

python strobe_match.py --queries data/sirv_transcripts.fasta \
                       --references data/ONT_sirv_cDNA_seqs.fasta \
                       --outfile strobemer_hits.txt --k 30 --kmer_index

# Reads vs reads matching

python strobe_match.py --queries data/sirv_transcripts.fasta \
                       --references data/sirv_transcripts.fasta \
                       --outfile strobemer_hits.txt --k 15  --strobe_w_size 50
```

## Output

Output format is similar to MUMmer:

```
>query_accession
ref_id  ref_pos query_pos   match_length_on_query
```

Small example output from aligning sirv reads to transcripts (data above) which also highlights the stobemers strength compared to kmers. While kmers can give a more nuanced differentiation (compare read hits to `SIRV606` and `SIRV616`) both the sequences are good candidates for downstream processing. In this small example, the strobemers produce fewer hits/less output needed for post clustering of matches, e.g., for downstream clustering/alignment/mapping.


**Strobemers (2,15,50)**
```
>41:650|d00e6247-9de6-485c-9b44-806023c51f13
SIRV606 33      90      485
SIRV616 33      90      485
>55:474|ed91d713-ef66-4526-9f8f-3ca0206564ab
SIRV704 2       0       409
>56:842|9d29a5f2-3704-4d73-bcf2-8a3f3d8b66ce
SIRV309 19      19      220
SIRV309 297     297     136
SIRV309 454     452     331
```

**kmers (k=30)**
```
>41:650|d00e6247-9de6-485c-9b44-806023c51f13
SIRV606 33      90      46
SIRV606 92      150     124
SIRV606 219     275     81
SIRV606 349     408     70
SIRV606 420     479     47
SIRV606 481     540     42
SIRV616 33      90      46
SIRV616 92      150     124
SIRV616 219     275     81
SIRV616 349     408     60
SIRV616 409     482     44
SIRV616 467     540     42
>55:474|ed91d713-ef66-4526-9f8f-3ca0206564ab
SIRV704 2       0       92
SIRV704 95      91      80
SIRV704 224     218     58
SIRV704 331     321     88
>56:842|9d29a5f2-3704-4d73-bcf2-8a3f3d8b66ce
SIRV309 43      45      107
SIRV309 187     189     50
SIRV309 297     297     136
SIRV309 470     465     59
SIRV309 530     523     33
SIRV309 586     579     38
SIRV309 660     650     133
```

## What is a MAM?

Currently it not yet properly defined but I'm working on it as the development progresses. The aim is to output regions that approxiamtely match eachother that can be used as candidate regions for alignment, clustering, etc. What is currently output can be thought of as regions with overlapping matches (overlap both in query and reference required). 

For kmers, this means that any two k-mer matches spanning positions `(q_1, q_1+k)` and `(q_2, q_2+k`) on the query and positions `(r_1, r_1+k)` and `(r_2, r_2+k)` where `q_1 <= q_2 <= q_1+k <= q_2+k` and `r_1 <= r_2 <= r_1+k <= r_2+k` are merget into one match of length `q_2+k - q_1`. Any chain of such overlapping matches are merged into one match. 


For strobemers, `strobe_match.py` saves the positions for both the first and second strobe. Two strobemers with start positions `(q_1, p_1)` and `(q_2, p_2)` on the query and length `k` strobes overlap if `q_1 <= q_2 <= p_1` and the same holds on the reference. Notice that because of the random length between the strobes, we can either have `q_1 <= q_2 <= p_1 <= p_2` or `q_1 <= q_2 <= p_2 <= p_1`. They are merged into one match of length `max(p_2 + k, p_1+k) - q_1`. Any chain of such overlapping matches are merged into one match. 


The tool currently have a known bug of not being able to merge matches when there exist a repeat occuring at least twice within _both_ the query and reference sequence. In this case the matches may become fragmented, i.e., not merged into MAMs.


CREDITS
----------------

Kristoffer Sahlin, Strobemers: an alternative to k-mers for sequence comparison, bioRxiv 2021.01.28.428549; doi: https://doi.org/10.1101/2021.01.28.428549

Paper found [here](https://www.biorxiv.org/content/10.1101/2021.01.28.428549v1)