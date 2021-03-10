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
                       --outfolder strobemer_output/  --k 15  --strobe_w_size 50


# Generate kmer matches (k=30) 
# between ONT SIRV reads and the true reference sequences

python strobe_match.py --queries data/sirv_transcripts.fasta \
                       --references data/ONT_sirv_cDNA_seqs.fasta \
                       --outfolder kmer_output/  --k 30 --kmer_index

# Reads vs reads matching

python strobe_match.py --queries data/sirv_transcripts.fasta \
                       --references data/sirv_transcripts.fasta \
                       --outfolder strobemer_output/ --k 15  --strobe_w_size 50
```

This will produce a file `matches.tsv` in the output folder. You can se a custom outfile name with the parameter `--prefix`.

## Output

Output format is similar to MUMmer:

```
>query_accession
ref_id  ref_pos query_pos   match_length_on_reference
```

Small example output from aligning sirv reads to transcripts (data above) which also highlights the stobemers strength compared to kmers. While kmers can give a more nuanced differentiation (compare read hits to `SIRV606` and `SIRV616`) both the sequences are good candidates for downstream processing. In this small example, the strobemers produce fewer hits/less output needed for post clustering of matches, e.g., for downstream clustering/alignment/mapping. Notice that randstobe hit positions are currently not deterministic due to hash seed is set at each new pyhon instantiation. I will fix the hash seed in future implementations.


**Strobemers (2,15,50)**
```
>41:650|d00e6247-9de6-485c-9b44-806023c51f13
SIRV606 35      92      487
SIRV616 35      92      473
>56:954|a23755a1-d138-489e-8efb-f119e679daf4
SIRV509 3       3       515
SIRV509 520     529     214
SIRV509 762     767     121
>106:777|0f79c12f-efed-4548-8fcc-49657f97a126
SIRV404 53      131     535
```

**kmers (k=30)**
```
>41:650|d00e6247-9de6-485c-9b44-806023c51f13
SIRV606 33      90      46
SIRV606 92      150     125
SIRV606 219     275     81
SIRV606 349     408     70
SIRV606 420     479     47
SIRV606 481     540     42
SIRV616 33      90      46
SIRV616 92      150     125
SIRV616 219     275     81
SIRV616 349     408     60
SIRV616 409     482     44
SIRV616 467     540     42
>56:954|a23755a1-d138-489e-8efb-f119e679daf4
SIRV509 68      72      141
SIRV509 230     233     100
SIRV509 331     335     105
SIRV509 435     442     40
SIRV509 475     483     36
SIRV509 579     585     41
SIRV509 621     627     46
SIRV509 695     701     44
SIRV509 812     815     53
>106:777|0f79c12f-efed-4548-8fcc-49657f97a126
SIRV404 53      131     58
SIRV404 128     208     127
SIRV404 283     364     30
SIRV404 422     494     142
```

## What is a MAM?

I have not yet rigorously defined a MAM (maximal approximate match) but I'm working on it as the development progresses. The aim is to output regions that _approximately match eachother_ on the reference and query (just like MEMs are exact matches between query and reference). The MAM regions can then be used to detect as candidate regions for alignment, clustering, or any other downstream analysis. 

The 'approximate' part is that there is a strobemer match, and the maximal is that we will merge matches that overlap on both query and reference sequence (if the order of the matches is the same on both the query and reference sequence). The 'maximal' _seems_ to be well defined if there are no repetitive matches in either the query or the reference. However, is is not trivial to compute if there are nested repeats in _both_ the query and reference. Currently, this is what is implemented:

For kmers, any two k-mer matches spanning positions `(q_1, q_1+k)` and `(q_2, q_2+k`) on the query and positions `(r_1, r_1+k)` and `(r_2, r_2+k)` on the reference where `q_1 <= q_2 <= q_1+k <= q_2+k` and `r_1 <= r_2 <= r_1+k <= r_2+k` are merged into one match of length `r_2+k - r_1`. Any chain of such overlapping matches are merged into one match. 


For strobemers, `strobe_match.py` saves the positions for both the first and second strobe. Two strobemers with start positions `(q_1, q'_1)` and `(q_2, q'_2)` on the query and `(r_1, r'_1)` and `(r_2, r'_2)` on the reference with length `k` strobes _overlap_ if `q_1 <= q_2 <= q'_1 +k` and `r_1 <= r_2 <= r'_2+k`. If there is an overlap the two strobes are merged into one match of length `max(q'_1+k, q'_2 + k) - q_1`. Notice that because of the random length between the strobes, we can either have `q_1 <= q_2 <= q'_1 <= q'_2` or `q_1 <= q_2 <= q'_2 <= q'_1`, hence we need the `max` function. Any chain of such overlapping matches are merged into one match. 


The tool currently have a known bug of not being able to merge matches when there exist a repeat occuring at least twice within _both_ the query and reference sequence. In this case the matches may become fragmented, i.e., not merged into MAMs.

## Proof of concept

I aligned ONT cDNA reads (meadian error rate 7.0%) from [this synthetic RNA dataset](https://www.ebi.ac.uk/ena/browser/view/PRJEB34849) to SIRV transcripts [available here](https://github.com/ksahlin/strobemers/blob/main/data/sirv_transcripts.fasta) using minimap2 with parameters `-k 10 -w 1` providing very sensitive/accurate alignment. The ONT reads have been processed into full length reads using pychopper. Thus, ideally all reads should span and align to the full transcript. I selected 100 reads aligning to each SIRV transcript (primary alignment), and compared match coverage and number of hits between the 100 reads and their reference transcript using kmers (k=30) and strobemers (n=2,k=15,w50) giving the same subsequence length of 30nt each.

The aim is that the matching should provide candidate regions/sequences to perform exact alignment against. Match coverage and number of hits are two important features for sequence matching. Match coverage in this experiment is the fraction of reference sequence covered. Number of hits is the number of MAMs per read. We want the match coverage to be high in this experiment since we know that the reads align well to the their respective SIRV reference. However, together with a high coverage, we want the number of matches to be as low as possible (where 1 is best) in order for fast post prosessing/clustering of matches (a.k.a. seeds) and low disk space. 

Below I show the match coverage and number of hits for strobemers and kmers in this experiment separated for each of the 63 SIRVs with more than 100 primary alignment. The line shows the mean and the shaded area around the line is the standard deviation of the data (i.e., coverage/nr matches) for each SIRV.


![match coverage](data/plot_coverage.png)
![number of hits](data/plot_nr_hits.png)

The two above metrics could be studied from another angle, which is the match length normalized with the SIRV transcript length. The plot below shows the mean normalized match length for kmers and strobemers.

![norm match length](data/plot_normalized_match_length.png)

CREDITS
----------------

Kristoffer Sahlin, Strobemers: an alternative to k-mers for sequence comparison, bioRxiv 2021.01.28.428549; doi: https://doi.org/10.1101/2021.01.28.428549

Paper found [here](https://www.biorxiv.org/content/10.1101/2021.01.28.428549v1)