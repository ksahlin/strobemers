Strobemers
===========

A repository for generating strobemers and evalaution.  

# Using the library

First off, this is a prototype implementation created for the analysis in the paper describing strobemers. As Python is terribly inefficient (and memory requiring) with string manipulation, any real implementation and usage of strobemers should probably happen in a low level language. However, the Python code is sufficient as proof of concept. 


The `indexing.py` module located in the `modules` folder contains functions for generating k-mers, spaced k-mers, minimizers, and strobemers (both minstrobe and randstrobe) of order 2 and 3. For randstrobes, there are two ways to create them. The first way is with the function `randstrobes`, which takes a string, k-mer size, and window size and returns a dictionary with positions as keys and the randstrobe sequences (strings) as values. The second way is to call `randstrobes_iter` which is a generator. Similarly to `randstrobes`, `randstrobes_iter` takes a string, k-mer size, and window size, but instead yields randstrobes from the sequence. `randstrobes_iter` can be used as follows

```
from modules import indexing
all_mers = defaultdict(int)
for s in indexing.randstrobes_iter(seq, k_size, order = 2, N_1 = 50 ):
    all_mers[s] += 1
```

minstobes has the same interface as randstrobes with functions `minstrobes` and `minstrobes_iter`

The script `analysis1and1.py` at the top level in this repository contains code for generating experiments described in section "strobemer vs k-mer matching" in the paper. The script `analysis3.py` at the top level in this repository contains code for generating experiments described in section "strobemer vs k-mer uniqueness" in the paper. The  script `analysis4.py` is simply a copy of analysis3, where I had planned to implement som proof of concept analysis of biological data if I can find some reasonable experiment using this prototype library. 



CREDITS
----------------

Kristoffer Sahlin, Strobemers: an alternative to k-mers for sequence comparison, bioRxiv 2021.01.28.428549; doi: https://doi.org/10.1101/2021.01.28.428549

Paper found [here](https://www.biorxiv.org/content/10.1101/2021.01.28.428549v1)