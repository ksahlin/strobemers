Evaluation of methods for constructing randstrobes
===========

This evaluation will compare combinations of hashing and linking k-mers into randstrobes.


### Hashing

1. No hash - simple 2bit encoding of nucleotides
2. Thomas Wang hash
3. [xxhash](https://github.com/Cyan4973/xxHash)
4. [wyhash](https://github.com/wangyi-fudan/wyhash)
Let the hash function be denoted by h.

### Linking

1. Sahlin ( h(k_1) + h(k') ) % p 
2. Shen ( (h(k_1) + h(k') ) & p 
3. Sahlin bitcount( h(k_1) ^ h(k') ) 
4. Guo-Pibri h(k_1) ^ h(k') (global XOR), (also a variant: h(k_1 ^ k'))
5. Liu-Patro-Li h( k_1 || k' ) 
6. Liu-Patro-Li h( k_1 || k' ) (using wyhash for linking)

Hashing with method 1 and 2 is not compatible with linking method 5 and 6 as the total strobemers length can be larger than 32.


### Metrics

1. Construction time
2. "Randomness" selecting strobes:
    - Distribution over positions sampled for the second strobe k_2 (/Number of uniquely sampled positions)
    - Distance distribution between the `k_1` and `k_2` (goes between `w_min` to `w_max`)
    - W-spread of coordinates of p2, the difference in position of p2 compared to the p2's of previousluy W sampled strobes. Calculated as abs( min_{j-W < i < j}(p2(j) - p2(i))).
3. Uniqueness of final hash value:
    - How many unique matches between two sequences
4. Downstream metric
    - How it affects strobealign (accuracy, final runtime, mapping sites tried, calls to ksw2)

For randomness, neither of the first two measures suffices. Imagine having a an approach that resulted in a fixed sampling distance `d` between `k_1` and `k_2` (i.e., a spaced k-mer). This will give an optimal result on the first measures. 

Instead, imagine having a an approach that resulted in an approach where position `p+w_max` was sampled for `k_2` for all `w_max-w_min` consecutive `k_1`'s in `[p, p+w_max-w_min]`, then the approach would "jump" forward and sample `k_2` at position `p+2w_max-w_min` for the next `w_max-w_min` consecutive positions of `k_1` in `[p+w_max-w_min, q+2w_max-w_min]`. This would give an optimal uniform distribution over the second measure, but the positional correlation between consecutive strobes is extremely high.

The third measure is designed to measure local dispersity and would give a low score to both the scenarious above.

## Run benchmark

Compile

```
g++-11 -std=c++14 main.cpp index.cpp xxhash.c -lz -fopenmp -o randstrobe_benchmark -O3 -mavx2
```

Run

```
./randstrobe_benchmark -x [1,2,3] -l [link method 1-6] refs.fasta queries.fasta
```

If `-l 5` is specified, `-x` will use xxhash128.

```
$ ./randstrobe_benchmark

Randstrobe evaluation

StrobeMap [options] <references.fasta> <queries.fast[a/q]>
options:
  -k INT strobe length, limited to 32 [20]
  -v INT strobe w_min offset [21]
  -w INT strobe w_max offset [120]
  -t INT number of threads [3]
  -o name of output tsv-file [output.tsv]
  -x Choice of hash function to use; 1: nohash, 2: Thomas Wang hash, 3:xxhash [1]. 
  -l Choice of link function to use; 1: method1 (sahlin modulo), 2: method2 (shen bitwise AND), 3: method3 (guo_pibri XOR), 4: method4 (sahlin bitcount XOR), 5: method5 (Liu-Patro-Li, concatenation), 6: method6 (Liu-Patro-Li, concatenation using wyhash for linking).
```


## Current timings on an E coli genome

Below the `ecoli.fasta` is [this](https://www.ncbi.nlm.nih.gov/nuccore/NC_000913.3?report=fasta) E. coli genome.

The second argument can be any file, as it is no currently used. This file will be used when implementing matching statistics under the methods. Currently, only the positions genereated on the genome are logged and output in a file `positions.tsv` the run directory. `positions.tsv` contains a row for each genomic position. Therefore I do not recomment running this analysis on large genomes yet.

```
xxhash - Guo-Pibri
----------------
./randstrobe_benchmark -x 3 -l 3 -v 21 -w 120  ecoli.fasta dummy.fa

Total time hashing: 0.068509 s
Total time linking: 0.766273 s


xxhash - Liu-Patro-Li
---------------
./randstrobe_benchmark -l 5 -v 21 -w 120  ~/Downloads/sequence.fasta  dummy.fa
  
Total time hashing: 0.039633 s (hashing is just reading in k-mers, xxhash performed after concatenation)
Total time linking: 3.50736 s


xxhash - Sahlin2
----------------
 ./randstrobe_benchmark -x 3 -l 4 -v 21 -w 120  ~/Downloads/sequence.fasta  dummy.fa

Total time hashing: 0.065009 s
Total time linking: 0.798029 s
```

## Pseudo-randomness evaluation

For intuition about the problems we encounter, see [this figure](https://github.com/ksahlin/strobemers/blob/main/randstrobe_implementations/figures/clumpings_motivation.pdf).

To measure the pseudorandomness, the script [evaluate_sampling.py](https://github.com/ksahlin/strobemers/tree/main/randstrobe_implementations/evaluation) can be run on the `positions.tsv` file. It currently takes a long time (around 1min) for the ecoli genome. 

```
python evaluate_sampling.py positions.tsv outpath/to/randstrobe_evaluation
```

The output is seen below for two measures

```
Hash,Link,total_unique,most_repetitive,distance_nonuniformity,d_min2,d_min3,d_min4,d_min5,d_max2,d_max3,d_max4,d_max5,d_max6-10
xxh64,Guo-Pibri,2987438,11,198.52100000000007,1.19,1.17,1.16,1.14,1.12,1.5,3.42,16.16,377.76
xxh64,Liu-Patro-Li,2942552,9,153.82339999999996,1.01,1.01,1.01,1.01,1.01,1.01,0.95,0.91,1.18
```

### The metrics

Below is a brief description ofwhat the metrics mean. See [the figure](https://github.com/ksahlin/strobemers/blob/main/randstrobe_implementations/figures/clumpings_motivation.pdf) to understand why we measure this.


1. **total_unique**: unique positions position of strobe2 (higher is better - as long as not scenario B i fig)
2. **most_repetitive**: position that was most sampled (may indicate a degenerate region, see scenario C in fig) 
3. **distance_nonuniformity**: Should be close to 1 for perfect uniform (if we want strobes sampled uniform from downstream window)
4. **d_minN**: Observed number of collisions (position sampled twice) divided by expected number of collisions in a group of N consecutive strobes. A value close to 1 is good. >1 indicates bias. Expected collisions derived experimentally using pythons random number generator. 
5. **d_maxN**: Tight groupings. A tight grouping is when the maximum distance between any two strobes within a gropup of N consecutive strobes is N (see scenario E in fig). d_maxN measure observed groupings divided by expected groupings. A value close to 1 is good. >1 indicates bias. 
```

