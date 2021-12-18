Evaluation of methods for constructing randstrobes
===========

This evaluation will compare

For documentation, I will implement a small bake-off as soon as I have time, which will include the following

### Hashing

1. No hash - simple 2bit encoding of nucleotides
2. Thomas Wang hash
3. xxhash3 (64bit)

Let the hash function be denoted by h.

### Linking

1. Sahlin ( h(k_1) + h(k') ) % p (method1 above)
2. Shen ( (h(k_1) + h(k') ) & p (method2 above)
3. Sahlin bitcount( h(k_1) ^ h(k') ) (method3 above)
4. Guo-Pibri h(k_1) ^ h(k') (global XOR), (also a variant: h(k_1 ^ k'))
5. Liu-Patro-Li h( k_1 || k' ) (concatenation)

Viable combinations of (hashing, linking) seem to be (1,1)-(1,4), (2,1)-(2,4), (3,1)-(3,5) as the total length can be larger than 32.


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
./randstrobe_benchmark -x [1,2,3] -l [link method 1-5] refs.fasta queries.fasta
```

If `-l 5` is specified, `-x` will use xxhash128.
