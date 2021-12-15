Evaluation of methods for constructing randstrobes
===========

This evaluation will compare

For documentation, I will implement a small bake-off as soon as I have time, which will include the following

### Hashing

1. No hash - simple 2bit encoding of nucleotides
2. Thomas Wang hash
3. xxhash3

Let the hash function be denoted by h.

### Linking

1. Sahlin ( h(k_1) + h(k') ) % p (method1 above)
2. Shen ( (h(k_1) + h(k') ) & p (method2 above)
3. Sahlin bitcount( h(k_1) ^ h(k') ) (method3 above)
4. Pibri h(k_1) ^ h(k') (global XOR), (also a variant: h(k_1 ^ k'))
5. Liu-Patro-Li h( k_1 || k' ) (concatenation)

Viable combinations of (hashing, linking) seem to be (1,1)-(1,4), (2,1)-(2,4), (3,1)-(3,5) as the total length can be larger than 32.


### Metrics

1. Construction time
2. Uniqueness:
    - Distribution of the number of times the position p2 is selected for k_2 (repetition distribution)
    - Fraction of unique positions p2. 
    - Downstream metrics such as how it affects strobealign (accuracy, final runtime, mapping sites tried, calls to ksw2)

