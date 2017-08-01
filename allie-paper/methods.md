methods
=======

# goal
- we want to be able to identify divergent behavior between samples taken under different conditions
    - *e.g.:*
        - stages of disease/treatment
        - patients with different strains of a disease
- goals:
    1. if a robust predictive model for some divergent behavior is obtained, it can be used for diagnosis
    2. if the model admits some ready biological interpretation, it can be understood and modelled for better diagnosis, and maybe treatment
- difficulty:
    - relationships between sample may be the product of multiple systems, leading to complex and/or multivariate relationships
        - cytof specifically can measure many variables at once at single-cell resolution with constantly increasing accuracy
            - *allows* more complex relationships, but *requires* different tools to analyze
            - sparsity and dimensionality of cytof data makes computation slow
            - obscures well-known relationships such as euclidean distance (hypersphere etc)
            - different markers measured may have wildly varying levels of variance, and typically escape classification by simple parametric models (e.g. gaussian)
        - different transformations of cytof data will typically aim to relax one of these constraints
            - may obscure some relationships, but make others more obvious
            - asinh (sparse) / viSNE (euclid) / MEM (varying variance)
        - really useful transformations admit some strong biological interpretation
            - if so, a model which can robustly predict or distinguish between outcomes in transformed data can then be converted into a testable hypothesis on the original data
                - if this is successful, the transformation becomes extremely useful for generating treatments for diseases

## viSNE / EMD
- a viSNE analysis was performed across the dataset of 47 files
    - output dimensions = 2
    - equal sampling, 5000 events per file
    - iterations = 1000
    - perplexity = 30
    - theta = 0.5
- the earth mover's distance was calculated between each pair of files using the `transport` R library
    - sampled 1000 events out of each 5000 randomly without replacement
        - same 1000 events used for all comparisons
        - **TODO: is this kosher?**
    - convert each `1000x2` matrix into a set of 1000 points in two dimensions, each with weight 1 (`transport::wpp`)
    - calculate the wasserstein distance between each point set using the `"shortsimplex"` method with no special parameters
    - optimizations:
        - the `parallel` library was used to speed the computation
        - each population was automatically given a `0` EMD compared to itself, and the already-computed EMD scores were reused across the diagonal

## heatmap
- given square triangular distance matrix
- scale assigns 100 to pairs of populations with 0 distance between them
    - e.g. population compared to itself
- 0 to the pair of populations with the maximum distance between them
    - among entries in the given distance matrix

1. [x] make emd like mem rmsd
2. [ ] write description of current working version
    - don't forget parallel library
    - "bullet points ok"
    - we're sampling and we don't know how that affects result
3. [ ] try binning instead of sampling and see if we can make something that runs faster without going haywire

### random unrelated links (really)
- https://www.math.hmc.edu/~su/papers.dir/metrics.pdf

## MEM / RMSD
- there's a paper on this so less necessary
