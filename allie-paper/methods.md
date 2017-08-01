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
- a viSNE analysis was performed across the dataset
    - output dimensions = 2
    - equal sampling, 5000 events per file
    - iterations = 1000
    - perplexity = 30
    - theta = 0.5
    - final KL divergence = 5.034414
- the earth mover's distance was calculated between each pair of files using the `transport` R library
    - sampled 1000 events randomly without replacement
        - same 1000 events used for all comparisons
        - **TODO: is this kosher?**
    - either "revsimplex" or "shortsimplex", determine which
        - mention multiscale parameters (`trcontrol()`)
        - the `transport` library's documentation has a detailed explanation of how the multiscale computation is performed
    - the `parallel` library was also used to speed the computation
        - automatically determined number of cores with `parallel::detectCores()`
    - no specific seed was used
    - each population was automatically given a `0` EMD compared to itself, and the already-computed EMD scores were reused across the diagonal


### random unrelated links (really)
- https://www.math.hmc.edu/~su/papers.dir/metrics.pdf

## MEM / RMSD
- there's a paper on this so less necessary
