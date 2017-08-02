methods
=======

# outline

## citations
- EMD and usage in bioinformatics: http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0151859
- prior work using EMD on viSNE axes: http://onlinelibrary.wiley.com/doi/10.1002/cyto.a.22906/full
- the R transport package: https://cran.r-project.org/web/packages/transport/citation.html
- the shortsimplex method (authors of transport package wrote this paper): http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0110214
- the MEM paper describes MEM and MEM RMSD, obviously

## written out

### Pairwise MEM RMSD Calculation
The MEM vectors for each non-reference population were calculated over 25 phenotype channels which were shared across all 140 non-reference populations and the single reference population. Each MEM vector contained the population's MEM score, calculated for each of the 25 common phenotype channels, in reference to the single reference population. The MEM RMSD between pairs of non-reference populations was then calculated using the Euclidean distance between these MEM vectors.

The 25 phenotype channels used in the MEM RMSD analysis were: [ICOS, CD19, TIM3, CCR5, CD4, CD20, CCR4, TCRGD, CD45RA, CD45, CXCR3, CCR7, CD28, CD69, HLA-DR, CD45RO, CD44, CD3, CXCR5, PD-1, CD56, CD16, CD38, CD8, CD27]. Some channel names were spelled differently across the files in the dataset, so all channel names were normalized using R before the MEM calculation. Each of the 25 channels or a variant appeared exactly once as a column name in each population, including the reference.

### Pairwise EMD Calculation
The Earth Mover's Distance (EMD) was calculated between each pair of populations using the "transport" library for R.  The dataset of 47 fcs files used for this analysis was uploaded to Cytobank. A viSNE analysis with two output dimensions was performed, equally sampling 5000 events per file, with 1000 iterations, perplexity equal to 30, and theta equal to 0.5.

The events with their viSNE axes were then downloaded from Cytobank, and the Earth Mover's Distance (EMD) was calculated between each pair of files using the "transport" library for R. The "wpp" object was used to represent each set of points in the two viSNE axes, and the "wasserstein" function was called on each pair of point sets to produce a distance matrix. Each point was assigned unit weight.

Because calculating a matrix with the EMD between each set of 5000 events from the viSNE analysis is computationally expensive, four optimizations were performed.

1. Each file was further downsampled to 1000 out of the original 5000 events per file in the viSNE analysis. Each event was still assigned unit weight, and each point set therefore still had an equal total mass of 1000. This did not affect the expected value of the EMD calculation because each subset of 1000 points per file is equally likely to be picked.
2. The "shortsimplex" method was used for the "wasserstein" function in the "transport" library, which accepted no other parameters besides the pair of weighted point sets.
3. Each population was automatically assigned a zero EMD compared to itself, and EMD scores already computed across the diagonal were simply copied, since EMD is a metric.
4. The "parallel" library was used to parallelize the computation of each row of the matrix in addition to the above, using the number of cores detected from the "detectCores" function in the "parallel" library.

### Heatmap Generation
Heatmaps representating population similarity were generated from each distance matrix using the "heatmap.2" function of the "gplots" library for R. The distance matrix was normalized by the maximum non-normalized distance d_max between any pair of populations, then multiplied by 100, then subtracted from 100. So for a square distance matrix D with entries d_{i,j}:

1. d_{i,j} >= 0
2. d_max = max_{i,j}(d_{i,j})
3. similarity_{i,j} = 100 - (d_{i,j} / d_max * 100)

The result was that zero entries in the original distance matrix would receive a similarity score of 100, while the pair of populations with greatest distance in the original distance matrix would receive a similarity score of 0.

## goal
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

### viSNE / EMD
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
    - convert each `1000x2` matrix into a set of 1000 discrete points in two dimensions, each with weight 1 (`transport::wpp`)
    - calculate the wasserstein distance between each point set using the `"shortsimplex"` method with no special parameters
    - optimizations:
        - the `parallel` library was used to speed the computation
        - each population was automatically given a `0` EMD compared to itself, and the already-computed EMD scores were reused across the diagonal

### heatmap
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

### MEM / RMSD
- there's a paper on this so less necessary
