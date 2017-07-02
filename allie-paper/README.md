allie-paper
===========

- **[EXPERIMENT LINK HERE](https://irishlab.cytobank.org/cytobank/experiments/22899)**
    - experiment id **22899**
- 4 parts total, we'll just do first three now
- review methodology and script some stuff that's not
    - also modify script to run off of binary fcs, not just csv
        - **DONE**
    - also fix MEM auto-scaling, if that's not done / easy
        - kirsten *did* add code to set a scale now
            - check this out
            - want to be able to reliably compare MEM across experiments
                - so just: don't scale it
        - try to understand why the scaling was done like that in the first place?

# part 1
- **NO POPULATIONS CONSIDERED IN THIS PART**
- run visne on all files from all patients
    - **for each file, make a new one with visne analysis id which has tsne axes for each event from that file**
    - starting pop is called "live leukocytes"
    - on all protein markers
        - not numeric ones
    - using some existing visne analysis in some existing experiment on cytobank
        - so each file's t-sne axes are just for the events in that file from that total visne
- install and run pairwise emd on all of the resulting files with tsne axes (this is an R script)
    - *emd compares events, with each file as its own "cluster", using only t-SNE axes for distance*
        - emd script does this already (produces the heatmap)
        - and maybe fix some small inconsistencies in column naming
    - this will produce a csv
    - heatmap created from csv (heatmap creation is also an R script)
        - requires manual editing of csv first
        - fix this
        - there is some excel code to do this
    - fix the emd script to allow not running it on cmdline
    - see allie's email "data analysis stuff" for scripts/excel file

# part 2
- "quantifying change in population frequency"
    - freq of subpops in parent pop
    - change of that freq
- getting freq is done in cytobank
    - e.g. freq NK cells out of all peripheral blood is a desired query, but requires some manual work with the cytobank interface
    - figure out how cytobank lets you export gates
        - would be great if it could have that info in the fcs file
        - but likely we'll only be able to download gates as a selection of rows from parent fcs
        - and identify them in parent by `event#`
    - see `Population frequency_populations of interest.xlsx` in Downloads for subpops and parent pops to download from illustration
        - [illustration here](https://irishlab.cytobank.org/cytobank/experiments/22899/illustrations/52053)
        - also check if there are any duplicated event#s
- quantifying change
    - sometimes populations go to 0, so currently we add a little offset so numbers don't blow up
        - see if there's an alternative, this is probs fine
- line graph of change as well
- result of getPopStats magic is like: `list(df_pre, df_3wk, df_12wk, df_6m)`
    - get pair of population and some ancestor of it (e.g. B_cells  -> Live Leukocytes) which has greatest change from pre -> 6m (meaning just abs or pos value)
    - get line graph of relative concentrations for each `pop -> immediate parent` pair over x = (pre, 3wk, 12wk, 6m) time span
        - by patient, also median of all

# part 3
- *same thing as part 1 but with MEM instead of EMD*
- run MEM on specified populations
    - with "all cells all patients" as ref pop
    - run it with same scale for all patients so they're comparable

# part 4???
- rmsd comparison of MEM is challenging
    - need to have same markers
        - **RMSD ON ALL SHARED MARKERS, OR SHARED MARKERS IN EACH PAIRWISE COMPARISON??**
    - making script to produce copies of all input files but only with channels that are shared
        - and maybe fix some small inconsistencies in column naming
- more info later

# TODO
- can/should we clean cyto data?
    - would do this in `read_clean_fcs` in [fcs.R](R/fcs.R)
    - ideas:
        - `stats::na.omit()`
        - `dplyr::distinct()`
    - when would these happen in the data?
        - **find out what these mean instead of just blindly removing them**
        - if they point to some error or inconsistency in dataset, removal may not make sense!
    - allie says nobody does any cleaning to their fcs files, so may not make sense
