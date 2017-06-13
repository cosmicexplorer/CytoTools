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
- run visne on all files from all patients (output should be a visne plot per file)
    - starting pop is called "live leukocytes"
    - on all protein markers
        - not numeric ones
        - run this on *markers shared between all pops*
            - and maybe fix some small inconsistencies in column naming
    - with all cells as ref pop for visne
        - so each file's visne is just its subset of that total visne
- install and run pairwise emd on all visne plots (this is an R script)
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

# part 3
- run MEM on specified populations
    - with "all cells all patients" as ref pop
    - run it with same scale for all patients so they're comparable

# part 4???
- rmsd comparison of MEM is challenging
    - need to have same markers
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
