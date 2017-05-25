allie-paper
===========

- 4 parts total, we'll just do first three now
- review methodology and script some stuff that's not
    - also modify script to run off of binary fcs, not just csv

# part 1
- run visne on all files from all patients (output should be a visne plot per file)
    - starting pop is called "live leukocytes"
    - on all protein markers
        - not numeric ones
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
        - illustration at https://irishlab.cytobank.org/cytobank/experiments/22899/
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
        - and maybe fix some small idiosyncracies in column namingneed
- more info later
