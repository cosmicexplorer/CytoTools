cytobank-api
============

right now these are just links on how to use the cytobank api / cytobank to accomplish the goals listed in the [README](README.md)

# what is happening in cytobank

- [parse cytobank gate/pop xml and split fcs into those sections](https://bioconductor.org/packages/devel/bioc/vignettes/CytoML/inst/doc/HowToParseGatingML.html#01_use_cytobank2gatingset)
- [description of "% in gate" stat](https://support.cytobank.org/hc/en-us/articles/206061197-Important-notes-on-the-Percent-in-Gate-statistic)
    - i just think this sounds like an interesting thing to be able to show
- [description of how to run a (manual) visne analysis](https://support.cytobank.org/hc/en-us/articles/206439707-How-to-Configure-and-Run-a-viSNE-Analysis)
- [analysis of visne walkthrough](https://support.cytobank.org/hc/en-us/articles/223521168-Analysis-and-Interpretation-of-viSNE-Results)
    - ^this has links showing how to make a heatmap on a subset of channels
        - [how to choose channels](https://support.cytobank.org/hc/en-us/articles/206147637#choosing_channels)
        - [how to scale heatmap](https://support.cytobank.org/hc/en-us/articles/206147637-How-to-create-and-configure-a-Heatmap#practical-example-data-scaling-heatmaps)

# how to use cytobank api

## REST api
- [http api ref](https://developer.cytobank.org/api-endpoint-reference.html)

## r language api
- [quickstart](https://cran.r-project.org/web/packages/CytobankAPI/vignettes/cytobank-quickstart.html)
- ["advanced" analyses](https://cran.r-project.org/web/packages/CytobankAPI/vignettes/cytobank-advanced-analyses.html)
