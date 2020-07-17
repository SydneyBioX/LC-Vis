# LC-N2G


Shiny app(LC-N2G) explores the relationship between nutrition and its corresponding gene expression data. The default dataset comes from the mouse nutrition study (GSE85998)[1]. The overall workflow of LC-N2G is as follows and for a full description, we refer to [our paper](Our-paper).

## Installation

### Packages Requirements

- shiny ≥ 1.4.0.2
- shinythemes ≥ 1.1.2
- shinyjs ≥ 1.1
- DT ≥ 0.13
- directlabels ≥ 2020.1.31
- GA ≥ 3.2
- tidyverse ≥ 1.3.0
- visNetwork ≥ 2.0.9
- dynamicTreeCut ≥ 1.63-1
- fields ≥ 10.3
- reshape2 ≥ 1.4.3
- ggdendro ≥ 0.1-20
- gridExtra ≥ 2.3

### Get started

To use this shiny app you can either:

 - visit our webpage http://shiny.maths.usyd.edu.au/LC-N2G/ , or
 
 - install it through
 
	``` r
	remotes::install_github("SydneyBioX/LCN2G")
	library(LCN2G)
	run_App()
	```
## Vignette

You can find the vignette at our website: https://sydneybiox.github.io/LCN2G/.


## Our paper

Xu, X.N., Solon-Biet, S.M. et al: LC-N2G: A Local Consistency Approach for Nutrigenomics Data Analysis. (Submitted to BMC Bioinformatics) 

## Reference

[1] Solon-Biet, S.M., Cogger, V.C. et al: Defining the nutritional and metabolic context of fgf21 using the geometric framework. Cell Metabolism 24, 555–565 (2016)

[2] Raubenheimer, D., Simpson, S.J.: Nutritional ecology and human health. Annual Review of Nutrition 36, 603–626 (2016)
