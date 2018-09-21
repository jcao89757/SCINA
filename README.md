![SCINA_logo](QBRC.jpg)
# SCINA
A semi-supervised category identification and assignment tool.
## Introduction
SCINA is an automatic cell type detection and assignment algorithm for single cell RNA-Seq (scRNA-seq) and Cytof/FACS data. SCINA is capable of assigning cell type identities to a pool of cells profiled by scRNA-Seq or Cytof/FACS data with prior knowledge of identifiers, such as genes and protein symbols, that are highly-expressed (or lowly-expressed) in each category. 

Please refer to our paper for more details of SCINA:
["SCINA: Semi-Supervised Typing of Single Cells *in silico*"](url pending) Zhang Z, Luo D, et al., 2018

Or please check our web server to run SCINA on the clould: http://lce.biohpc.swmed.edu/scina 

For more details on eTME signatures please refer to our paper:
["An Empirical Approach Leveraging Tumorgrafts to Dissect the Tumor Microenvironment in Renal Cell Carcinoma Identifies Missing Link to Prognostic Inflammatory Factors."](http://cancerdiscovery.aacrjournals.org/content/early/2018/06/08/2159-8290.CD-17-1246) Wang T, Lu R, et al., 2018

## Getting started with SCINA
The SCINA algorithm is constructed with R, and is originally developed as an R package. Users who are not famaliar with basic R programming are suggested to use our web server to run SCINA with a user-friendly GUI. 
### Installation Instructions
#### Install from CRAN
```{r}
install.packages('SCINA')
library('SCINA')
```
#### Install from GitHub
```{r}
library('devtools')
install_github('jcao89757/SCINA')
library('SCINA')
```
### Dependencies
R (version 2.15.0 or later)

**R Packages**

MASS, gplots, devtools(installation from GitHub only)
## Guided Tutorial
### Prepare input data
### Standard pre-processing workflow
### Set model parameters
### Predict object categories with SCINA
### Result visualization
## Version update
1.0.0: First release. (09-20-2018)
