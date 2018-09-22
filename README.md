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

MASS, gplots, devtools (installation from GitHub only)
## Guided Tutorial
For this tutorial, we will be predicting cell types from a trimmed single-cell RNA sequencing (scRNA-seq) data matrix freely available from 10X Genomics. The signatures we used as prior knowledge are extracted from published eTME signatures (Wang T, et al., 2018). All the data we need for this tutorial is availiable at [here](https://github.com/jcao89757/SCINA/tree/master/inst/extdata).
### Prepare input data
The SCINA model takes at least two input data matracies to predict categories.
1. A normalized matrix representing the target dataset. Columns correpond to objects (cell barcodes for example), rows correspond to attributes or variables (gene symbols for example). Please find the [.RData example](https://github.com/jcao89757/SCINA/tree/master/inst/extdata/example_expmat.RData) and the [.csv example](https://github.com/jcao89757/SCINA/tree/master/inst/extdata/example_expmat.csv) for the target expression matrix.
2. A list contains multiple signature identifier lists. Each signature identifier list (genes for example) represents prior knowledge for one category (cell type for example), containing genes or protein symbols with high degree of detection. Please find the [.RData example](https://github.com/jcao89757/SCINA/tree/master/inst/extdata/example_signatures.RData) and the [.csv example](https://github.com/jcao89757/SCINA/tree/master/inst/extdata/example_signatures.csv) for the signature lists.

Both matrices can be uploaded from .Rdata files or .csv files. If the target dataset is uploaded with .csv files, the format requirements are the same as the descriptions above **(Fig.1)**. If the signature identifier list is uploaded with .csv files, each column contains one signature list, and its column name should be the name of the category. Each signature identifier list contains gene or protein symbols. The identifier lists do not need to have the same length **(Fig.2)**.

![exp_example](exp_example.jpg)

**Fig.1 |** An example of a target dataset in .csv format.

![exp_signature](exp_signature.jpg)

**Fig.2 |** An example of signature lists in .csv format.

After downloading, please load the matrices to your R environment.
```{r}
#.Rdata examples
load(system.file('extdata','example_expmat.RData', package = "SCINA"))
load(system.file('extdata','example_signatures.RData', package = "SCINA"))
exp = exp_test$exp_data

# Or .csv examples
exp=read.csv('your/path/to/example_expmat.csv',row.names=1,stringsAsFactors = F)
signatures=preprocess.signatures('your/path/to/example_signatures.csv')
```
### Standard pre-processing workflow
The example expression matrix we provided here is a normalized example. In most scenarios, users are encouraged to preprocess raw count outputs of their sequencing data. Considering the features of scRNA-seq data, we suggest that users may follow the pre-processing code below to achieve a best performance on their scRNA-seq raw counts.
```{r}
#Install preprocessCore if required
source("http://bioconductor.org/biocLite.R")
biocLite("preprocessCore")
library('preprocessCore')
#Read data
exp_raw=read.csv('your/path/to/raw/expression_matrix.csv',row.names=1,stringsAsFactors = F)
#Avoid NAs
exp_raw[is.na(exp_raw)]=0
#Log scale and quantile normalization
exp_raw=log(exp_raw+1)
exp=normalize.quantiles(exp_raw)
```
### Set model parameters
The SCINA algorithm has multiple parameters that users may tune to achieve a better performance. The table below contains description of those parameters.

|Parameters|Details|
|----------|-------|
|max_iter|An integer > 0. Default is 100. Max iterations allowed for SCINA algorithm.|
|convergence_n|An integer > 0. Default is 10. Stop SCINA if during the last n rounds of iterations, cell type assignment keeps steady above the convergence_rate.|
|convergence_rate|A float between 0 and 1. Default is 0.99. Percentage of cells for which the type assignment remains stable for the last n rounds.|
|sensitivity_cutoff|A float between 0 and 1. Default is 1. The cutoff to remove signatures whose cells types are deemed as non-existent at all in the data by SCINA.|
|rm_overlap|A binary value, default 1 (TRUE), denotes that shared symbols between signature lists will be removed. If 0 (FALSE) then allows different cell types to share the same identifiers.|
|allow_unknown|A binary value, default 1 (TRUE). If 0 (FALSE) then no cell will be assigned to the 'unknown' category.|

In addition to the parameters, SCINA may generate a log file to record the running status and errors. Please specify a string as the log file's name (and input the string as 'log_file'), path may be included in the string. The default log_file is 'SCINA.log'

### Predict object categories with SCINA
SCINA can generate two output lists in a result class for processed data matrix. 
```{r}
results = SCINA(exp, signatures, max_iter = 100, convergence_n = 10, 
    convergence_rate = 0.999, sensitivity_cutoff = 0.9, rm_overlap=TRUE, allow_unknown=TRUE, log_file='SCINA.log')
View(results$cell_labels)
View(results$probabilities)
```
More detail of outputs is described in the below table.

|Output|Details|
|------|-------|
|cell_labels|A vector contains cell type predictive results for each cell.|
|probabilities|A probability matrix indicating the predicted probability for each cell belongs to each cell type respectively.|
### Result visualization
## Version update
1.0.0: First release. (09-20-2018)
