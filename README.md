# PhenoExploreR
PhenoExploreR

This package implements an exploratory model that single cells might use to adapt to microenvironment stressors. 

The approach and exploration of the [IVY Glioblastoma Atlas](http://glioblastoma.alleninstitute.org/) data using this package, 
and is described in O. Celiku, M. Gilbert, & O. Lavi (2019).

## Installation

```r
# Install the current version from GitHub:
# install.packages("devtools")
devtools::install_github("oricel/PhenoExploreR")
```

## Implementation

Source code is available in the R directory. 

 * PhenoExploreRClasses contains class specifications
 * PhenoExploreRDistances contains functions related to pathway distribution distance and global distribution distance
 * PhenoExploreRWalkers contains implementation of the intrinsic walker ("wander") and exploratory walker ("explore")
 * PhenoExploreRPlots contains functions for graphical inspection of the results.
 
## Data

The data are available in the data directory. 

 * IVY.rda is an object of "phenotype" class, containing expr, design, motif activity, pathways, locations, colors attributes.
 * IVY_expr.tsv is a text version of the expression.
 * IVY_full.expr.tsv is a text version of the IVY expression before selection of the gene list of interest.
 * IVY_sample_info.tsv contains sample information. More information is available in the [IVY Glioblastoma Atlas](http://glioblastoma.alleninstitute.org/) site.
 * TCGA_GBM_vs_LGG_FC.tsv contains the fold changes of the genes as computed by differential expression analysis between 
 the TCGA GBMs and LGGs.
 * c2.cp.kegg.v6.1.symbols.gmt contains the KEGG ontology as used in this study; these were obtained from [MSigDB](http://software.broadinstitute.org/gsea/msigdb/index.jsp)
 
## Examples

Examples are given in the examples directory:

  * IVY_overview.Rmd shows how to explore the initial data; the result should look as the corresponding .html
  * IVY_runs.Rmd shows how to run intrinsic and exploratory simulations of the behavior; the 
  results of running this file should be saved in the same directory as intrinsic_run_example.RDS and exploratory_run_example.RDS; 
  these are not included in the directory due to the large size.
  * IVY_runs_inspection.Rmd shows how to inspect and plot various aspects of the simulated data; the results 
  should look similar to the corresponding .html, although the results may differ as new simulations produce somewhat different results.
