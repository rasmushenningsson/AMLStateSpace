# R plots
The plots generated in R (https://www.r-project.org/) were produced using either ggplot2 version 3.5.1 (https://ggplot2.tidyverse.org/) or Seurat version 4.2.0 (https://satijalab.org/seurat/). 

## Installation
First, R needs to be installed. Follow the instructions for your platform at https://www.r-project.org/.

Then, start R and install the required packages:

```r
> install.packages(c("ggplot2", "Seurat", "tidyverse", "patchwork", "gridExtra", "ggpubr", "cowplot", "scales", "RColorBrewer", "survminer", "survival"))
```

## Running
With R running in the current folder, the figures can be created by running:

```r
> source("Fig1A.R")
> source("Fig2B.R")
> source("Fig3A.R")
> source("Fig3C.R")
> source("Fig3D.R")
> source("Fig4A.R")
> source("Fig4B.R")
> source("Fig4C.R")
> source("Fig5A.R")
> source("Fig5B.R")
```

Note, the scripts Fig2B.R, Fig3A.R, Fig4B.R, and Fig5B.R require the file **Lilljebjorn\_etal\_scRNA-data\_AML\_NBM_Seurat4.rds**. This file can be downloaded from https://doi.org/10.17044/scilifelab.23715648 and then placed in the **data/** directory. When the scripts load this file into memory >32 GB of RAM is required.