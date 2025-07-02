library(Seurat)
library(tidyverse)
library(patchwork)

#Warning, the full Seurat object requires >32GB RAM to load.
AML.merged <- readRDS("data/Lilljebjorn_etal_scRNA-data_AML_NBM_Seurat4.rds")

clusterplot <- DimPlot(AML.merged, group.by="seurat_clusters", label=TRUE, raster=FALSE) + 
  coord_fixed() + theme_void() + theme(legend.position="none", plot.title = element_blank())

cellcolors_orig <- c("HSC"="#66b266", 
                "Monocytes"="#fa9fb5", 
                "GMP"="#008000", 
                "Megakaryocytic cells"="#e5687e",
                "LMPP"="#756bb1", 
                "Erythroid"= "#a7241d",
                "B-cells"= "#a64ca6",
                "T-cells"="#7fb4e5",
                "NK-cells"="#196fbe",
                "Dendritic cells"= "#ff4d00")

celltypeplot <- DimPlot(AML.merged, group.by="celltype.l4", label=FALSE, raster=FALSE, cols = cellcolors_orig) + 
  coord_fixed() + theme_void() + theme(legend.position="bottom", plot.title = element_blank())

AMLcolors <- c("NBM-CD34" = "gray60",
               "NBM-MNC" = "gray80", 
               "NPM1"="#377eb8",
               "AML-MR"="#4daf4a",
               "TP53"="#e41a1c",
               "In two subgroups" = "black", 
               "No class-defining"= "#fa9fb5", 
               "CBFB::MYH11" = "#984ea3", 
               "RUNX1::RUNX1T1" = "#ff7f00")


AMLtypeplot <- DimPlot(AML.merged, group.by="AML.type", label=FALSE, raster=FALSE, cols = AMLcolors) + 
  coord_fixed() + theme_void() + theme(legend.position="bottom", plot.title = element_blank())

Fig2B <- clusterplot | celltypeplot | AMLtypeplot

options(bitmapType = 'cairo')

ggsave("Figures/Fig2B.pdf", Fig2B,height=8, width=21)