library(Seurat)
library(tidyverse)
library(RColorBrewer)

#Warning, the full Seurat object requires >32GB RAM to load.
AML.merged <- readRDS("data/Lilljebjorn_etal_scRNA-data_AML_NBM_Seurat4.rds")

C3a <- c("AML24","AML25","AML27","AML28","AML79","AML104")
C3b <- c("AML7","AML9","AML33","AML97","AML105")
NBMMNC <- c("NBM4-MNC", "NBM7-MNC", "NBM8-MNC", "NBM10-MNC", "NBM11-MNC") # 5
NBMCD34 <- c("NBM8-CD34","NBM10-CD34","NBM11-CD34") # 3

NPM1.plotcells <- subset(AML.merged, subset = (orig.ident %in% c(NBMMNC,NBMCD34) | (orig.ident %in% c(C3a,C3b) & celltype.aml == "AML Immature")))

plotcells.Meta <- FetchData(NPM1.plotcells, vars = c("orig.ident","celltype.aml")) %>% rownames_to_column() %>% 
  mutate(sampletype = case_when(
      orig.ident %in% NBMMNC ~ "NBM-MNC",
      orig.ident %in% NBMCD34 ~ "NBM-CD34",
      orig.ident %in% c(C3a,C3b) ~ as.character(orig.ident)),
    sampletype = factor(sampletype, levels = c("NBM-CD34","NBM-MNC",C3a,C3b) )
  ) %>% column_to_rownames()

NPM1.plotcells <- AddMetaData(NPM1.plotcells,plotcells.Meta)

custom.palette <- c("gray50","gray80",brewer.pal(6, "Blues"),brewer.pal(5, "Reds"))
Fig4b <- DimPlot(NPM1.plotcells, group.by="sampletype", label=FALSE, raster=FALSE, cols=custom.palette) + 
  coord_fixed() + theme_void() + theme(legend.position="right", plot.title = element_blank())

ggsave("Figures/Fig4B.pdf", Fig4b,height=8, width=7)