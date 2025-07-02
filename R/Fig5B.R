library(Seurat)
library(tidyverse)
library(patchwork)
library(grid)
library(ggpubr)

#Warning, the full Seurat object requires >32GB RAM to load.
AML.merged <- readRDS("data/Lilljebjorn_etal_scRNA-data_AML_NBM_Seurat4.rds")


MHCclass1 <- c("HLA-A","HLA-B","HLA-C")
MHCclass2 <- scan("data/HLA_MHC_Class2.txt", what="character")

AML.merged <- AddModuleScore(AML.merged, features=list(MHCclass1,MHCclass2), name="MHC_class_")

Imcells <- subset(AML.merged, subset = (celltype.aml %in% c("AML Immature")) | ((AML.type %in% c("NBM-MNC","NBM-CD34") ) & celltype.aml %in% c("HSC")) )

Imcells.Meta <- FetchData(Imcells, vars=c("orig.ident","celltype.aml")) %>% rownames_to_column() %>% mutate(
  celltype.aml = factor(x=celltype.aml, levels=c("HSC","AML Immature"))) %>% 
  column_to_rownames()

Imcells <- AddMetaData(Imcells,Imcells.Meta)

p <- VlnPlot(Imcells, features="MHC_class_2", group.by="AML.type.NPM1class",pt.size = 0,split.by="celltype.aml", cols=c("#66b266","#fec44f")) + 
  geom_boxplot(width=0.2,outlier.shape=NA,fill="white") + 
  theme_classic(base_size=5) + 
  theme(legend.key.size = unit(0.4, 'cm'), legend.position = "inside" , legend.position.inside = c(.85,.9), legend.text = element_text(size=8), axis.text.x = element_text(angle=45, hjust=1, size = 8), plot.title = element_blank(), axis.title.x = element_blank(),
        axis.text.y = element_text(size = 8), line = element_line(linewidth = 0.6)) + 
  stat_compare_means(comparisons = list( c("NBM-CD34","NPM1 class I"),c("NPM1 class I","NPM1 class II") ), label = "p.signif", size = 7, fontface = "bold") + 
  scale_y_continuous(limits=c(NA,2), expand=c(NA,0))

ggsave("Figures/Fig5B.pdf",p,width=90,height=97.5, unit="mm")
