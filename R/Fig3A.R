library(Seurat)
library(tidyverse)

#Warning, the full Seurat object requires >32GB RAM to load.  
AML.merged <- readRDS("data/Lilljebjorn_etal_scRNA-data_AML_NBM_Seurat4.rds")


n_cells <- FetchData(AML.merged, 
                     vars = c("orig.ident", "seurat_clusters")) %>%
  group_by(orig.ident) %>%
  dplyr::count(seurat_clusters) %>% 
  spread(seurat_clusters, n)

n_cells <- n_cells %>% mutate(sampletype = case_when(
  str_detect(orig.ident,"MNC") ~ "NBM.MNC",
  str_detect(orig.ident,"CD34") ~ "NBM.CD34",
  str_detect(orig.ident,"AML") ~ "AML"))

n_sum <- n_cells %>% group_by(sampletype) %>% summarise(across(matches("\\d$"),sum, na.rm = TRUE)) %>% 
  mutate(sum = rowSums(across(where(is.numeric)))) %>% 
  column_to_rownames("sampletype")

cluster_pct <- as_tibble(t(n_sum), rownames = "rowname") %>% transmute(
  AMLpc = AML/(AML+n_sum$sum[1]/n_sum$sum[2]*NBM.CD34+n_sum$sum[1]/n_sum$sum[3]*NBM.MNC), 
  NBM.CD34pc = NBM.CD34/(NBM.CD34+n_sum$sum[2]/n_sum$sum[1]*AML+n_sum$sum[2]/n_sum$sum[3]*NBM.MNC),
  NBM.MNCpc = NBM.MNC/(NBM.MNC+n_sum$sum[3]/n_sum$sum[2]*NBM.CD34+n_sum$sum[3]/n_sum$sum[1]*AML),
  rowname = str_replace(rowname,"X","SCT.")) 

n_celltype_cells <- FetchData(AML.merged, 
                              vars = c("celltype.l4", "seurat_clusters")) %>%
  group_by(celltype.l4) %>%
  dplyr::count(seurat_clusters) %>% 
  spread(seurat_clusters, n)

n_celltype_pct <- n_celltype_cells %>% ungroup() %>% mutate(across(where(is.numeric), ~ replace_na(.,0))) %>% 
  summarise(across(where(is.numeric), ~ ./sum(.)*100)) %>% 
  mutate(celltype = n_celltype_cells$celltype.l4, .before=`0`)

tn_celltype_pct <- as.data.frame(t(n_celltype_pct %>% column_to_rownames("celltype"))) %>% rownames_to_column()

cluster_pct <- cluster_pct %>% right_join(tn_celltype_pct, by="rowname")


AMLIm <- cluster_pct %>% filter(AMLpc > 0.7, HSC + LMPP > 40) %>% mutate(rowname = as.numeric(rowname)) %>% pull(rowname)

new.anno <- FetchData(AML.merged, vars = c("orig.ident2", "seurat_clusters", "celltype.l4")) %>%   rownames_to_column() %>% 
  mutate(celltype.aml = case_when(seurat_clusters %in% c(AMLIm) & !grepl("NBM",orig.ident2) ~ "AML Immature", 
                                  TRUE ~ celltype.l4)) %>% 
  column_to_rownames()

cellorder <- c("AML Immature","HSC","LMPP","GMP","Monocytes","Megakaryocytic cells","Erythroid","Dendritic cells","NK-cells","T-cells", "B-cells")
new.anno$celltype.aml <- factor(new.anno$celltype.aml, levels = cellorder)
new.anno$celltype.l4 <- factor(new.anno$celltype.l4, levels = cellorder) %>% droplevels()
AML.merged <- AddMetaData(AML.merged,new.anno)

cellcolors <- c("HSC"="#66b266", 
                "Monocytes"="#fa9fb5", 
                "GMP"="#008000", 
                "Megakaryocytic cells"="#e5687e",
                "LMPP"="#756bb1", 
                "Erythroid"= "#a7241d",
                "B-cells"= "#a64ca6",
                "T-cells"="#7fb4e5",
                "NK-cells"="#196fbe",
                "Dendritic cells"= "#ff4d00",
                "AML Immature" = "#fec44f")

UMAPplot <- DimPlot(AML.merged, group.by="celltype.aml", label=FALSE, raster=FALSE, cols = cellcolors) + 
  coord_fixed() + theme_void() + theme(legend.position="right", plot.title = element_blank())

options(bitmapType = 'cairo')

ggsave("Figures/Fig3A.pdf", UMAPplot,height=8, width=7)