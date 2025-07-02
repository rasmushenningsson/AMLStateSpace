library(tidyverse)
library(scales)


mut_vaf <- read.csv("data/single_cell_mutation_data/mutations_bulkvaf.csv", row.names = 1)
annotations <- read_csv("data/SingleCell_Metadata.csv.gz")
data = NULL

filenames <- list.files(path = "data/single_cell_mutation_data", pattern = "*_genotype_table.csv", full.names = TRUE)

for (file in filenames) {
  print(file)
  case = str_extract(file,"AML[0-9]*[^_]")
  
  caseLong <- case_when(case == "AML85D" ~ "AML85-D",
                        case == "AML85R" ~ "AML85-R",
                        TRUE ~ case)
  
  print(case)
  print(caseLong)

  
  table <- read.csv(file, row.names = 1)
  colnames(table) <- gsub("_WT","_wt" ,colnames(table))
  colnames(table) <- gsub("_Mut","_mut" ,colnames(table))
  colnames(table) <- gsub(".*\\.","",colnames(table))
  
  mut_table <- mut_vaf %>% filter(AML == case)
  
  variantcases <- gsub("_mut", "", colnames(table %>% select(-cell, -ends_with("_wt"))))
  for (variant in variantcases) {
    print(variant)
    class <- mut_table$Class[mut_table$variant == variant]
    colnames(table)[grep(variant,colnames(table))] <- paste0(class,"_",colnames(table)[grep(variant,colnames(table))])
  }
  
  AnnData <- annotations %>% filter(orig.ident == caseLong) %>% select(cell,celltype.aml)
  table$cell <- paste(caseLong,"_",table$cell,"-1",sep="" )
  
  AnnData <- AnnData %>% inner_join(table)
  #AnnData <- AnnData %>% mutate(all.mutSum = rowSums(select(., ends_with("_mut"))), all.wtSum = rowSums(select(., ends_with("_wt"))), all.mut.status = case_when(all.mutSum > 0 ~ "mut", all.wtSum > 4 ~ "wt", TRUE ~ "unknown"))
  #AnnData <- AnnData %>% mutate(n.mutSum = (all.mutSum - rowSums(select(., matches("^MT.*_mut$")))), n.wtSum = (all.wtSum - rowSums(select(., matches("^MT.*_wt$")))), nuclear.mut.status = case_when(n.mutSum > 0 ~ "mut", n.wtSum > 4 ~ "wt", TRUE ~ "unknown"))
  
  VafTable <- AnnData %>% group_by(celltype.aml) %>% summarise(n=n(),across(contains(c("_")),~sum(.x)))
  
  
  ################################################################
  # Will only use data from heterozygous mutations (ignore others)
  #
  
  heterozygous <- grepl('heterozygous_',colnames(AnnData)) 
  if(any(heterozygous == TRUE)) {  
    
    AnnData <- AnnData %>% 
      mutate(heterozygous.mutSum = rowSums(select(., matches("^heterozygous_.*_mut$"))), 
             heterozygous.wtSum = rowSums(select(., matches("^heterozygous_.*_wt$"))), 
             heterozygous.mut.status = case_when(heterozygous.mutSum > 0 ~ "mut", 
                                                 heterozygous.wtSum > 4 ~ "wt", 
                                                 TRUE ~ "unknown"), 
             heterozygous.transcript = case_when(heterozygous.mutSum > 0 ~ "informative", 
                                                 heterozygous.wtSum > 0 ~ "informative", 
                                                 TRUE ~ "unknown"))   
    # If there is more than one
    
    Nb <-length(grep('heterozygous_',colnames(AnnData)))
    
#    if(Nb == 2){
#      VafTable2 <- AnnData %>% 
#        group_by(celltype.aml) %>% 
#        summarise(n=n(),across(contains(c("_")),~sum(.x)), 
#                  heterozygous.prop=(sum(heterozygous.mutSum)/(sum(heterozygous.mutSum)+sum(heterozygous.wtSum))/0.5), 
#                  heterozygous.prop2= sum(heterozygous.mut.status == "mut")/sum(n), 
#                  heterozygous.prop3= sum(heterozygous.mut.status == "mut")/sum(heterozygous.transcript == "informative"))
#      if(any(VafTable2$heterozygous.prop == "NaN")) { VafTable2$heterozygous.prop2[is.na(VafTable2$heterozygous.prop)] <- NaN } 
#      if(any(VafTable2$heterozygous.prop == "NaN")) { VafTable2$heterozygous.prop3[is.na(VafTable2$heterozygous.prop)] <- NaN } 
#    } else {
      VafTable2 <- AnnData %>% 
        group_by(celltype.aml) %>% 
        summarise(n=n(),across(contains(c("_")),~sum(.x)), 
                  heterozygous.mut=sum(heterozygous.mutSum),
                  heterozygous.wt=sum(heterozygous.wtSum), 
                  heterozygous.prop=(sum(heterozygous.mutSum)/(sum(heterozygous.mutSum)+sum(heterozygous.wtSum))/0.5), 
                  heterozygous.prop2= sum(heterozygous.mut.status == "mut")/sum(n), 
                  heterozygous.prop3= sum(heterozygous.mut.status == "mut")/sum(heterozygous.transcript == "informative"))
      if(any(VafTable2$heterozygous.prop == "NaN")) { VafTable2$heterozygous.prop2[is.na(VafTable2$heterozygous.prop)] <- NaN } 
      if(any(VafTable2$heterozygous.prop == "NaN")) { VafTable2$heterozygous.prop3[is.na(VafTable2$heterozygous.prop)] <- NaN } 
#    }
  
  
    VafTable <- VafTable %>% inner_join(VafTable2)
    for (i in 1:length(VafTable$heterozygous.prop)) {
      if(VafTable$heterozygous.prop[i] == "NaN"){print("NaN exsist")} 
      else if (VafTable$heterozygous.prop[i] > 1) { VafTable$heterozygous.prop[i] <- 1}
      }
    
    # heterozygous.prop = proportion of reads with heterozygous mutation divided by 0.5. 
    # This value is an estimation of heterozygously mutated cells in population based on overall vaf
    #
    # heterozygous.prop2 = proportion of all cells where an actual mutation has been detected.
    # This value represents the minimum number of mutated cells in the population based on the available data.
    #
    # heterozygous.prop3 = proportion of cells with read information where an actual mutation has been detected.
    # This value represents another estimate of the proportion of mutated cells in the population, but this is probably
    # less reliable than both the values above
    #
    
    VafTable$heterozygous.prop <- round(VafTable$heterozygous.prop, 2)
    VafTable$heterozygous.prop2 <- round(VafTable$heterozygous.prop2, 2)
    VafTable$heterozygous.prop3 <- round(VafTable$heterozygous.prop3, 2)
    
    # structure table 
    #hetero <- colnames(VafTable %>% select(starts_with("heterozygous")) )
    #VafTable <- VafTable %>% relocate(all_of(hetero), .after = last_col())
    
    VafShort <- VafTable %>% 
      select(celltype.aml,n,heterozygous.prop,heterozygous.prop2,heterozygous.mut,heterozygous.wt) %>% 
      rename(c("celltype" = "celltype.aml", "ncell" = "n", "PredHet" = "heterozygous.prop", "PropAll" = "heterozygous.prop2")) %>% 
      rowwise() %>%
      mutate(AML = max(c(PredHet,PropAll)),
             Case = case)
    
    data <- rbind(data,VafShort)

    }
}


# Prepare the data:

cellevels <- c("AML Immature", "HSC", "LMPP", "GMP", "Monocytes", "Megakaryocytic cells", "Erythroid", 
               "Dendritic cells", "NK-cells", "T-cells", "B-cells")

I2order <- c("AML7", "AML9","AML24", "AML25", "AML27", 
             "AML28", "AML33", "AML79","AML97", "AML104", "AML105","AML136",  
             "AML10", "AML21", "AML32", "AML55", "AML110", 
             "AML111", "AML126", "AML138", "AML151", "AML157", "AML172", "AML34", 
             "AML37", "AML48", "AML62", "AML80", "AML83", "AML85D", "AML85R", 
             "AML155", "AML61", "AML117", "AML161", "AML23", "AML123", "AML124")


NPM1 <- c("AML24","AML25","AML27","AML28","AML79","AML104","AML136","AML7","AML9","AML33","AML97","AML105")
CScases <- c("AML21","AML32", "AML10", "AML110", "AML126", "AML151", "AML55", "AML111",  "AML157", "AML172", "AML138")
TP53 <- c("AML37","AML48","AML62", "AML80", "AML83", "AML85D", "AML85R")
Other <- c("AML34", "AML155")
RUNX1f <- c("AML61","AML117","AML161")
CBFBf <- c("AML23","AML123","AML124")



data <- as_tibble(data) %>% 
  mutate(celltype = base::factor(celltype, levels = cellevels), 
         Case = base::factor(Case,levels=rev(I2order)),
         group = case_when(Case %in% NPM1 ~ "NPM1",
                            Case %in% CScases ~ "AML-MR",
                            Case %in% TP53 ~ "TP53",
                            Case %in% RUNX1f ~ "RUNX1::\nRUNX1T1",
                            Case %in% CBFBf ~ "CBFB::\nMYH11",
                            TRUE ~ "Other"),
         group = factor(group,levels = c("NPM1","AML-MR","TP53","Other","CBFB::\nMYH11","RUNX1::\nRUNX1T1")))


data <- data %>% mutate(AMLprop = case_when( (`heterozygous.mut` + `heterozygous.wt`) >3 ~ AML, TRUE ~ NA_real_))


# Filter out AML37, AML62 and AML85D since it is unclear if their mutations are really heterozygous, and AML151 which is lacking data 

data <- data %>% filter(! (Case %in% c("AML37","AML62","AML85D", "AML151")) )


# Plot the heatmap:

hmap <- 
  ggplot(data,aes(celltype,Case,fill=AMLprop,color='NA')) + geom_tile(color="black") + 
  theme_classic(base_size = 5) + 
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5), 
        strip.text.y = element_text(size = 3),
        legend.position = "none",
        axis.title.x = element_blank()) +
  facet_grid("group",scales="free",space="free") + 
  scale_fill_continuous(labels=label_percent(), na.value="#A08968") +
  labs(fill="Cells with\nAML-mutations", color = "", x="Celltype") +
  #scale_shape_manual(values = c("Not enough data" = 0, "No cells" = 0)) +
  #  scale_fill_manual(values = c("Not enough data" = "#A08968", "No cells" = "#fffff")) +
  scale_x_discrete(labels=c("Erythroid" = "Erythroid cells"))

# Save the heatmap:

ggsave("Figures/Fig3D.pdf",hmap,width=65,height = 85, units = "mm")
