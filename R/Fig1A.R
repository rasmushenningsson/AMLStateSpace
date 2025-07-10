library(tidyverse)
library(gridExtra)
library(cowplot)
library(patchwork)

Sample_data <- read.csv("data/sample_data.csv")

# Main Landscape plot

# All AML cases, facet by Papaemmanuil group.

gene_order_long <- c("NPM1", "FLT3.ITD", "DNMT3A", "TET2", "IDH2", "IDH1", "WT1", 
                         "FLT3", "SPTBN2", "SRCAP", "PDS5B", "KRAS", "MED15", "ELMSAN1", "TBL1XR1", 
                         "RUNX1", "SRSF2", "NRAS", "ASXL1", "STAG2", "BCOR", "CEBPA", "U2AF1", "PTPN11",
                         "ZRSR2", "SMC1A", "BCORL1", "CBL", "EZH2", "GRM4", "PHF6", "ACAN", "ATM", 
                         "MED12", "SF3B1","TP53", "NF1", "TAS1R2",   "NCOR1", 
                         "JAK2", "APOB", "UNC79", "UNC13C", "TENM3", "SLC26A3", "RAD21", 
                         "PTPRB",  "MECOM", "KMT2C", "KIT",  "GATA2", "FHOD3", 
                         "ERC2", "DNAH17", "DDX41", "CSMD3", "BPTF", "ZFHX4", "ZC3H4", 
                         "SMC3", "SETD2", "SETBP1", "PCSK5", "ETV6", "DNAH12", "DHX29", 
                         "CSMD1", "CACNA1E", "ASXL2", "AHCTF1")

gene_order_short <- c("NPM1", "FLT3.ITD", "DNMT3A", "TET2", "IDH2", "IDH1", "WT1", 
                          "FLT3",  "SRCAP", "PDS5B", "RUNX1", "SRSF2", "NRAS", "ASXL1", 
                          "STAG2", "BCOR", "CEBPA", "U2AF1", "PTPN11", "ZRSR2", "SMC1A", 
                          "BCORL1", "TP53", "NF1", "JAK2")


n <- length(Sample_data$SamplingEventID)

Sampleinfo  <- paste("Samples (n=",n,")",sep="")
lines <- seq(from=5.5,to=25.5, by=5)
grouplevels <- c("NPM1","AML-MR","TP53","In two subgroups","No class-defining","No detected driver","biallelic CEBPA","CBFB::MYH11","PML::RARA","KMT2A fusion", "RUNX1::RUNX1T1","DEK::NUP214","GATA2::MECOM")

#############################
# Mutation part of the plot

darkcells1 <-  Sample_data %>% 
  arrange(DNMT3A,NPM1,FLT3.ITD,TET2,RUNX1,SRSF2,TP53,ASXL1,IDH2,NRAS,IDH1,STAG2,SamplingEventID) %>% 
  mutate(SamplingEventID = factor(SamplingEventID,unique(SamplingEventID))) %>% 
  pivot_longer(cols = str_replace(gene_order_long,"FLT3_ITD","FLT3.ITD"), names_to = "Gene",values_to="status") %>% 
  select(SamplingEventID,Gene,status,PapaemClass_short) %>% 
  filter(Gene %in% c("IDH1","WT1","FLT3","SRCAP","PDS5B","BCOR","CEBPA","U2AF1","PTPN11","ZRSR2") & is.na(status)) %>% 
  mutate(status = "NA")

p1 <- Sample_data %>% 
  arrange(DNMT3A,NPM1,FLT3.ITD,TET2,RUNX1,SRSF2,TP53,ASXL1,IDH2,NRAS,IDH1,STAG2,SamplingEventID) %>% 
  mutate(SamplingEventID = factor(SamplingEventID,unique(SamplingEventID))) %>% 
  pivot_longer(cols = str_replace(gene_order_long,"FLT3_ITD","FLT3.ITD"), names_to = "Gene",values_to="status") %>% 
  select(SamplingEventID,Gene,status,PapaemClass_short) %>% 
  filter(Gene %in% gene_order_short) %>% 
  ggplot(aes(SamplingEventID,factor(Gene,levels = rev(gene_order_short)))) + 
  geom_tile(aes(fill=status),colour = "white") + 
  facet_grid(~factor(PapaemClass_short, levels = grouplevels),scales="free_x",space="free") + 
  theme_gray(base_size = 6) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none", 
        axis.ticks = element_blank(), 
        axis.title.x = element_blank(), 
        strip.text.x = element_text(angle=90), 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(colour = "black", face = "italic"),
        panel.spacing = unit(1, "pt"),
        plot.margin = margin(t = 1.5,r = 1.5,b=0, l=1.5,unit = "pt") ) +
  scale_fill_manual(values=c("CHET" = "#7b3294","HET" = "#2b83ba","HOM" = "#d7191c", "NA" = "gray80"),na.value="transparent") + 
  #  scale_fill_manual(values=c("#7b3294","#2b83ba","#d7191c"),na.value=c("transparent","gray80","black")) + 
  labs(y="Mutated genes") + 
  #geom_hline(yintercept=lines, color="white") + 
  scale_y_discrete(labels=c("FLT3.ITD" = "FLT3-ITD")) +
  geom_tile(aes(fill=status), color="white", data = darkcells1)

p1grob <- ggplotGrob(p1)

#wrap_plots(p1grob)

# Get a usable label:

p1Label <- Sample_data %>% 
  arrange(DNMT3A,NPM1,FLT3.ITD,TET2,RUNX1,SRSF2,TP53,ASXL1,IDH2,NRAS,IDH1,STAG2,SamplingEventID) %>% 
  mutate(SamplingEventID = factor(SamplingEventID,unique(SamplingEventID))) %>% 
  pivot_longer(cols = str_replace(gene_order_long,"FLT3_ITD","FLT3.ITD"), names_to = "Gene",values_to="status") %>% 
  select(SamplingEventID,Gene,status,PapaemClass_short) %>% 
  filter(Gene %in% gene_order_short) %>% 
  ggplot(aes(SamplingEventID,factor(Gene,levels = rev(gene_order_short)))) + 
  geom_tile(aes(fill=status),colour = "black",size=0.3) + 
  facet_grid(~factor(PapaemClass_short, levels = grouplevels),scales="free_x",space="free") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "right", axis.ticks = element_blank(), axis.title.x = element_blank(), strip.text.x = element_text(angle=90), axis.text.x = element_blank(), axis.text.y = element_text(colour = "black")) + 
  scale_fill_manual(labels=c("CHET" = "Compund heterozygous", "HET" = "Heterozygous/Present", "HOM" = "Homozygous", "NA" = "No mutation/Fusion/CNA") ,values=c("HET" = "#2b83ba","CHET" = "#7b3294","HOM" = "#d7191c", "NA" = "gray80"),na.value="transparent") + 
  labs(y="Mutated genes", fill="Mutation type/Fusions/CNAs") + 
  geom_hline(yintercept=lines, color="white") + 
  scale_y_discrete(labels=c("FLT3.ITD" = "FLT3-ITD")) 

L1 <- get_legend(p1Label+ theme_gray(base_size = 6) + theme(legend.justification=c(0,1), legend.key.size = unit(0.1,"cm"))) 

#wrap_plots(L1)

#############################
# Fusion gene part of the plot

Fusionorder <- c("CBFB::MYH11","PML::RARA","KMT2A fusion","RUNX1::RUNX1T1","DEK::NUP214","GATA2::MECOM","other MECOM-r","KAT6A::NCOA2","KAT6A::EP300","HOXA9::NUP98","NUP214::SET")

darkcells2 <- Sample_data %>% 
  arrange(DNMT3A,NPM1,FLT3.ITD,TET2,RUNX1,SRSF2,TP53,ASXL1,IDH2,NRAS,IDH1,STAG2,SamplingEventID) %>% 
  mutate(temp=TRUE,SamplingEventID = factor(SamplingEventID,unique(SamplingEventID))) %>% 
  pivot_wider(names_from=FusionGene,values_from="temp") %>% select(!"NA") %>%
  pivot_longer(cols = c("PML::RARA":"CBFB::MYH11"), names_to = "FusionGene",values_to="Value") %>% 
  select(SamplingEventID,PapaemClass_short,FusionGene,Value) %>% 
  filter(FusionGene %in% c("CBFB::MYH11","PML::RARA","KMT2A fusion","RUNX1::RUNX1T1","DEK::NUP214") & is.na(Value)) %>% 
  mutate(Value = "NA")

lines <- c(4.5,9.5)
p2 <- Sample_data %>% 
  arrange(DNMT3A,NPM1,FLT3.ITD,TET2,RUNX1,SRSF2,TP53,ASXL1,IDH2,NRAS,IDH1,STAG2,SamplingEventID) %>% 
  mutate(temp=TRUE,
         SamplingEventID = factor(SamplingEventID,unique(SamplingEventID)), 
         FusionGene = str_replace_all(FusionGene, c("(ETV6::MECOM)|(LINC00824::MECOM)" = "other MECOM-r"))) %>% 
  pivot_wider(names_from=FusionGene,values_from="temp") %>% select(!"NA") %>%
  pivot_longer(cols = c("PML::RARA":"CBFB::MYH11"), names_to = "FusionGene",values_to="Value") %>% 
  select(SamplingEventID,PapaemClass_short,FusionGene,Value) %>% 
  ggplot(aes(SamplingEventID,factor(FusionGene,levels = rev(Fusionorder)))) + 
  geom_tile(aes(fill=Value),colour = "white") + 
  facet_grid(~factor(PapaemClass_short, levels = grouplevels),scales="free_x",space="free") + 
  theme_gray(base_size = 6) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none", 
        axis.ticks = element_blank(), 
        axis.title.x = element_blank(), 
        strip.text.x = element_text(angle=90), 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(colour = "black",face = "italic"),
        panel.spacing = unit(1, "pt"),
        plot.margin = margin(t = 0,r = 1.5,b=0, l=1.5,unit = "pt") ) +
  scale_fill_manual(values=c("TRUE" = "#2b83ba", "NA" = "gray80"),na.value="transparent") + 
  labs(y="Fusions") + 
  #geom_hline(yintercept=lines, color="white") + 
  geom_tile(aes(fill=Value), color="white", data = darkcells2)


p2grob <- ggplotGrob(p2)
panels <- grep("panel", p2grob$layout$name)
top <- unique(p2grob$layout$t[panels])
p2grob = p2grob[-(top-1), ]

#wrap_plots(p2grob)

#############################
# CNA part of the plot

darkcells3 <- Sample_data %>% 
  arrange(DNMT3A,NPM1,FLT3.ITD,TET2,RUNX1,SRSF2,TP53,ASXL1,IDH2,NRAS,IDH1,STAG2,SamplingEventID) %>% 
  mutate(SamplingEventID = factor(SamplingEventID,unique(SamplingEventID))) %>% 
  pivot_longer(cols = c(del5q:dup8q), names_to = "CNA",values_to="status") %>%
  select(SamplingEventID,PapaemClass_short,CNA,status) %>% 
  filter(CNA %in% c("del7q","dup8q","del5q","del17p","del12p") & is.na(status)) %>% 
  mutate(status = "NA")

lines <- c(5.5)
p3 <- Sample_data %>% 
  arrange(DNMT3A,NPM1,FLT3.ITD,TET2,RUNX1,SRSF2,TP53,ASXL1,IDH2,NRAS,IDH1,STAG2,SamplingEventID) %>% 
  mutate(SamplingEventID = factor(SamplingEventID,unique(SamplingEventID))) %>% 
  pivot_longer(cols = c(del5q:dup8q), names_to = "CNA",values_to="status") %>% 
  ggplot(aes(SamplingEventID,fct_rev(reorder(CNA,is.na(status),sum)))) + 
  geom_tile(aes(fill=status),colour = "white") + 
  facet_grid(~factor(PapaemClass_short, levels = grouplevels),scales="free_x",space="free") + 
  theme_gray(base_size = 6) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none", 
        axis.ticks = element_blank(), 
        axis.title.x = element_blank(), 
        strip.text.x = element_text(angle=90), 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(colour = "black"),
        panel.spacing = unit(1, "pt"),
        plot.margin = margin(t = 0,r = 1.5,b=0, l=1.5,unit = "pt") ) +
  scale_fill_manual(values=c("TRUE"= "#2b83ba", "NA" = "gray80"),na.value="transparent") + 
  labs(y="CNAs") + 
  #geom_hline(yintercept=lines, color="white") + 
  scale_y_discrete(labels=c("del7q" = "-7/7q", "dup8q" = "+8/8q", "del17p" = "-17/17p","del5q" = "-5/5q", "del9q" = "-9q", "del12p" = "-12/12p", "minusY" = "-Y", "del18q" = "-18/18q", "del20q" = "-20/20q", "dup11q" = "+11/11q")) +
  geom_tile(aes(fill=status),colour="white", data = darkcells3)



p3grob <- ggplotGrob(p3)
panels <- grep("panel", p3grob$layout$name)
top <- unique(p3grob$layout$t[panels])
p3grob = p3grob[-(top-1), ]

#wrap_plots(p3grob)


#############################
# Sample information part of the plot

p4 <- Sample_data %>% 
  arrange(DNMT3A,NPM1,FLT3.ITD,TET2,RUNX1,SRSF2,TP53,ASXL1,IDH2,NRAS,IDH1,STAG2,SamplingEventID) %>% 
  mutate(SamplingEventID = factor(SamplingEventID,unique(SamplingEventID))) %>% 
  pivot_longer(cols = c(AML.history,AML.timepoint,RNA_Tissue,FAB,ELN2022,Sex), names_to = "Clinical",values_to="status") %>% 
  mutate(Clinical = str_replace_all(Clinical,c("AML.timepoint" = "AML timepoint", "AML.history" = "AML history", "ELN2022" = "Risk (ELN2022)", "RNA_Tissue" = "Tissue")),
         status = case_when(Clinical == "AML timepoint" & status == "R" ~ "Relapse",
                            Clinical == "AML timepoint" & status == "D" ~ "Diagnosis",
                            Clinical == "AML history" & status == "APL" ~ "AML",
                            TRUE ~ status )) %>% 
  ggplot(aes(SamplingEventID,Clinical)) + 
  geom_tile(aes(fill=status),colour = "white") + 
  facet_grid(~factor(PapaemClass_short, levels = grouplevels),scales="free_x",space="free") + 
  theme_gray(base_size = 6) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        strip.text.x = element_text(angle=90),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_text(colour = "black"),
        panel.spacing = unit(1, "pt"),
        plot.margin = margin(t = 0,r = 1.5,b=1.5, l=1.5,unit = "pt") ) + 
  scale_fill_manual(values=c("Favorable" = "#31a354","Intermediate" = "#ff7f00", "Adverse" = "#d7191c", "No category" = "gray80", "Not Classified" = "gray80",
                             "Low-to-intermediate risk APL" = "#31a354", "High-risk APL" = "#d7191c",
                             "AML"="#2b83ba","APL"= "#4daf4a","BM"="#2b83ba", "Diagnosis"="#2b83ba","During disease" ="#ff7f00","F"= "#d7191c", "M" = "#2b83ba","M0"="#e41a1c","M1"="#377eb8","M2"="#4daf4a", "M3" = "#984ea3", "M4" = "#ff7f00", "M5" = "#ffff33", "M7"= "#a65628","PB" = "#d7191c","R" = "darkgrey","Relapse"="#d7191c","sAML"="#984ea3","tAML"= "#ff7f00"),na.value="gray80", 
                    labels=c("AML"="De novo AML/APL", "sAML" = "Secondary AML/APL", "tAML" = "Therapy-related AML","BM" = "Bone marrow","PB" = "Peripheral blood","F" = "Female", "M" = "Male")) + 
  #labs(x=Sampleinfo, y="Sample info")
  labs(y="Sample info")

p4grob <- ggplotGrob(p4)
panels <- grep("panel", p4grob$layout$name)
top <- unique(p4grob$layout$t[panels])
p4grob = p4grob[-(top-1), ]

#wrap_plots(p4grob)

#############################
# Labels for sample information part of the plot

p4Label1 <- Sample_data %>% 
  arrange(DNMT3A,NPM1,FLT3.ITD,TET2,RUNX1,SRSF2,TP53,ASXL1,IDH2,NRAS,IDH1,STAG2,SamplingEventID) %>% 
  mutate(SamplingEventID = factor(SamplingEventID,unique(SamplingEventID))) %>% 
  pivot_longer(cols = c(RNA_Tissue), names_to = "Clinical",values_to="status") %>% 
  mutate(Clinical = str_replace_all(Clinical,c("RNA_Tissue" = "Tissue"))) %>% 
  ggplot(aes(SamplingEventID,Clinical)) + 
  geom_tile(aes(fill=status),colour = "black",size=0.3) +  
  facet_grid(~factor(PapaemClass_short, levels = grouplevels),scales="free_x",space="free") + 
  theme(legend.position = "right",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks = element_blank(), strip.text.x = element_text(angle=90), axis.text.x = element_blank(), axis.text.y = element_text(colour = "black")) + 
  #  scale_fill_manual(values=c("#2b83ba","#fec44f","BM" = "#2b83ba","#2b83ba","#a1d99b","#a1d99b","#756bb1","#d7191c","#d7191c","#fec44f","#31a354","#2b83ba","#2b83ba","#ffffcc", "#ffeda0", "#fed976", "#feb24c", "#fd8d3c", "#fc4e2a", "#b10026","darkgrey","darkgrey","black","#d7191c","#bcbddc","#d7191c","#d7191c","#e5f5e0","#a1d99b","#31a354","#e31a1c","#deebf7","#deebf7","#9ecae1","#efedf5"),na.value="transparent") + 
  scale_fill_manual(values=c("BM" = "#2b83ba","PB" = "#d7191c"),na.value="transparent", 
                    labels = c("BM" = "Bone marrow","PB" = "Peripheral blood")) + 
  labs(x=Sampleinfo, fill="Tissue", y="Sample information")

p4Label2 <- Sample_data %>% 
  arrange(DNMT3A,NPM1,FLT3.ITD,TET2,RUNX1,SRSF2,TP53,ASXL1,IDH2,NRAS,IDH1,STAG2,SamplingEventID) %>% 
  mutate(SamplingEventID = factor(SamplingEventID,unique(SamplingEventID))) %>% 
  pivot_longer(cols = c(Sex), names_to = "Clinical",values_to="status") %>% 
  ggplot(aes(SamplingEventID,Clinical)) + 
  geom_tile(aes(fill=status),colour = "black",size=0.3) + 
  facet_grid(~factor(PapaemClass_short, levels = grouplevels),scales="free_x",space="free") + 
  theme(legend.position = "right",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks = element_blank(), strip.text.x = element_text(angle=90), axis.text.x = element_text(angle = 90, vjust = 0.1, hjust = 1.0, colour = "black"), axis.text.y = element_text(colour = "black")) + 
  scale_fill_manual(values=c("F"= "#d7191c", "M" = "#2b83ba" ),na.value="transparent", 
                    labels=c("F" = "Female", "M" = "Male")) + 
  labs(x=Sampleinfo, fill="Sex",y="Sample information")

p4Label3 <- Sample_data %>% 
  arrange(DNMT3A,NPM1,FLT3.ITD,TET2,RUNX1,SRSF2,TP53,ASXL1,IDH2,NRAS,IDH1,STAG2,SamplingEventID) %>% 
  mutate(SamplingEventID = factor(SamplingEventID,unique(SamplingEventID))) %>% 
  pivot_longer(cols = c(ELN2022), names_to = "Clinical",values_to="status") %>% 
  mutate(Clinical = str_replace_all(Clinical,c("ELN2022" = "Risk group (ELN2022)")),
         status = factor(status, levels = c("Favorable","Intermediate","Adverse","No category","Low-to-intermediate risk APL", "High-risk APL"))) %>% 
  ggplot(aes(SamplingEventID,Clinical)) + 
  geom_tile(aes(fill=status),colour = "black",size=0.3) + 
  facet_grid(~factor(PapaemClass_short, levels = grouplevels),scales="free_x",space="free") + 
  theme(legend.position = "right",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks = element_blank(), strip.text.x = element_text(angle=90), axis.text.x = element_blank(), axis.text.y = element_text(colour = "black")) + 
  scale_fill_manual(values=c("Favorable" = "#31a354","Intermediate" = "#ff7f00", "Adverse" = "#d7191c", "No category" = "gray80"),na.value="transparent",
                    labels=c("No category" = "Not Classified" )) +   
  labs(x=Sampleinfo, fill="Risk (ELN 2022)", y="Sample information")

p4Label4 <- Sample_data %>% 
  arrange(DNMT3A,NPM1,FLT3.ITD,TET2,RUNX1,SRSF2,TP53,ASXL1,IDH2,NRAS,IDH1,STAG2,SamplingEventID) %>% 
  mutate(SamplingEventID = factor(SamplingEventID,unique(SamplingEventID))) %>% 
  pivot_longer(cols = c(FAB), names_to = "Clinical",values_to="status") %>% 
  ggplot(aes(SamplingEventID,Clinical)) + 
  geom_tile(aes(fill=status),colour = "black",size=0.3) + 
  facet_grid(~factor(PapaemClass_short, levels = grouplevels),scales="free_x",space="free") + 
  theme(legend.position = "right",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks = element_blank(), strip.text.x = element_text(angle=90), axis.text.x = element_blank(), axis.text.y = element_text(colour = "black")) + 
  scale_fill_manual(values=c("M0"="#e41a1c","M1"="#377eb8","M2"="#4daf4a", "M3" = "#984ea3", "M4" = "#ff7f00", "M5" = "#ffff33", "M7"= "#a65628"),
                    na.value="transparent",
                    guide = guide_legend(ncol=2)) + 
  labs(x=Sampleinfo, y="Sample information", fill="FAB")

p4Label5 <- Sample_data %>% 
  arrange(DNMT3A,NPM1,FLT3.ITD,TET2,RUNX1,SRSF2,TP53,ASXL1,IDH2,NRAS,IDH1,STAG2,SamplingEventID) %>% 
  mutate(SamplingEventID = factor(SamplingEventID,unique(SamplingEventID))) %>% 
  pivot_longer(cols = c(AML.timepoint), names_to = "Clinical",values_to="status") %>% 
  mutate(Clinical = str_replace_all(Clinical,c("AML.timepoint" = "AML timepoint"))) %>%  
  ggplot(aes(SamplingEventID,Clinical)) + 
  geom_tile(aes(fill=status),colour = "black",size=0.3) + 
  facet_grid(~factor(PapaemClass_short, levels = grouplevels),scales="free_x",space="free") + 
  theme(legend.position = "right",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks = element_blank(), strip.text.x = element_text(angle=90), axis.text.x = element_blank(), axis.text.y = element_text(colour = "black")) + 
  scale_fill_manual(values=c("Diagnosis"="#2b83ba","During disease" ="#ff7f00","Relapse"="#d7191c"),na.value="transparent") + 
  labs(x=Sampleinfo, y="Sample information", fill="AML timepoint")

p4Label6 <- Sample_data %>% 
  arrange(DNMT3A,NPM1,FLT3.ITD,TET2,RUNX1,SRSF2,TP53,ASXL1,IDH2,NRAS,IDH1,STAG2,SamplingEventID) %>% 
  mutate(SamplingEventID = factor(SamplingEventID,unique(SamplingEventID))) %>% 
  pivot_longer(cols = c(AML.history), names_to = "Clinical",values_to="status") %>% 
  mutate(Clinical = str_replace_all(Clinical,c("AML.history" = "AML history"))) %>% 
  ggplot(aes(SamplingEventID,Clinical)) + 
  geom_tile(aes(fill=status),colour = "black",size=0.3) + 
  facet_grid(~factor(PapaemClass_short, levels = grouplevels),scales="free_x",space="free") + 
  theme(legend.position = "right",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks = element_blank(), strip.text.x = element_text(angle=90), axis.text.x = element_blank(), axis.text.y = element_text(colour = "black")) + 
  scale_fill_manual(values=c("AML"="#2b83ba","Relapse"="#d7191c","sAML"="#984ea3","tAML"= "#ff7f00"),na.value="transparent", 
                    labels=c("AML"="De novo AML/APL", "sAML" = "Secondary AML/APL", "tAML" = "Therapy-related AML")) + 
  labs(x=Sampleinfo, y="Sample information", fill="AML history")

L4 <- get_legend(p4Label1+theme_gray(base_size = 6) + theme(legend.justification=c(0,1), legend.key.size = unit(0.1,"cm"))) 
L5 <- get_legend(p4Label2+ theme_gray(base_size = 6) +theme(legend.justification=c(0,1), legend.key.size = unit(0.1,"cm"))) 
L6 <- get_legend(p4Label3+ theme_gray(base_size = 6) +theme(legend.justification=c(0,1), legend.key.size = unit(0.1,"cm"))) 
L7 <- get_legend(p4Label4+ theme_gray(base_size = 6) +theme(legend.justification=c(0,1), legend.key.size = unit(0.1,"cm"))) 
L8 <- get_legend(p4Label5+ theme_gray(base_size = 6) +theme(legend.justification=c(0,1), legend.key.size = unit(0.1,"cm"))) 
L9 <- get_legend(p4Label6+ theme_gray(base_size = 6) +theme(legend.justification=c(0,1), legend.key.size = unit(0.1,"cm")))

#wrap_plots(L4,L5,L6,L7,L8,L9,nrow = 1)


##############################################
# Stitch it all together:

Landplot <- arrangeGrob(rbind(p1grob,p2grob,p3grob,p4grob))

#wrap_plots(Landplot)

ggsave("Figures/Fig1A.pdf",Landplot, width = 170, height = 100,units = "mm")

Legends <- wrap_plots(A=L1,B=L4,C=L5,D=L6,E=L7,F=L8,G=L9,widths = c(1.3,0.75,0.4,0.75,0.5,0.5,1))

ggsave("Figures/Fig1A_Legend.pdf",Legends, height = 13, width = 154, units = "mm")


