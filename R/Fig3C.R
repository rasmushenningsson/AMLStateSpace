library(tidyverse)
library(gridExtra)
library(cowplot)
library(patchwork)

##############################
#
# Prepare data:


Sample_data <- read.csv("data/sample_data.csv")

Sample_data <- Sample_data %>% 
  add_row(SamplingEventID = c("NBM4-MNC", "NBM7-MNC", "NBM8-MNC", "NBM10-MNC", "NBM11-MNC", "NBM8-CD34", "NBM10-CD34", "NBM11-CD34"), 
          PapaemClass = rep(c("NBM-MNC","NBM-CD34"), each=5, len=8))

grouplevels <- c("NPM1","AML-MR","TP53","In two subgroups","No class-defining","CBFB::MYH11","RUNX1::RUNX1T1", "NBM-MNC","NBM-CD34")

scAMLs <- c("AML79-AML-1","AML27-AML-1","AML28-AML-1","AML34-AML-1","AML32-AML-1","AML21-AML-1","AML23-AML-1","AML25-AML-1",
            "AML110-AML-1","AML123-AML-1","AML138-AML-1","AML7-AML-1","AML9-AML-1","AML24-AML-1","AML33-Relapse-1","AML97-AML-1",
            "AML104-AML-1","AML105-AML-1","AML136-AML-1","AML10-AML-1","AML55-AML-1","AML111-AML-1","AML126-AML-1","AML151-AML-1",
            "AML157-AML-1","AML172-AML-1","AML37-AML-1","AML48-AML-1","AML62-AML-1","AML80-AML-1","AML83-AML-1","AML85-AML-1",
            "AML155-AML-1","AML61-AML-1","AML85-Relapse-1","AML117-AML-1","AML124-AML-1","AML161-AML-1", "NBM4-MNC", "NBM7-MNC", 
            "NBM8-MNC", "NBM10-MNC", "NBM11-MNC", "NBM8-CD34", "NBM10-CD34", "NBM11-CD34")

caseorder <- 
  Sample_data %>% arrange(DNMT3A,NPM1,FLT3.ITD,TET2,RUNX1,SRSF2,TP53,ASXL1,IDH2,NRAS,IDH1,STAG2,SamplingEventID) %>% 
  mutate(SamplingEventID = factor(SamplingEventID,unique(SamplingEventID))) %>%
  filter(SamplingEventID %in% scAMLs) %>% 
  droplevels() %>% pull(SamplingEventID) %>% levels()

n_cells <- read_csv("data/cell_proportions.csv")

n_cells$PapaemClass <- factor(n_cells$PapaemClass, levels = grouplevels)
cellorder <- c("AML Immature","HSC","LMPP","GMP","Monocytes","Megakaryocytic cells","Erythroid cells","Dendritic cells","NK-cells","T-cells", "B-cells")
n_cells$celltype.aml <- factor(n_cells$celltype.aml, levels = cellorder)


##########################################
#
# Make the cell proportion barplot:

cellcolors <- c("AML Immature" = "#fec44f",
                "HSC"="#66b266", 
                "LMPP"="#756bb1",
                "GMP"="#008000",
                "Monocytes"="#fa9fb5",
                "Megakaryocytic cells"="#e5687e", 
                "Erythroid cells"= "#a7241d", 
                "Dendritic cells"= "#ff4d00", 
                "NK-cells"="#196fbe",
                "T-cells"="#7fb4e5", 
                "B-cells"= "#a64ca6" )

barplot <- 
  ggplot(n_cells, aes(y = prop, x=factor(SamplingEventID, levels = caseorder), fill = celltype.aml)) + 
  geom_bar(width = 1, stat = "identity", color = "white", size = 0.1) +
  scale_y_continuous(labels=scales::percent,expand = expansion(mult = c(0.04, 0))) +
  scale_x_discrete(expand = c(0,0)) +
  labs(y="Cell type", x=element_blank()) + 
  facet_grid( ~PapaemClass,scales="free", space="free_x") + 
  scale_fill_manual(name = "Celltype", values=cellcolors) +
  theme_gray(base_size = 5) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none", 
        axis.ticks = element_blank(), 
        axis.title.x = element_blank(), 
        strip.text.x = element_text(angle=90), 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(colour = "black"),
        panel.spacing = unit(1, "pt"),
        plot.margin = margin(t = 2,r = 2,b=0, l=2,unit = "pt") )

barplotgrob <- ggplotGrob(barplot)

#wrap_plots(barplotgrob)


#############################################
#
# Plot mutation information


gene_order_bar <- c("NPM1","DNMT3A","FLT3.ITD", "IDH2", "IDH1",
                         "RUNX1","SRSF2","ASXL1","BCOR","STAG2",
                         "TP53", "FLT3","TET2","NRAS","APOB",
                         "CEBPA","KRAS","BCORL1","EZH2","JAK2",
                         "PTPN11", "SETBP1","WT1","ZRSR2"
)

darkcells2 <-  Sample_data %>% 
  arrange(DNMT3A,NPM1,FLT3.ITD,TET2,RUNX1,SRSF2,TP53,ASXL1,IDH2,NRAS,IDH1,STAG2,SamplingEventID) %>% 
  mutate(SamplingEventID = factor(SamplingEventID,unique(SamplingEventID))) %>% 
  filter(SamplingEventID %in% scAMLs) %>% 
  pivot_longer(cols = str_replace(gene_order_bar,"FLT3_ITD","FLT3.ITD"), names_to = "Gene",values_to="status") %>% 
  select(SamplingEventID,Gene,status,PapaemClass_short) %>% 
  filter(Gene %in% c("RUNX1","SRSF2","ASXL1","BCOR","STAG2","CEBPA","KRAS","BCORL1","EZH2","JAK2") & is.na(status)) %>% 
  mutate(status = "NA")

lines <- seq(from=4.5,to=19.5, by=5)

p2 <- 
  Sample_data %>% arrange(DNMT3A,NPM1,FLT3.ITD,TET2,RUNX1,SRSF2,TP53,ASXL1,IDH2,NRAS,IDH1,STAG2,SamplingEventID) %>% 
  mutate(SamplingEventID = factor(SamplingEventID,unique(SamplingEventID))) %>% 
  filter(SamplingEventID %in% scAMLs) %>% 
  pivot_longer(cols = str_replace(gene_order_bar,"FLT3_ITD","FLT3.ITD"), names_to = "Gene",values_to="status") %>% 
  droplevels() %>%  
  ggplot(aes(SamplingEventID,factor(Gene,levels = rev(gene_order_bar)))) + 
  geom_tile(aes(fill=status),colour = "white") + 
  facet_grid(~factor(PapaemClass_short, levels = grouplevels),scales="free_x",space="free") +
  theme_gray(base_size = 5) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none", 
        axis.ticks = element_blank(), 
        axis.title.x = element_blank(), 
        strip.text.x = element_text(angle=90), 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(colour = "black", face = "italic"),
        panel.spacing = unit(1, "pt"),
        plot.margin = margin(t = 0,r = 2,b=0, l=2,unit = "pt") ) +
  scale_fill_manual(values=c("CHET" = "#7b3294","HET" = "#2b83ba","HOM" = "#d7191c", "NA" = "gray80"),na.value="transparent") +
  labs(y="Mutated genes") + 
  #geom_hline(yintercept=lines, color="white") +
  scale_y_discrete(labels=c("FLT3.ITD" = "FLT3-ITD"),expand = c(0,0)) +
  scale_x_discrete(expand=c(0,0)) +
  geom_tile(aes(fill=status), color="white", data = darkcells2)

p2grob <- ggplotGrob(p2)
panels <- grep("panel", p2grob$layout$name)
top <- unique(p2grob$layout$t[panels])
p2grob = p2grob[-(top-1), ]


#wrap_plots(p2grob)

#############################################
#
# Plot gene fusion information

p3 <- 
  Sample_data %>% arrange(DNMT3A,NPM1,FLT3.ITD,TET2,RUNX1,SRSF2,TP53,ASXL1,IDH2,NRAS,IDH1,STAG2,SamplingEventID) %>% 
  mutate(SamplingEventID = factor(SamplingEventID,unique(SamplingEventID)), temp=TRUE) %>%
  filter(SamplingEventID %in% scAMLs) %>% 
  pivot_wider(names_from=FusionGene,values_from="temp") %>% 
  select(!"NA") %>% 
  pivot_longer(cols = c("RUNX1::RUNX1T1":"CBFB::MYH11"), names_to = "FusionGene",values_to="Value")%>%
  ggplot(aes(SamplingEventID,fct_rev(reorder(FusionGene,is.na(Value),sum)))) +
  geom_tile(aes(fill=Value),colour = "white") +
  facet_grid(~factor(PapaemClass_short, levels = grouplevels),scales="free_x",space="free") + 
  theme_gray(base_size = 5) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    legend.position = "none", 
    axis.ticks = element_blank(), 
    axis.title.x = element_blank(), 
    strip.text.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1.0, colour = "black"), 
    axis.text.y = element_text(colour = "black", face="italic"),
    panel.spacing = unit(1, "pt"),
    plot.margin = margin(t = 0,r = 2,b=2, l=2,unit = "pt"),
    panel.background = element_rect(fill = "gray80")) +
  scale_fill_manual(values=c("#2b83ba"),na.value="transparent") + 
  labs(y="Fusions") + 
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(labels=function(x) {
    lbl <- str_replace_all(x,c("-AML-1"="","-Relapse-1" = "R"))
    lbl},
    expand=c(0,0)
  )



p3grob <- ggplotGrob(p3)
panels <- grep("panel", p3grob$layout$name)
top <- unique(p3grob$layout$t[panels])
p3grob = p3grob[-(top-1), ]


#wrap_plots(p3grob)


#################################################
# Stitch it all together:

full_plot <- plot_grid(barplot, NULL,p2grob, NULL,p3grob, ncol=1, align="v", axis = "lr", rel_heights = c(2,-0.022,1.8,-0.024,0.95))

ggsave("Figures/Fig3C.pdf", full_plot, width = 78.5, height = 85,units = "mm")

#################################################
# Make Legend

bplotLegend <- 
  ggplot(n_cells, aes(y = prop, x=factor(SamplingEventID, levels = caseorder), fill = celltype.aml)) + 
  geom_bar(width = 1, stat = "identity", color = "black", size = 0.3) +
  scale_y_continuous(labels=scales::percent) + 
  labs(y="Cell type", x=element_blank()) + 
  facet_grid( ~PapaemClass,scales="free", space="free_x") + 
  scale_fill_manual(name = "Celltype", values=cellcolors) +
  theme_gray(base_size = 5) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "right", 
        axis.ticks = element_blank(), 
        axis.title.x = element_blank(), 
        strip.text.x = element_text(angle=90), 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(colour = "black"),
        panel.spacing = unit(1, "pt"),
        plot.margin = margin(t = 2,r = 2,b=0, l=2,unit = "pt") )


Leg1 <- get_legend(bplotLegend+ theme_gray(base_size = 5) + theme(legend.justification=c(0,1), legend.key.size = unit(0.2,"cm"))) 

p2legend <- 
  Sample_data %>% arrange(DNMT3A,NPM1,FLT3.ITD,TET2,RUNX1,SRSF2,TP53,ASXL1,IDH2,NRAS,IDH1,STAG2,SamplingEventID) %>% 
  mutate(SamplingEventID = factor(SamplingEventID,unique(SamplingEventID))) %>% 
  filter(SamplingEventID %in% scAMLs) %>% 
  pivot_longer(cols = str_replace(gene_order_bar,"FLT3_ITD","FLT3.ITD"), names_to = "Gene",values_to="status") %>% 
  droplevels() %>%  
  ggplot(aes(x=factor(SamplingEventID, levels = caseorder),factor(Gene,levels = rev(gene_order_bar)))) + 
  geom_tile(aes(fill=status),colour = "black",size=0.3) + 
  theme_gray(base_size = 5) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "right", 
        axis.ticks = element_blank(), 
        axis.title.x = element_blank(), 
        strip.text.x = element_text(angle=90), 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(colour = "black", face = "italic"),
        panel.spacing = unit(1, "pt"),
        plot.margin = margin(t = 0,r = 2,b=0, l=2,unit = "pt") ) +
  scale_fill_manual(labels=c("CHET" = "Compund heterozygous", "HET" = "Heterozygous/Present", "HOM" = "Homozygous", "NA" = "No mutation/fusion"),
                    values=c("HET" = "#2b83ba","CHET" = "#7b3294","HOM" = "#d7191c", "NA" = "gray80"),na.value="transparent") +
  labs(y="Mutated genes",fill="Mutation type/fusion") + 
  #geom_hline(yintercept=lines, color="white") +
  scale_y_discrete(labels=c("FLT3.ITD" = "FLT3-ITD"),expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) 


Leg2 <- get_legend(p2legend+ theme_gray(base_size = 5) + theme(legend.justification=c(0,1), legend.key.size = unit(0.2,"cm"))) 

#wrap_plots(Leg2)



Full_legend <- wrap_plots(Leg1,Leg2,ncol=1,heights = c(1,0.5))

ggsave("Figures/Fig3C_Legend.pdf",Full_legend,height = 45, width = 25, units="mm")




