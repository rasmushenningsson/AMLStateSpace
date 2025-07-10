library(tidyverse)
library(gridExtra)
library(cowplot)
library(patchwork)

##############################
#
# Prepare data:


Sample_data <- read.csv("data/sample_data.csv")

scAMLs <- c("AML126-AML-1","AML10-AML-1","AML151-AML-1","AML21-AML-1","AML110-AML-1",
            "AML62-AML-1","AML83-AML-1","AML48-AML-1","AML157-AML-1","AML32-AML-1",
            "AML138-AML-1","AML172-AML-1","AML85-Relapse-1", "AML85-AML-1", "AML111-AML-1", 
            "AML55-AML-1","AML80-AML-1","AML161-AML-1","AML61-AML-1","AML117-AML-1",
            "AML123-AML-1","AML23-AML-1","AML124-AML-1","AML136-AML-1","AML34-AML-1",
            "AML104-AML-1","AML79-AML-1","AML27-AML-1","AML25-AML-1","AML24-AML-1",
            "AML28-AML-1","AML97-AML-1","AML9-AML-1","AML7-AML-1","AML33-Relapse-1",
            "AML105-AML-1","AML37-AML-1","AML155-AML-1")

caseorder <- scAMLs

ncaseorder <- caseorder %>% gsub("-AML-1","",.) %>% gsub("AML85-Relapse-1","AML85-R",.) %>% gsub("-Relapse-1","",.) %>% gsub("AML85$","AML85-D",.)

n_cells <- read_csv("data/cell_proportions.csv")

n_cells_small <- n_cells %>% subset(str_detect(SamplingEventID,"^AML"))


##########################################
#
# Make the cell proportion barplot:

barplot <- n_cells_small %>% filter(celltype.aml == "AML Immature") %>% 
  ggplot(aes(y = prop, x=factor(SamplingEventID, levels = caseorder), fill = celltype.aml)) + 
  geom_bar(width = 1, stat = "identity", color = "white", size = 0.1) + 
  scale_y_continuous(labels=scales::percent,expand = expansion(mult = c(0.04, 0)),limits = c(0,1)) + 
  labs(y="AML Immature contents", x=element_blank()) + 
  scale_fill_manual(name = "Celltype", values=c("AML Immature" = "#fec44f")) + 
  theme_gray(base_size = 5) +
  theme(strip.text.x = element_text(angle = 90),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.margin = margin(t = 4,r = 2,b=0, l=2,unit = "pt") ) +
  scale_x_discrete(expand = c(0,0)) +
  geom_vline(xintercept= seq_along(unique(ncaseorder)) + 0.5, color = "white", size=.1)
#  geom_vline(xintercept= c(0.5,17.5,23.5,36.5,38.5), color = "black") + 
#  geom_vline(xintercept= c(5.5,11.5,20.5,25.5,31.5), color = "black", linetype = "longdash")


barplotgrob <- ggplotGrob(barplot)

#wrap_plots(barplotgrob)


#############################################
#
# Plot mutation information


selectedGeneOrder <- c("NPM1","DNMT3A","FLT3.ITD", "IDH2", "IDH1",
                       "RUNX1","SRSF2","ASXL1","BCOR", "TET2",
                       "TP53", "FLT3","NRAS","KRAS","CEBPA")

darkcells2 <-  Sample_data %>% 
  arrange(DNMT3A,NPM1,FLT3.ITD,TET2,RUNX1,SRSF2,TP53,ASXL1,IDH2,NRAS,IDH1,STAG2,SamplingEventID) %>% 
  mutate(SamplingEventID = factor(SamplingEventID,unique(SamplingEventID))) %>% 
  filter(SamplingEventID %in% scAMLs) %>% 
  pivot_longer(cols = str_replace(selectedGeneOrder,"FLT3_ITD","FLT3.ITD"), names_to = "Gene",values_to="status") %>% 
  select(SamplingEventID,Gene,status,PapaemClass_short) %>% 
  filter(Gene %in% c("RUNX1","SRSF2","ASXL1","BCOR","TET2") & is.na(status)) %>% 
  mutate(status = "NA")

lines <- seq(from=4.5,to=19.5, by=5)

p2 <- 
  Sample_data %>% arrange(DNMT3A,NPM1,FLT3.ITD,TET2,RUNX1,SRSF2,TP53,ASXL1,IDH2,NRAS,IDH1,STAG2,SamplingEventID) %>% 
  mutate(SamplingEventID = factor(SamplingEventID,unique(SamplingEventID))) %>% 
  filter(SamplingEventID %in% scAMLs) %>% 
  pivot_longer(cols = str_replace(selectedGeneOrder,"FLT3_ITD","FLT3.ITD"), names_to = "Gene",values_to="status") %>% 
  droplevels() %>%  
  ggplot(aes(x=factor(SamplingEventID, levels = caseorder),factor(Gene,levels = rev(selectedGeneOrder)))) + 
  geom_tile(aes(fill=status),colour = "white") + 
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
  scale_x_discrete(expand = c(0,0)) +
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
  ggplot(aes(x = factor(SamplingEventID, levels = caseorder),fct_rev(reorder(FusionGene,is.na(Value),sum)))) +
  geom_tile(aes(fill=Value),colour = "white") +
  theme_gray(base_size = 5) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    legend.position = "none", 
    axis.ticks = element_blank(), 
    axis.title.x = element_blank(), 
    strip.text.x = element_blank(),
    #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1.0, colour = "black"), 
    axis.text.x = element_blank(),
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

#############################################
#
# Plot subtype information

p4 <- Sample_data %>% 
  arrange(DNMT3A,NPM1,FLT3.ITD,TET2,RUNX1,SRSF2,TP53,ASXL1,IDH2,NRAS,IDH1,STAG2,SamplingEventID) %>% 
  mutate(SamplingEventID = factor(SamplingEventID,unique(SamplingEventID))) %>% 
  pivot_longer(cols = c(PapaemClass_short), names_to = "Clinical",values_to="status") %>% 
  mutate(Clinical = str_replace_all(Clinical,"PapaemClass_short","AML subtype")) %>% 
  filter(SamplingEventID %in% scAMLs) %>% 
  ggplot(aes(x=factor(SamplingEventID, levels = caseorder),Clinical)) + 
  geom_tile(aes(fill=status),colour = "white") + 
  theme_gray(base_size = 5) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        strip.text.x = element_text(angle=90),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1.0, colour = "black"), 
        axis.title.y = element_blank(),
        axis.text.y = element_text(colour = "black"),
        panel.spacing = unit(1, "pt"),
        plot.margin = margin(t = 0,r = 2,b=2, l=2,unit = "pt") ) + 
  scale_fill_manual(values=c("NPM1"="#4e9eee","AML-MR"="#36e548",
                             "TP53"="#e31a1c","In two subgroups" = "#000000",
                             "No class-defining"= "#b2b2b2","CBFB::MYH11" = "#c261eb", 
                             "RUNX1::RUNX1T1" = "#ffaa00"    ),
                    na.value="transparent") + 
  labs(y="") +
  scale_x_discrete(expand = c(0,0),
                   labels=function(x) {
                     lbl <- str_replace_all(x,c("-AML-1"="","-Relapse-1" = "R"))
                     lbl}
  ) + 
  scale_y_discrete(expand = c(0,0))

p4grob <- ggplotGrob(p4)
panels <- grep("panel", p4grob$layout$name)
top <- unique(p4grob$layout$t[panels])
p4grob = p4grob[-(top-1), ]


#wrap_plots(p4grob)





#################################################
# Stitch it all together:


full_plot <- plot_grid(barplotgrob, NULL,p2grob, NULL,p3grob,NULL, p4grob, 
                  ncol=1, align="v", axis = "lr",
                  rel_heights = c(1.5,-0.066,1.7,-0.024,0.6,-0.063,0.5))

ggsave("Figures/Fig4A.pdf", full_plot, width = 107.1, height = 75,units = "mm")

