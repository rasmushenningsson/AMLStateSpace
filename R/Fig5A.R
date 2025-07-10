library(tidyverse)
library(patchwork)
library(RColorBrewer)
library(scales)

NPM1classI <- c("AML24","AML25","AML27","AML28","AML79","AML104","AML136")
NPM1classII <- c("AML7","AML9","AML33","AML97","AML105")

NBMMNC <- c("NBM4-MNC","NBM7-MNC","NBM8-MNC","NBM10-MNC","NBM11-MNC") 
NBMCD34 <- c("NBM8-CD34","NBM10-CD34","NBM11-CD34") 
Other <- c("AML34", "AML155") 

CBFBMYH11 <- c("AML23","AML123","AML124") 
RUNX1RUNX1T1 <- c("AML61","AML117","AML161") 
chrsplice <- c("AML10","AML21","AML32","AML55","AML110","AML111","AML126","AML138","AML151","AML157","AML172")
TP53 <- c("AML37","AML48","AML62","AML80","AML83","AML85-D","AML85-R") 

AverageExp <- read_csv("data/Average_expression_SampleCelltype.csv.gz") %>% column_to_rownames("...1")
MetaData <- colnames(AverageExp) %>% as_tibble() %>% 
  separate(value,into = c("Sample","Celltype"),sep = "_",remove = F) %>% 
  mutate(Sample = str_remove_all(Sample,"SCT\\."),
         Sample = str_replace_all(Sample,c("\\."="-")),
         Celltype = str_replace_all(Celltype,c("\\."=" ")),
         AML.type = case_when(Sample %in% NBMMNC ~ "NBM-MNC",
                              Sample %in% NBMCD34 ~ "NBM-CD34",
                              Sample %in% NPM1classI ~ "NPM1 class I",
                              Sample %in% NPM1classII ~ "NPM1 class II",
                              Sample %in% Other ~ "Other",
                              Sample %in% TP53 ~ "TP53",
                              Sample %in% CBFBMYH11 ~ "CBFB::MYH11",
                              Sample %in% RUNX1RUNX1T1 ~ "RUNX1::RUNX1T1",
                              Sample %in% chrsplice ~ "AML-MR"),
         AML.type = factor(AML.type, levels = c("NBM-CD34","NBM-MNC","NPM1 class I","NPM1 class II","AML-MR","TP53","Other","CBFB::MYH11","RUNX1::RUNX1T1")),
         Celltype = factor(Celltype, levels = c("AML Immature","HSC","LMPP","GMP","Monocytes", "Megakaryocytic cells", "Erythroid", "Dendritic cells", "NK cells", "T cells", "B cells"))) %>% 
  rename("name" = "value")

LogNormAverage <- t(scale(t(log(AverageExp+0.01))))
ExprDataSC_LogNorm <- as.data.frame(LogNormAverage) %>% rownames_to_column("Gene") %>% pivot_longer(SCT.NBM4.MNC_HSC:SCT.AML124_B.cells) %>% left_join(MetaData, by = "name") %>% rename("AverageSCT" = "value")

MHCclass1 <- c("HLA-A","HLA-B","HLA-C")
MHCclass2 <- scan("data/HLA_MHC_Class2.txt", what="character")
MHCnonclass <- c("HLA-E","HLA-F","HLA-G")

AMLs <- c(NPM1classI, NPM1classII, Other, CBFBMYH11, RUNX1RUNX1T1, chrsplice, TP53)
NBMs <- c(NBMCD34)

Anno <- tibble(Gene = c(MHCclass1,MHCclass2,"CIITA", MHCnonclass), List = factor(c(rep("MHC-I", length(MHCclass1)),rep("MHC-II", length(MHCclass2)), "CIITA", rep("MHC-nonclassical", length(MHCnonclass)) ), levels = c("MHC-nonclassical","MHC-I","MHC-II","CIITA")))

fixlabels <- function(string) {str_replace_all(string, c("CBFB::MYH11" = "CBFB:: MYH11", "RUNX1::RUNX1T1" = "RUNX1:: RUNX1T1", "MHC-nonclassical" = "HLA non-classical"))}

MHCplot <- ExprDataSC_LogNorm %>% 
  filter(Gene %in% c(MHCclass1,MHCclass2,MHCnonclass,"CIITA") & ( Celltype == "AML Immature" | ( Sample %in% NBMCD34 & Celltype == "HSC")  ) ) %>% left_join(Anno, by = "Gene")  %>% 
  ggplot(aes(Sample,Gene,fill=AverageSCT)) + 
  geom_tile() + 
  scale_fill_gradient2(low = "#2b83ba", high = "#d7191c", mid="white",limits = c(-2,2), oob = scales::squish) +
  facet_grid(List~AML.type, scales = "free", space = "free", labeller = as_labeller(fixlabels, default = label_wrap_gen(width=10) ) ) +
  theme_classic(base_size = 5) +
  theme(axis.text.x = element_text(angle=90,hjust = 1), strip.text.y.right = element_text(angle=0), axis.title = element_blank()) 

ggsave("Figures/Fig5A.pdf",MHCplot,width=170,height = 50, units = "mm")
