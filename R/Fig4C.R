library(tidyverse)
library(survminer)
require("survival")
library(patchwork)


Lund_survival <- read_csv("data/NPM1_survival_data_Lund.csv")
TCGA_survival <- read_csv("data/NPM1_survival_data_TCGA.csv")
Beat2_survival <- read_csv("data/NPM1_survival_data_Beat-AML2.csv")
Clinseq_survival <- read_csv("data/NPM1_survival_data_Clinseq.csv")


merged_survival <- rbind(Lund_survival,TCGA_survival,Beat2_survival,Clinseq_survival)

merged_survival_intensive_NPM1class <- merged_survival %>% filter(Treatment == "Chemo" & NPM1group %in% c("NPM1 class I", "NPM1 class II")) %>% 
  mutate(NPM1group = factor(NPM1group, levels = c("NPM1 class I","NPM1 class II")))

os_merged <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ NPM1group + SCT.group, data = merged_survival_intensive_NPM1class)

os_merged_plot <- ggsurvplot(os_merged,
                              pval = TRUE,
                              legend = "none",
                              legend.labs = c("NPM1 class I/No SCT","NPM1 class I/SCT","NPM1 class II/No SCT","NPM1 class II/SCT"),
                              palette = c("lightblue","#377eb8", "lightcoral", "#e41a1c"),
                              risk.table = TRUE,
                              xlim = c(0,100),
                              #censor.shape = 124,
                              size = 0.6,
                              xlab = "Time (years)",
                              ylab = "Probability of Survival",
                              xscale = "m_y",
                              break.time.by = 24,
                              break.y.by = 0.2
)



KMplot <- os_merged_plot$plot + 
  annotate("text", x=rep(62,4), y = c(0.65,0.34,0.23,0.11), label = c("NPM1 class I/SCT","NPM1 class I/No SCT","NPM1 class II/No SCT","NPM1 class II/SCT"), hjust=0,
           color = c("#377eb8","lightblue","lightcoral","#e41a1c"))

table <- os_merged_plot$table + theme(axis.text.y.left = element_text(face = "bold"))

arranged <- KMplot / table + plot_layout(heights = c(0.85,0.15))

ggsave("Figures/Fig4C.pdf", arranged, width = 8, height = 8, device = cairo_pdf)


