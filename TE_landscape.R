library(reshape)
library(ggplot2)
library(viridis)
library(hrbrthemes)
library(tidyverse)
library(gridExtra)


#divergence_year <- read.table("/Users/yiding/Desktop/Pma_data/repeatmasker/bga_3_rounds_masked.align.landscape.My.Rclass.tab", sep="\t", 
#                              row.names = 1)


divergence_year <- read.table("/Users/yiding/Desktop/Pma_data/repeatmasker/bga_3_rounds_masked.align.landscape.Div.Rclass.tab", sep="\t", 
                              row.names = 1)


divergence_year <- t(divergence_year)
rownames(divergence_year) <- 1:50
kd_melt2 = melt(divergence_year[, 2:6], id="Rclass")

kd_melt2$X2 <- factor(kd_melt2$X2, levels = c("DNA", "LTR", "MITE", "LINE", "Unknown"))
ggplot(kd_melt2, aes(y=value, x=X1, fill = X2)) + 
  geom_bar(position="stack", stat="identity",color="white") +
  scale_fill_manual(values = c("#995d81","#FE5E41","#f6f740","#FFE88A","#9BC2DC")) +
  theme_classic() +
  xlab("Divergence(%)") +
  ylab("Length of sequence (bp)") + 
  labs(fill = "") +
  coord_cartesian(xlim = c(0, 55)) +
  theme(axis.text=element_text(size=11),axis.title =element_text(size=12)) 


data <- readxl::read_xlsx("/Users/yiding/Desktop/Pma_data/repeatmasker/TE_analysis_results/AT_vs_TAD_TE.xlsx", sheet = 2)

data <- readxl::read_xlsx("/Users/yiding/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)(2)/Pma_data/repeatmasker/TE_analysis_results/AT_vs_TAD_TE.xlsx", sheet = 2)

data <- data[, c(1,2,3,5)]
colnames(data) <- c("type", "subtype", "subsubtype", "pvalue")
data$logP <- -log10(data$pvalue)

ggplot(data, aes(y = logP, x = type, fill = subtype)) +
  geom_boxplot() +
  facet_wrap(~subtype, nrow = 1) +
  geom_hline(yintercept=1.30103, linetype="dashed", 
             color = "red", size=2) +
  scale_fill_manual(values = c("#995d81","#FE5E41","#f6f740","#FFE88A","#9BC2DC")) +
  ylab("log10(p-value)") +
  geom_jitter(alpha = 0.2)

data <- readxl::read_xlsx("/Users/yiding/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)(2)/Pma_data/repeatmasker/TE_analysis_results/Obs_vs_random.xlsx")

data <- data[, c(1,3,6,9)] 
colnames(data) <- c("trancript_type", "overlap_category", "percentage", "type")

data <- data %>% 
  filter(overlap_category != "transcript")

ggplot(data, aes(fill = overlap_category, y= type, x=percentage)) + 
  geom_bar(position="stack", stat="identity", color="black") +
  facet_wrap(~trancript_type, ncol = 1, strip.position = "left") + 
  scale_fill_manual(values = c("#995d81","#FE5E41","#f6f740","#FFE88A","#9BC2DC","#5C95FF","#DAD6D6","#71F79F","#BF211E")) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        strip.background = element_blank(), strip.placement = "outside")


data <- read.table("/Users/yiding/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)(2)/Pma_data/orthofinder_output/Pm_data/Syntnet_network/Bga_clusters_id_and_number.txt")

data <- data %>% 
  group_by(V2) %>% 
  summarise(n = n())

ggplot(data, aes(x = n, y = V2)) +
  geom_point() 






