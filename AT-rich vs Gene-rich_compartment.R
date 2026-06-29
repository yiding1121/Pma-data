HiC_eigs1_scores=read.csv('~Desktop//Pma_data/HiC/data/coolers-OPM/gene_rich_AT_rich_eigs_1_100k_scores.csv',sep = ",",header=T)

HiC_eigs1_scores

tmp <- data.frame(Eigs1_scores = c(HiC_eigs1_scores$Gene_rich, HiC_eigs1_scores$AT_rich),
                  group = c(rep("gene", length(HiC_eigs1_scores$Gene_rich)), rep("AT", length(HiC_eigs1_scores$AT_rich))))

tmp

p <- ggplot(tmp, aes(x=group, y=value, fill = group)) + 
  geom_violin(trim=TRUE) + stat_summary(fun = 'median',
                                 geom = "crossbar", 
                                 width = 0.5,
                                 colour = "black") + stat_compare_means(comparisons = test_sign, label = "p.signif")

p

library(ggstatsplot)
library(ggplot2)
ggbetweenstats(bf.message = TRUE,centrality.plotting = median,
  data  = tmp,
  x     = group,
  y     = Eigs1_scores,
  title = "AT-rich vs Gene-rich"
)

ggbetweenstats(
  bf.message = TRUE,
  centrality.plotting = median,
  data  = tmp,
  x     = group,
  y     = Eigs1_scores,
  title = "AT-rich vs Gene-rich"
) +
  theme_bw(base_size = 16) +   
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold")
  )


