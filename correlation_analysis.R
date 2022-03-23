library(tidyverse)
library(ggplot2)
library(reshape2)
library(dplyr)
library(egg)
library(ggpubr)


#read TCGA RSEM genes normalized into DATAFRAME

count.matrix.genes.rsem.norm <- read.csv(file = "<Path to RSEM Normalized Level 3 Gene Counts file>",
                                         header = T,
                                         sep = "\t",
                                         stringsAsFactors = F)

gene.names <- sapply(strsplit(count.matrix.genes.rsem.norm$Hybridization.REF, "\\|"), `[`, 1)
count.matrix.genes.rsem.norm$Hybridization.REF <- gene.names

#ATF4 gene targets (produced by microarray datasets)

filtered.genes.rsem.norm <- count.matrix.genes.rsem.norm[(count.matrix.genes.rsem.norm$Hybridization.REF == "COL1A2" |
                                                            count.matrix.genes.rsem.norm$Hybridization.REF == "TEKT2" |
                                                            count.matrix.genes.rsem.norm$Hybridization.REF == "THAP8" |
                                                            count.matrix.genes.rsem.norm$Hybridization.REF == "TLL1" |
                                                            count.matrix.genes.rsem.norm$Hybridization.REF == "TMCC1" |
                                                            count.matrix.genes.rsem.norm$Hybridization.REF == "GLMN" |
                                                            count.matrix.genes.rsem.norm$Hybridization.REF == "CDON" |
                                                            count.matrix.genes.rsem.norm$Hybridization.REF == "HDAC6" |
                                                            count.matrix.genes.rsem.norm$Hybridization.REF == "HTR3A" |
                                                            count.matrix.genes.rsem.norm$Hybridization.REF == "KBTBD6" |
                                                            count.matrix.genes.rsem.norm$Hybridization.REF == "KIR3DL3" |
                                                            count.matrix.genes.rsem.norm$Hybridization.REF == "ERLEC1" |
                                                            count.matrix.genes.rsem.norm$Hybridization.REF == "IQCB1" |
                                                            count.matrix.genes.rsem.norm$Hybridization.REF == "LYRM2" |
                                                            count.matrix.genes.rsem.norm$Hybridization.REF == "MTBP" |
                                                            count.matrix.genes.rsem.norm$Hybridization.REF == "NUF2" |
                                                            count.matrix.genes.rsem.norm$Hybridization.REF == "PMS2L3" |
                                                            count.matrix.genes.rsem.norm$Hybridization.REF == "RGS11" |
                                                            count.matrix.genes.rsem.norm$Hybridization.REF == "SCTR" |
                                                            count.matrix.genes.rsem.norm$Hybridization.REF == "KIAA1644" |
                                                            count.matrix.genes.rsem.norm$Hybridization.REF == "SNORD89" |
                                                            count.matrix.genes.rsem.norm$Hybridization.REF == "STX11" |
                                                            count.matrix.genes.rsem.norm$Hybridization.REF == "ZNHIT2" |
                                                            count.matrix.genes.rsem.norm$Hybridization.REF == "ABHD14B" |
                                                            count.matrix.genes.rsem.norm$Hybridization.REF == "LUZP6" |
                                                            count.matrix.genes.rsem.norm$Hybridization.REF == "ECSIT" |
                                                            count.matrix.genes.rsem.norm$Hybridization.REF == "LOC645166" |
                                                            count.matrix.genes.rsem.norm$Hybridization.REF == "GLP1R" |
                                                            count.matrix.genes.rsem.norm$Hybridization.REF == "GOLPH3" |
                                                            count.matrix.genes.rsem.norm$Hybridization.REF == "HLA-G" |
                                                            count.matrix.genes.rsem.norm$Hybridization.REF == "LHFPL5" |
                                                            count.matrix.genes.rsem.norm$Hybridization.REF == "MSL3" |
                                                            count.matrix.genes.rsem.norm$Hybridization.REF == "RHEBL1"), ]


#if melanoma samples, use primary + metastatic samples - else use only primary samples
#change gene name below (i.e., COL1A2) to check for other genes

primary.metastatic.melanomas <- filtered.genes.rsem.norm %>% dplyr::select(contains("Hybridization.REF") | contains(".01A.") | contains(".06A."))
row.names(primary.metastatic.melanomas) <- primary.metastatic.melanomas$Hybridization.REF
primary.metastatic.melanomas <- primary.metastatic.melanomas[-1]
primary.metastatic.melanomas <- t(primary.metastatic.melanomas)
primary.metastatic.melanomas <- apply(primary.metastatic.melanomas, 2, as.numeric)
primary.metastatic.melanomas <- primary.metastatic.melanomas + 0.1
primary.metastatic.melanomas <- as.data.frame(log2(primary.metastatic.melanomas))
primary.metastatic.melanomas.four <- data.frame(primary.metastatic.melanomas[, 
                    names(primary.metastatic.melanomas) %in% c("COL1A2")])
primary.metastatic.melanomas.four$mean_exp <- rowMeans(primary.metastatic.melanomas.four)
primary.metastatic.melanomas.targets <- data.frame(primary.metastatic.melanomas[, 
                    !names(primary.metastatic.melanomas) %in% c("COL1A2")])
primary.metastatic.melanomas.targets$mean_exp <- rowMeans(primary.metastatic.melanomas.targets)
primary.metastatic.melanomas.targets$mean_exp_four <- primary.metastatic.melanomas.four$mean_exp
primary.metastatic.melanomas.targets$COL1A2 <- primary.metastatic.melanomas$COL1A2


primary.melanomas <- filtered.genes.rsem.norm %>% dplyr::select(contains("Hybridization.REF") | contains(".01A."))
row.names(primary.melanomas) <- primary.melanomas$Hybridization.REF
primary.melanomas <- primary.melanomas[-1]
primary.melanomas <- t(primary.melanomas)
primary.melanomas <- apply(primary.melanomas, 2, as.numeric)
primary.melanomas <- primary.melanomas + 0.1
primary.melanomas <- as.data.frame(log2(primary.melanomas))
primary.melanomas.four <- data.frame(primary.melanomas[, 
                      names(primary.melanomas) %in% c("COL1A2")])
primary.melanomas.four$mean_exp <- rowMeans(primary.melanomas.four)
primary.melanomas.targets <- data.frame(primary.melanomas[, 
                      !names(primary.melanomas) %in% c("COL1A2")])
primary.melanomas.targets$mean_exp <- rowMeans(primary.melanomas.targets)
primary.melanomas.targets$mean_exp_four <- primary.melanomas.four$mean_exp
primary.melanomas.targets$COL1A2 <- primary.melanomas$COL1A2

#calculate Pearson correlation coefficient using the above gene (i.e., COL1A2) against ATF4-targets for both primary samples alone and primary + metastatic samples

p_coef_both8 <- cor.test(x = primary.melanomas.targets$COL1A2, y = primary.melanomas.targets$mean_exp, method = "pearson")
pm_coef_both8 <- cor.test(x = primary.metastatic.melanomas.targets$COL1A2, y = primary.metastatic.melanomas.targets$mean_exp, method = "pearson")

#produce and save figures

tiff("<Path to save figure>", units = "in", width = 5, height = 3, res = 300, compression = "lzw")
ggplot(primary.melanomas.targets, aes(x = mean_exp, y = COL1A2)) +
  geom_point(size = 1, color = "cornflowerblue", ) +
  geom_smooth(method = "lm", color = "cornflowerblue") +
  ggtitle(sprintf("Pearson's r = %0.2f, P-value = %0.2e", p_coef_both8$estimate, p_coef_both8$p.value)) + 
  xlab("Random genes (mean log2 expression)") + 
  ylab("COL1A2 (log2 expression)") +
  theme_classic() +
  theme(axis.line.x = element_line(color = "black", size = 0.3),
        axis.line.y = element_line(color = "black", size = 0.3),
        axis.text.x = element_text(color = "black", size = 18),
        axis.text.y = element_text(color = "black", size = 18))
dev.off()

tiff("<Path to save figure>", units = "in", width = 5, height = 3, res = 300, compression = "lzw")
ggplot(primary.metastatic.melanomas.targets, aes(x = mean_exp, y = COL1A2)) +
  geom_point(size = 1, color = "cornflowerblue", ) +
  geom_smooth(method = "lm", color = "cornflowerblue") +
  ggtitle(sprintf("Pearson's r = %0.2f, P-value = %0.2e", pm_coef_both8$estimate, pm_coef_both8$p.value)) + 
  xlab("Random genes (mean log2 expression)") + 
  ylab("COL1A2 (log2 expression)") +
  theme_classic() +
  theme(axis.line.x = element_line(color = "black", size = 0.3),
        axis.line.y = element_line(color = "black", size = 0.3),
        axis.text.x = element_text(color = "black", size = 18),
        axis.text.y = element_text(color = "black", size = 18))
dev.off()

# END OF SCRIPT