## Figures of mutations distribution in chromosomes
library(dplyr)
library(ggplot2)
gs_list <- read.table("filt258.tabular", header = T)
genes <- c("ADA", "ABCB1", "ABCG2", "ASXL1", "BCOR", "BCORL1",
           "BRAF", "CBL", "CDA", "CDKN2A", "CEBPA", "CSF3R",
           "CSFR1", "DAXX", "DCK", "DCLK1", "DIS3", "DNMT3A",
           "ETV6", "EZH2", "FAM5C", "FBXW7", "FLT3", "GATA1",
           "GATA2", "GATA3", "HNRNPK", "HRAS", "IDH1", "IDH2",
           "IL7R", "JAK1", "JAK2", "JAK3", "KDM6A", "KMT2A",
           "KIT", "KRAS", "MIR-142", "MPL", "MYD88", "NOTCH1",
           "NPM1", "NRAS", "NT5C2", "PHF6", "PTEN", "PTPN11",
           "PTPRT", "RAD21", "RUNX1", "SETBP1", "SF1", "SF3A1", "SF3B1",
           "SMC1A", "SMC3", "SRSF2", "STAG2", "TERC", "TERT", "TET2", 
           "TP53", "U2AF1", "U2AF2", "WAC", "WT1", "ZRSR2")

newGs2 <- gs_list %>% filter(gene %in% genes)

newGs2 <- newGs2 %>% mutate(chr = gsub("chr", "", chrom), chr = ifelse(chr == "X", "23", chr), chr = as.character(chr), chr = factor(chr, levels = sort(unique(as.numeric(chr)))))

tmp <- newGs2 %>% group_by(chr, gene) %>% count() %>% ungroup()

tmp <- tmp %>% arrange(chr, gene) %>% mutate(Gene = paste('chr', chr, ',', gene), Gene = factor(Gene, Gene))

library(RColorBrewer)

brewer.pal(n=12, name = "Paired")
col1 <- c("#6A3D9A", "#FFFF99", "#CAB2D6", "#FDBF6F", "#33A02C", "#FB9A99")
tmp$col <- sapply(tmp$colindex, function(e){col1[e]})

## Plot the mutations distribution on chromosomes
xy <- ggplot(data = tmp, aes(x = chr, y = n)) + geom_bar(stat = "identity", aes(fill = Gene)) + 
  scale_fill_manual(values = tmp$col) +
  scale_x_discrete(breaks = 1:23, labels = paste0('chr', c(1:22, 'X'))) +
  ggtitle("Distribution of mutations in chromosomes after filtering for 68 genes") + xlab("Chromosomes in order") + ylab("Number of mutations") +
  theme_bw()
pdf('tmpg.pdf', height = 10, width = 6)
print(xy)
dev.off()


## Arrange genes in descending order
tmp1 <- tmp %>% arrange(gene) %>% mutate(gene = factor(gene, gene[order(n, decreasing = T)]))

## Plot the distribution of mutations in the 68 most recurrently mutated genes in AML 
ggplot(data = tmp1, aes(x = gene, y = n)) + geom_bar(stat = "identity", fill=colsC2) +
  # theme_pubr() +
  theme(axis.text.x = element_text(angle = 90, size = 8, vjust = 0)) + ggtitle("Distribution of mutations in the 68 genes") + xlab("Gene names") + ylab("Number of mutations") 
