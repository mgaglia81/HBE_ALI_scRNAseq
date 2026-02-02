library(dplyr)
library(Seurat)
library(cowplot)
library(ggpubr)


s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

cell.combined.cc <- CellCycleScoring(HBEcells, s.features = s.genes, g2m.features = g2m.genes, set.ident = T)

cell.combined.cc1 <- subset(cell.combined.cc, DPI == '1')

Mock1.cc <- subset(cell.combined.cc1, Condition == 'Mock')
WT1.cc <- subset(cell.combined.cc1, Condition == 'WT')
X1.cc <- subset(cell.combined.cc1, Condition == 'X')


phase.list <- c('G1','S','G2M')

#relevel factors -- cluster
# View current levels
levels(cell.combined.cc1@meta.data$seurat_clusters)
# Reorder levels manually
cell.combined.cc1@meta.data$seurat_clusters <- factor(cell.combined.cc1@meta.data$seurat_clusters, levels = c('7','8','0','3','9','12','2','11','5','10','14','13','1','4','6'))
# View new levels
levels(cell.combined.cc1@meta.data$seurat_clusters)

install.packages("tidyverse")
library(tidyverse)
library(dplyr)
CC <-table(cell.combined.cc1@meta.data$seurat_clusters, cell.combined.cc1@meta.data$Condition, cell.combined.cc1@meta.data$Phase)
CCdf <- as.data.frame(CC)

#cluster = Var1, Condition = Var2, Phase = Var3

levels(CCdf$Var3)
CCdf$Var3 <- factor(CCdf$Var3, levels = c('G1', 'S', 'G2M'))

percentCCdf  <- CCdf   %>%
  group_by(Var2, Var1, Var3) %>%
  summarise(n = Freq) %>%
  mutate(percentage = 100*(n / sum(n)))

p1 <- ggplot(percentCCdf, aes(x = Var2, y = percentage, fill = Var3, group_by(seurat_clusters)))+
  geom_bar(stat = "identity")+
  facet_wrap(~Var1) +
  geom_text(aes(label = paste(round(percentage,2),"%")), position = position_stack(vjust =  0.5)) +
  theme_classic() +
  #scale_fill_manual(values = c("G1" = "firebrick", "S" = "seagreen3", "G2M" = "lightblue1"))
  scale_fill_manual(values = c("G1" = "gray30", "S" = "#00a69c", "G2M" = "lightblue1"))

#day1 by condition and infection

HBEcells@meta.data$flu_infection[HBEcells@meta.data$percent.IAV <= 1] <- "Uninfected"
HBEcells@meta.data$flu_infection[HBEcells@meta.data$percent.IAV > 1] <- "Infected"

cell.combined.cc <- CellCycleScoring(HBEcells, s.features = s.genes, g2m.features = g2m.genes, set.ident = T)
cell.combined.cc1 <- subset(cell.combined.cc, DPI == '1')

Mock1.cc <- subset(cell.combined.cc1, Condition == 'Mock')
WT1.cc <- subset(cell.combined.cc1, Condition == 'WT')
X1.cc <- subset(cell.combined.cc1, Condition == 'X')

CC <-table(cell.combined.cc1@meta.data$flu_infection, cell.combined.cc1@meta.data$Condition, cell.combined.cc1@meta.data$Phase)
CCdf <- as.data.frame(CC)
CCdf$Var3 <- factor(CCdf$Var3, levels = c('G1', 'S', 'G2M'))
CCdf$Var1 <- factor(CCdf$Var1, levels = c('Uninfected', 'Infected'))

percentCCdf  <- CCdf   %>%
  #group_by(DPI, Donor, Condition, cluster) %>%
  group_by(Var2, Var1, Var3) %>%
  summarise(n = Freq) %>%
  mutate(percentage = 100*(n / sum(n)))


p1 <- ggplot(percentCCdf, aes(x = Var2, y = percentage, fill = Var3, group_by(seurat_clusters)))+
  geom_bar(stat = "identity")+
  facet_wrap(~Var1) +
  geom_text(aes(label = paste(round(percentage,2),"%")), position = position_stack(vjust =  0.5)) +
  theme_classic() +
  #scale_fill_manual(values = c("G1" = "firebrick", "S" = "seagreen3", "G2M" = "lightblue1"))
  scale_fill_manual(values = c("G1" = "gray30", "S" = "#00a69c", "G2M" = "lightblue1"))

