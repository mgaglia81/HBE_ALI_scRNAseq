library(dplyr)
library(Seurat)
library(cowplot)
library(ggplot2)

HBEcells <- readRDS("~/Desktop/HBE_single_cell/mt25_nfeat12k/seurat_obj_integrated_percentIAV_celltypes_assigned_mt25_nfeat12k_updated.RDS")

HBEcells@meta.data$flu_infection <- "Uninfected"
HBEcells@meta.data$flu_infection[HBEcells@meta.data$percent.IAV <= 1] <- "Uninfected"
HBEcells@meta.data$flu_infection[HBEcells@meta.data$percent.IAV > 1] <- "Infected"

HBEcells$infection.condition <- paste(HBEcells$flu_infection, HBEcells$Condition, sep="_")

HBEcells_1dpi <- subset(HBEcells, DPI == 1)
HBEcells_3dpi <- subset(HBEcells, DPI == 3)

Idents(HBEcells_1dpi) <- "infection.condition"
Idents(HBEcells_1dpi) <- "flu_infection"
VlnPlot(HBEcells_1dpi, features = ("IAV-HA"))

#genes to plot
MHCgenesall <- c("B2M","CALR","CANX","CD4","CD74","CD8A","CD8B","CIITA","CREB1","CTSB","CTSS","HLA-A","HLA-B","HLA-C","HLA-DMA","HLA-DMB","HLA-DOA","HLA-DOB","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQA2","HLA-DQB1","HLA-DQB2","HLA-DRA","HLA-DRB1","HLA-DRB3","HLA-DRB4","HLA-DRB5","HLA-E","HLA-F",
              "HLA-G","HSP90AA1","HSP90AB1","HSPA5","IFI30","IFNA1","IFNA10","IFNA13","IFNA14","IFNA16","IFNA17","IFNA2","IFNA21","IFNA4","IFNA5","IFNA6","IFNA7","IFNA8","KIR2DL1","KIR2DL2","KIR2DL3","KIR2DL4","KIR2DL5A","KIR2DS1","KIR2DS2","KIR2DS3","KIR2DS4","KIR2DS5","KIR3DL1","KIR3DL2",
              "KIR3DL3","KLRC1","KLRC2","KLRC3","KLRC4","KLRD1","LGMN","LTA","NFYA","NFYB","NFYC","PDIA3","PSME1","PSME2","RFX5","RFXANK","RFXAP","TAP1","TAP2","TAPBP")

IFNgenes <- c("IFNB1", "IFNL1", "IFNL2", "IFNL3", "IFNL4")

sorted_schogginsISG <- c("IFIT2","B2M","IFIT3","ISG15","CXCL10","MX1","IFI6","IFIT1","STAT1","CXCL11","IFITM1","IFI27","MX2","CFB","IFIH1","XAF1","TNFSF10","HERC6",
                         "EXT1","IFITM3","HLA-C","IFI16","IL1RN","SLFN5","LY6E","TXNIP","OAS2","EPSTI1","TYMP","ADAR","IFI44L","S100A8",
                         "EIF2AK2","BST2","HLA-E","GBP1","OASL","EPAS1","LAP3","PLSCR1","ETV6","SP110","PTMA","CD74","SMAD3","TIMP1","IFI44","OPTN",
                         "FNDC3B","SCARB2","GBP4","ISG20","C15orf48","SAMHD1","JAK2","APOL1","MCL1","DTX3L","PNPT1","LIPA","EHD4","GBP3","IFNGR1","LAMP3",
                         "PML","SAA1","CD38","IL6ST","IDO1","ELF1","TAP1","JUNB","APOL6","BAG1","UBE2L6","TNFAIP3","DYNLT1","LGALS3","CDKN1A","TRIM38",
                         "TRIM5","LGALS9","NMI","C1S","STAT2","IRF1","CHMP5","IFI30","MYD88","IFI35","TAP2","TRIM14","NCOA3","HLA-F","B4GALT5","PMAIP1",
                         "APOBEC3A","DDX3X","PARP12","SERPING1","MAP3K5","GLRX","DUSP5","TCF7L2","GALNT2","PLEKHA4","CXCL1","IFITM2","ARHGEF3","ACSL1","TNFSF13B","ANKFY1",
                         "HSH2D","TRIM25","GBP2","IRF7","NPAS2","PSMB8","APOL2","HK2","USP18","IFIT5","AIM2","CD9","TRAFD1","GBP5","PSMB9","NAPA",
                         "CYP1B1","PHF11","PNRC1","CX3CL1","IRF9","RBCK1","TDRD7","MAX","SLC16A1","DCP1A","SLC25A28","GAK","IRF2","SPTLC2","GTPBP1","RIPK2",
                         "IGFBP2","SECTM1","MARCKS","NUP50","LMO2","TRIM21","SAMD4A","GPX2","LGMN","RNF19B","GCA","PPM1K","BCL2L14","ODC1","TNFRSF10A","TLK2",
                         "STEAP4","BATF2","TLR3","CD274","ATF3","PI4K2B","CCND3","RTP4","PRKD2","ZBP1","MASTL","CASP7","AHNAK2","GCH1","PMM2","STARD5",
                         "ETV7","SPSB1","TMEM51","RPL22","EIF3L","CEBPD","DDIT4","SLC15A3","KIAA0040","MAFF","GMPR","MAFB","ADM","VEGFC","DHX58","UNC93B1",
                         "ALDH1A1","GTPBP2","SSBP3","OGFR","LRG1","PFKFB3","C4orf33","ABTB2","CCL5","BCL3","SCO2","ERLIN1","RAB27A","IL15RA","SLC1A1","ULK4",
                         "HES4","PIM3","LEPR","BLVRA","MT1X","TMEM140","FBXO6","DEFB1","CPT1A","MAP3K14","PXK","MTHFD2L","BTN3A3","HLA-G","HPSE","SIRPA",
                         "HEG1","ANKRD22","GK","ARG2","SERPINB9","VAMP5","SLC25A30","SOCS1","FKBP5","SOCS2","IL15","CRY1","TRIM34","CCNA1","CCDC92","PDGFRL",
                         "TAGAP","COMMD3","ABLIM3","IMPA2","NFIL3","HESX1","TNFAIP6","SERPINE1","CXCL9","PUS1","PDK1","FFAR2","THBD","ANGPTL1","TBX3","NCF1",
                         "CLEC2B","SNN","UPP2","CCL2","AKT3","MKX","BUB1","NDC80","P2RY6","FUT4","MICB","STAP1","NOD2","ZNF385B","CCL8","PADI2",
                         "MT1G","MCOLN2","RNASE4","MT1M","ENPP1","AMPH","MT1F","ADAMDEC1","IL17RB","CD80","MSR1","CD69","GEM","RASSF4","CCL4","CCR1",
                         "FNDC4","MS4A4A","FCGR1A","NRN1","CREB3L3","GJA4","AQP9","CES1","MT1H","MAB21L2","CTCFL","RGS1","CRP","CD163")


# Create a new column called "new_name" for infection status and condition

HBEcells_1dpi$new_name <- "NA"
HBEcells_1dpi$new_name <- paste(HBEcells_1dpi$Condition, HBEcells_1dpi$flu_infection, sep="_")
HBEcells_1dpi$new_name[HBEcells_1dpi$new_name == "Mock_Infected"] <- "Mock"
HBEcells_1dpi$new_name[HBEcells_1dpi$new_name == "Mock_Uninfected"] <- "Mock"

#checklevels
VlnPlot(HBEcells_1dpi, features = ("IAV-HA"), group.by = "new_name")

#infected vs uninfected

cell.averages1 <- AggregateExpression(HBEcells_1dpi, return.seurat = T, group.by = c('new_name'))
cell.averages1$names <- rownames(cell.averages1@meta.data)
cell.averages1 <- ScaleData(cell.averages1, features = rownames(cell.averages1))

Idents(cell.averages1) <- "names"

levels(cell.averages1)
# Reorder levels manually
my_levels = c('Mock','WT-Infected','X-Infected','WT-Uninfected','X-Uninfected')
cell.averages1@active.ident <- factor(x = cell.averages1@active.ident, levels = my_levels)

DoHeatmap(cell.averages1, features = MHCDEGs_1DPI, draw.lines = F, group.colors = c("black", "darkolivegreen1", "darkolivegreen4", "plum3", "purple4")) +  scale_fill_gradient2(
       low = rev(c("darkblue",'#67a9cf')),
       mid = "black",
       #high = rev(c('#ef8a62','#b2182b',"firebrick4")),
         high = rev(c('gold',"orange3")),
       midpoint = 0,
       guide = "colourbar",
       aesthetics = "fill", 
       na.value = "white")

#manually ordered 
MHCDEGs_1DPI <- c("CALR","HLA-E","HSP90AB1","PSME1","HSP90AA1","PDIA3","TAPBP","CANX","LGMN","PSME2","TAP1","CD74","HLA-A","HLA-B","HLA-C","B2M","HLA-DRA","HLA-F","CREB1")


## graph for 3 dpi
HBEcells_3dpi$new_name <- paste(HBEcells_3dpi$Condition, HBEcells_3dpi$flu_infection, sep="_")
HBEcells_3dpi$new_name[HBEcells_3dpi$new_name == "Mock_Infected"] <- "Mock"
HBEcells_3dpi$new_name[HBEcells_3dpi$new_name == "Mock_Uninfected"] <- "Mock"
VlnPlot(HBEcells_3dpi, features = ("IAV-HA"), group.by = "new_name")

cell.averages3 <- AggregateExpression(HBEcells_3dpi, return.seurat = T, group.by = c('new_name'))
cell.averages3$names <- rownames(cell.averages3@meta.data)
cell.averages3 <- ScaleData(cell.averages3, features = rownames(cell.averages3))

#manually ordered 
MHCDEGs_3DPI <- c("HLA-F","HLA-B","HLA-A","HLA-C","TAP1","CALR","CD74","PDIA3","HLA-E","IFI30","CANX","HSP90AA1","PSME2","PSME1","HSP90AB1","HLA-G")

levels(cell.averages3)
# Reorder levels manually
my_levels = c('Mock','WT-Infected','X-Infected','WT-Uninfected','X-Uninfected')
cell.averages3@active.ident <- factor(x = cell.averages1@active.ident, levels = my_levels)

DoHeatmap(cell.averages3, features = MHCDEGs_3DPI, draw.lines = F, group.colors = c("black", "darkolivegreen1", "darkolivegreen4", "plum3", "purple4")) +  scale_fill_gradient2(
  low = rev(c("darkblue",'#67a9cf')),
  mid = "black",
  #high = rev(c('#ef8a62','#b2182b',"firebrick4")),
  high = rev(c('gold',"orange3")),
  midpoint = 0,
  guide = "colourbar",
  aesthetics = "fill", 
  na.value = "white")
