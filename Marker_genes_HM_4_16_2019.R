# use the marker genes to demonstrate MR
# in our study

library(tidyverse)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(KEGGREST)
library(data.table)
library(ggpubr)
library(clusterProfiler)
library(dplyr)
library(circlize)
library(grid)

# gene.list
# gene.list <- data.table::fread("gene_list_Ying.txt", header = TRUE)
# colnames(gene.list) <- c("Pathway", "Gene")
# 
# gene.list2 <- strsplit(gene.list$Gene, split = ",")
# gene.list3<- data.frame(Pathway = rep(gene.list$Pathway, sapply(gene.list2, length)), Gene = unlist(gene.list2))
# gene.list <- gene.list3
# rm(gene.list2)
# rm(gene.list3)
# # convert the expression
# 
# eg <- bitr(gene.list$Gene, fromType = "SYMBOL", toType = "ENSEMBL", OrgDb = org.Hs.eg.db)
# 
# gene.list <- left_join(gene.list, eg, by =c("Gene" = "SYMBOL"))
setwd("~/storage/Metabolic_reprogramming/Proton_producing_reactions")
gene.list <- data.table::fread("gene_list_makers.txt", header = TRUE)

#projects
projects <- c("TCGA-BLCA","TCGA-BRCA","TCGA-COAD","TCGA-ESCA","TCGA-HNSC","TCGA-KIRC",
              "TCGA-KICH","TCGA-KIRP","TCGA-LIHC","TCGA-LUAD","TCGA-LUSC","TCGA-STAD", 
              "TCGA-THCA", "TCGA-PRAD")


rel.survival <- data.table::fread(
  "tcga_survival.csv"
)
colnames(rel.survival) <- c("Project", "TCGA", "SEER")
rel.survival <- rel.survival %>%
  mutate(
    Project = str_glue("TCGA-{Project}"),
    TCGA = str_extract(TCGA, "\\([\\d.]+"),
    SEER = str_extract(SEER, "\\([\\d.]+")
  ) %>%
  mutate(
    TCGA = as.numeric(str_remove(TCGA, "\\(")),
    SEER = as.numeric(str_remove(SEER, "\\("))
  ) %>%
  filter(Project %in% projects) %>%
  arrange(SEER, TCGA) %>%
  mutate(Project = str_remove(Project, "TCGA-"))


# clincal
clinical <- data.table::fread(
  "~/CSBL_shared/RNASeq/TCGA/annotation/counts_annotation.csv")
normal.samples <- clinical[sample_type == "Solid Tissue Normal"]
tumor.samples <- clinical[sample_type != "Solid Tissue Normal"]

# lets read in data 

proj <- projects[5]
setwd("~/CSBL_shared/RNASeq/TCGA/FPKM")
for(proj in projects) {
  print(proj)
  res <- data.table::fread(str_glue(
    "~/CSBL_shared/RNASeq/TCGA/DEA/tumor_vs_normal/{proj}.csv"),
    data.table = FALSE) %>%
    dplyr::rename(ENSEMBL = V1) %>%
    mutate(ENSEMBL = str_remove(ENSEMBL, "\\..*$")) %>%
    filter(padj < 0.01) %>%
    dplyr::select(ENSEMBL, log2FoldChange)
  
  gene.list <- left_join(gene.list, res, by =c("Ensembl" ="ENSEMBL"))
  
}

#   filter(padj < 0.01) %>%
#   dplyr::select(ENSEMBL, log2FoldChange) %>%
#   mutate(log2FoldChange = 2 ^ log2FoldChange)


colnames(gene.list) <- c("Pathway", "Gene", "ENSEMBL","BLCA","BRCA","COAD","ESCA","HNSC",
                         "KIRC","KICH","KIRP","LIHC","LUAD","LUSC","STAD","THCA","PRAD")

# replace NA with 0

gene.list[is.na(gene.list)] <- 0

# plot FC
FC <- gene.list
FC$Gene <- NULL


long_FC <- FC %>%
  gather(key = "Project", value = "FC", -c(ENSEMBL, Pathway)) %>%
  group_by(Pathway, Project) %>%
  summarise(
    meanFC = mean(FC)
  )
Wide_mean_FC <- long_FC %>%
  spread(key = Project, value = meanFC) %>%
  column_to_rownames("Pathway") %>%
  dplyr::select(one_of(rel.survival$Project))


HeatmapColors <- function() {
  colors <- RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")
  colors[4] <- "#FFFFFF"
  
  colorRampPalette(rev(colors))(100)
}

ComplexHeatmap::Heatmap(
  Wide_mean_FC,
  col = circlize::colorRamp2(c(-4, 0, 4), c("#4575B4", "white", "#D73027")),
  cluster_rows = F, cluster_columns = F, show_row_names = T,
  heatmap_legend_param = list(title = "FC"),
  column_names_gp = grid::gpar(fontsize = 8),
  row_names_gp = grid::gpar(fontsize = 8),
  show_row_dend = FALSE
)



df <- data.table::fread("~/CSBL_shared/RNASeq/TCGA/FPKM/TCGA-KICH.FPKM.csv")
df <- as.data.frame(df)
df <- df %>% column_to_rownames("Ensembl")
df <- FPKM_to_TPM(df)
df <- df %>%
  rownames_to_column(var = "Ensembl")

df2 <- data.table::fread("~/CSBL_shared/RNASeq/TCGA/FPKM/TCGA-HNSC.FPKM.csv")
df2 <- as.data.frame(df2)
df2 <- df2 %>% column_to_rownames("Ensembl")
df2 <- FPKM_to_TPM(df2)
df2 <- df2 %>%
  rownames_to_column(var = "Ensembl")




df$Ensembl <- substr(df$Ensembl, 1, 15)
df2$Ensembl <- substr(df2$Ensembl, 1, 15)

df <- filter(df, Ensembl == "ENSG00000110090")
df2 <- filter(df2, Ensembl == "ENSG00000110090")


clinical <- data.table::fread(
  "~/CSBL_shared/RNASeq/TCGA/annotation/counts_annotation.csv")
normal.samples <- clinical[sample_type == "Solid Tissue Normal"]
tumor.samples <- clinical[sample_type != "Solid Tissue Normal"]


rowMeans(df[ ,colnames(df) %in% normal.samples$barcode[normal.samples$project == "TCGA-KICH"]])
rowMeans(df[ ,colnames(df) %in% tumor.samples$barcode[tumor.samples$project == "TCGA-KICH"]])

rowMeans(df2[ ,colnames(df2) %in% normal.samples$barcode[normal.samples$project == "TCGA-HNSC"]])
rowMeans(df2[ ,colnames(df2) %in% tumor.samples$barcode[tumor.samples$project == "TCGA-HNSC"]])


