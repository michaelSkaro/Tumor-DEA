# this code will be used to generate a bird's 
# eye view of the metabolic reporgramming in 
# cancer with repsect to each gene product's 
# contribution to intracellular pH
# I will mkae a heatmap of each of the proton producing


library(tidyverse)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(heplots)
library(candisc)
library(KEGGREST)
library(data.table)
library(ggpubr)
library(pcxn)
library(vegan)
library(clusterProfiler)
library(dplyr)

# proton producing reactions in our back ground
dat <- data.table::fread("~/storage/Metabolic_reprogramming/uniprot_proton_rxn_in_MR.csv", data.table = FALSE)

# convert these Symbols to the ENSEMBL ids
eg <- bitr(dat$Symbol, fromType = "SYMBOL", toType = "ENSEMBL", OrgDb = org.Hs.eg.db, drop = TRUE)

# join
dat <- left_join(dat, eg, by = c("Symbol" = "SYMBOL"))


# fold change of each of these genes between normal and tumor
setwd("~/CSBL_shared/RNASeq/TCGA/DEA/tumor_vs_normal")
DESeq_res_objects<- list.files()
DESeq_res_objects <- DESeq_res_objects[2:15]

proj <- DESeq_res_objects[1]
for(proj in DESeq_res_objects){
  
  res <- data.table::fread(str_glue("{proj}"), data.table = FALSE)
  colnames(res) <- c("ENSEMBL", "baseMean", "log2FoldChange","lfcSE","stat", "pvalue","padj")
  res <- res[res$padj <1e-2,]
  res$ENSEMBL <- substr(res$ENSEMBL, 1, 15)
  
  res <- res %>%
    dplyr::select(ENSEMBL, log2FoldChange)
  
  res$log2FoldChange <- 2^(res$log2FoldChange)
 # merge table with current cancer type add column with  
  
  dat2 <- left_join(dat, res, by ="ENSEMBL")
  
  dat <- dat2
  
}

colnames(dat) <- c("Symbol", "NetProton", "Uniprot", "Reaction", "ENSEMBL", 
                   "BLCA","BRCA","COAD","ESCA","HNSC","KICH",
                   "KIRC","KIRP","LIHC","LUAD","LUSC","PRAD",
                   "STAD","THCA") 

dat <- dat %>%
  mutate_all(~replace(., is.na(.), 1))

# add pathways as a column to combine the genes
setwd("~/storage/Metabolic_reprogramming")
y <- data.table::fread("cleaned_gene_list.csv", data.table = FALSE)

dat <- dat %>%
  left_join(y, by = c("ENSEMBL" = "Ensembl"))

completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}


dat <- completeFun(dat, "Pathway")
dat <- dat[order(dat$Pathway),]

FC <- dat[,5:20]
long_FC <- FC %>%
  gather(key = "Project", value = "FC", -c(ENSEMBL, Pathway)) %>%
  group_by(Pathway, Project) %>%
  summarise(
    meanFC = mean(FC)
  )
Wide_mean_FC <- long_FC %>%
  spread(key = Project, value = meanFC)

rownames(Wide_mean_FC) <- Wide_mean_FC$Pathway
Wide_mean_FC$Pathway <- NULL

Wide_mean_FC <- as.matrix(Wide_mean_FC)

# plotting time

ComplexHeatmap::Heatmap(t(Wide_mean_FC), col = colorRamp2(c(0, 2, 4), c("blue", "white", "red")), 
                        cluster_rows = TRUE, cluster_columns = FALSE, 
                        show_row_names = TRUE, heatmap_legend_param = list(title = "FC"), 
                        column_names_gp = gpar(fontsize = 8), row_names_gp = gpar(fontsize = 8))



