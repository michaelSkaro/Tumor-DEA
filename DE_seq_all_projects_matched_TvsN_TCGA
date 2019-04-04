# Differenital expression between primary tumors of each tissue and the 
# primary tumors that develop in the location the primary lesions metastesize towards

# In this investigation we will coduct three analyses
  # 1:Tumor vs Tumor, this will establish that there is a bonfide difference between the two primary lesions
  # 2: Primary vs. Normal (Primary loc) VS Primary vs Normal (Met loc)

# tumor vs tumor DEseq2

# read in the data for each tumor type

library(DESeq2)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(WGCNA)




projects <- c("TCGA-BLCA","TCGA-BRCA","TCGA-COAD","TCGA-ESCA","TCGA-HNSC","TCGA-KICH","TCGA-KIRC","TCGA-KIRP","TCGA-LIHC","TCGA-LUAD","TCGA-LUSC","TCGA-PRAD","TCGA-STAD","TCGA-THCA")

annot <- data.table::fread("~/CSBL_shared/ID_mapping/Ensembl_symbol_entrez.csv")

clinical <- data.table::fread(
  "~/CSBL_shared/RNASeq/TCGA/annotation/counts_annotation.csv")
normal.samples <- clinical[sample_type == "Solid Tissue Normal"]
tumor.samples <- clinical[sample_type != "Solid Tissue Normal"]
tumor.samples <- tumor.samples[tumor.samples$caseID %in% normal.samples$caseID,]


proj <- projects[1]

for (proj in projects) {
  
  df.exp <- data.table::fread(str_glue("~/CSBL_shared/RNASeq/TCGA/counts/{proj}.counts.csv")) %>%
    as_tibble()
  rownames(df.exp) <- df.exp$Ensembl
  
  coldata.t <- tumor.samples[tumor.samples$project == proj,]
  coldata.n <- normal.samples[normal.samples$project == proj,]
  
  coldata <- rbind(coldata.n, coldata.t)
  rownames(coldata) <- coldata$barcode
  
  df.exp <- df.exp[ ,colnames(df.exp) %in% coldata$barcode]
  
  
  dds <- DESeqDataSetFromMatrix(countData = df.exp, colData = coldata, design = ~ sample_type)
  
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  
  
  dds$sample_type <- relevel(dds$sample_type, ref = "Solid Tissue Normal")
  
  
  dds <- DESeq(dds)
  res <- results(dds)
  resOrdered <- res[order(res$pvalue),]
  
  write.csv("stuffs.txt")
}



# using the results from the tumor vs. normal pathway enrichment we will look at the DE gene lists and create venn diagrams to disern the deifferences

