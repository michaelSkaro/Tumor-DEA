library(tidyverse)
library(psych)
library(clusterProfiler)
library(org.Hs.eg.db)



projects <- c("TCGA-BLCA","TCGA-BRCA","TCGA-COAD","TCGA-ESCA","TCGA-HNSC","TCGA-KICH","TCGA-KIRC","TCGA-KIRP","TCGA-LIHC","TCGA-LUAD","TCGA-LUSC","TCGA-PRAD","TCGA-STAD","TCGA-THCA")


# okay I need the gene lists for the FR and the metabolisms

gene.list <- data.table::fread("/home/mskaro1/storage/Co_expression_networks_for_fenton_reactions_in_cancer/cleaned_gene_list.csv") %>%
  as_tibble() %>%
  filter(Pathway != "UPR") %>%
  mutate(PathwayType = case_when(
    Pathway == "cytosol Fenton" ~ "Cytosolic Fenton",
    Pathway == "mitochondria Fenton" ~ "Mitochondrial Fenton",
    T ~ "Metabolic"
  ))

View(gene.list)

annot <- data.table::fread("~/CSBL_shared/ID_mapping/Ensembl_symbol_entrez.csv") %>%
  dplyr::select(Ensembl = ensembl_gene_id, Genes = external_gene_name) %>%
  filter(Genes %in% gene.list$Symbol)

clinical <- data.table::fread(
  "~/CSBL_shared/RNASeq/TCGA/annotation/counts_annotation.csv")
normal.samples <- clinical[sample_type == "Solid Tissue Normal"]
tumor.samples <- clinical[sample_type != "Solid Tissue Normal"]
tumor.samples <- tumor.samples[tumor.samples$caseID %in% normal.samples$caseID,]
#rm(clinical)

FPKMtoTPM <- function(x) {
  return(exp(log(x) - log(sum(x)) + log(1e6)))
}

proj <- projects[1]
for (proj in projects) {
  df.exp <- data.table::fread(
    str_glue("~/CSBL_shared/RNASeq/TCGA/FPKM/{proj}.FPKM.csv")) %>%
    mutate_if(is.numeric, FPKMtoTPM) %>%
    mutate(Ensembl = gsub("\\.\\d+$", "", Ensembl)) %>%
    inner_join(annot, by = "Ensembl") %>%
    distinct(Genes, .keep_all = T) %>%
    column_to_rownames("Genes") %>%
    dplyr::select(-Ensembl)
  
  # TODO: should we remove the genes with no expression?
  df.exp <- df.exp[rowSums(df.exp == 0) < ncol(df.exp) / 2, ]
  
  # TODO: sample numbers not matching
  df.normal <- df.exp[ ,colnames(df.exp) %in% normal.samples$barcode[normal.samples$project == proj]]
  df.tumor <- df.exp[ ,colnames(df.exp) %in% tumor.samples$barcode[tumor.samples$project == proj]] 
  
  
  
  df.normal <- df.normal %>%
    rownames_to_column("Genes") %>%
    right_join(gene.list, by = c("Genes" = "Symbol")) %>%
    drop_na() %>% distinct(Genes, .keep_all = T)
  
  df.tumor <- df.tumor %>%
    rownames_to_column("Genes") %>%
    right_join(gene.list, by = c("Genes" = "Symbol")) %>% 
    drop_na() %>% distinct(Genes, .keep_all = T)
  
  dfn.cor <- cor(t(df.normal[, str_starts(colnames(df.normal), "TCGA")]),
                 method = "spearman")
  dft.cor <- cor(t(df.tumor[, str_starts(colnames(df.tumor), "TCGA")]),
                 method = "spearman")
  
  
  p <- pheatmap::pheatmap(dfn.cor, cluster_rows = T, cluster_cols = T)
  ggplot2::ggsave(p, filename = str_glue("./{proj}_normal_heatmap_pathway_4_11_2019.png"), device = "png",
                  width = 60, height = 40, units = "in", dpi = "retina", limitsize = F)
  
  order_i_need <- p$tree_row$order
  
  dft.cor.new <- dft.cor[order_i_need, order_i_need]
  
  p <- pheatmap::pheatmap(dft.cor.new, cluster_rows = F, cluster_cols = F)
  ggplot2::ggsave(p, filename = str_glue("./{proj}_tumor_heatmap_pathway_4_11_2019.png"), device = "png",
                  width = 60, height = 40, units = "in", dpi = "retina", limitsize = F)
  
  
  
  
  plot.avg.exp <- function(df, alpha.level = 0.9, suffix) {
    ymax <- ceiling(max(df$MeanNormCounts))
    p <- df.exp %>%
      ggplot(aes(x = tumor_stage, y = MeanNormCounts,
                 group = Pathway, color = Pathway, label = Pathway))+
      geom_point()+
      geom_line()+
      ggtitle(proj)+
      labs(
        x = "Tumor Stage",
        y = "log(Normalized Counts)"
      )+
      # Plot labels in the right margin
      coord_cartesian(clip = 'off')+
      geom_label_repel(
        aes(
          label = ifelse(tumor_stage == "IV", Pathway, ''),
          fill = Pathway, color = Pathway, alpha = alpha.level
        ),
        fontface = 'bold', color = 'white',
        segment.color = 'grey50',
        direction = "y",
        xlim = c(4.2, 20),
        ylim = c(NA, ymax)
      )+
      theme_minimal()+
      theme(
        legend.position = "none",
        plot.margin = unit(c(1,7,0,0), "in")
      )
    
    ggsave(p, file = str_glue(
      "~/storage/Metabolic_reprogramming/average_pathway_expression/plot/{proj}_{suffix}.png"),
      device = "png", width = 20, height = 20,
      units = "in", dpi = "retina")
  }
  
  proj <- projects[2]
  for (proj in projects) {
    df <- data.table::fread(str_glue(
      "~/storage/data/metabolic_reprogramming/average_pathway_exp/TCGA-BRCA.csv"
    )) %>%
      mutate(MeanNormCounts = log(MeanNormCounts))
    
    pathways <- sort(unique(gene.list$Pathway))
    middle.index <- ceiling(length(pathways) / 2)
    
    pw1 <- pathways[1:middle.index]
    pw2 <- pathways[(middle.index+1):length(pathways)]
    
    plot.avg.exp(df[df$Pathway %in% pw1, ], suffix = "1", alpha.level = 1)
    plot.avg.exp(df[df$Pathway %in% pw2, ], suffix = "2", alpha.level = 1)
  }
  
  
}




