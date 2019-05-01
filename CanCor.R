# assess the top 5 most reprogrammed metabolisms and their relationship to
# fenton reactions in the cytosol.


# example data

data(Rohwer, package="heplots")
X <- as.matrix(Rohwer[,6:10])  # the PA tests
Y <- as.matrix(Rohwer[,3:5])   # the aptitude/ability variables

# visualize the correlation matrix using corrplot()
if (require(corrplot)) {
  M <- cor(cbind(X,Y))
  corrplot(M, method="ellipse", order="hclust", addrect=2, addCoef.col="black")
}



(cc <- cancor(X, Y, set.names=c("PA", "Ability")))

plot(cc)

# 
# # simulate dummy data
# 
# dat<- matrix(rnorm(10000, 50, 8), nrow = 100, ncol = 100)
# 
# # visualize the distribution
# plot(density(dat))
# 
# 
# # partition into two groups
# 
# dat <- as.data.frame(dat)
# colnames(dat) <- rep(1, 100) 
# colnames(dat)[51:100] <- rep(0, 50)
# 
# 
# X <- as.matrix(dat[,1:50])  
# Y <- as.matrix(dat[,51:100])  
# 
# # visualize the correlation matrix using corrplot()
# 
# if (require(corrplot)) {
#   M <- cor(cbind(X,Y))
#   corrplot(M, method = "color")
#   
# }
# 
# cc <- cancor(X,Y, set.names = c("1", "0"))
# plot(cc)


#Use on ranked pathways for input into 

projects <-c("TCGA-BRCA", "TCGA-COAD", "TCGA-HNSC", "TCGA-KIRC", "TCGA-KIRP", "TCGA-LUAD", "TCGA-STAD", "TCGA-THCA",
  "TCGA-BLCA", "TCGA-KICH", "TCGA-LIHC", "TCGA-LUSC")  
  
annot <- data.table::fread(
  "~/CSBL_shared/ID_mapping/Ensembl_symbol_entrez.csv"
) %>%
  dplyr::select(Ensembl = ensembl_gene_id, GeneName = external_gene_name)

gene.list <- data.table::fread("~/storage/Metabolic_reprogramming/cleaned_gene_list.csv") %>%
  #filter(str_detect(Pathway, "Fenton", negate = T)) %>%
  dplyr::select(-PathwayType) %>%
  # bind_rows(proteasomes) %>%
  # bind_rows(gangliosides, sialic.acids) %>%
  # Fix multiple mappings in symbol -> Ensembl
  left_join(annot, by = "Ensembl") %>%
  drop_na(GeneName) %>%
  distinct() %>%
  mutate(
    Pathway = str_replace_all(Pathway, "[\\s/-]+", "_"),
    Symbol = case_when(
      Symbol != GeneName ~ GeneName,
      T ~ Symbol
    )) %>%
  dplyr::select(-GeneName)

# I have the genes now I need the expression for these genes in each project. I then can perform a correlation alanlysis 
# on the genes from the pathways. 

pw.list <- unique(gene.list$Pathway)

clinical <- data.table::fread(
  "~/CSBL_shared/RNASeq/TCGA/annotation/counts_annotation.csv")
normal.samples <- clinical[sample_type == "Solid Tissue Normal"]
tumor.samples <- clinical[sample_type != "Solid Tissue Normal"]
#tumor.samples <- tumor.samples[tumor.samples$caseID %in% normal.samples$caseID,]


FPKMtoTPM <- function(x) {
  return(exp(log(x) - log(sum(x)) + log(1e6)))
}

pw <- pw.list[1]

proj <- projects[1]
# read in my TCGA RNA-seq data
for (proj in projects) {
  df.exp <- data.table::fread(
    str_glue("~/CSBL_shared/RNASeq/TCGA/FPKM/{proj}.FPKM.csv")) %>%
    mutate_if(is.numeric, FPKMtoTPM)%>%
    mutate(Ensembl = gsub("\\.\\d+$", "", Ensembl))
   
  # now we have the expression, we need to subset the expression to only the data in the genes we are interested in. 
  df.exp <- df.exp[df.exp$Ensembl %in% gene.list$Ensembl,] 
  
  rownames(df.exp) <- df.exp$Ensembl
    
  # TODO: should we remove the genes with no expression?
  df.exp <- df.exp[rowSums(df.exp == 0) < ncol(df.exp) / 2, ]
  
  # TODO: sample numbers not matching
  df.normal <- df.exp[ ,colnames(df.exp) %in% normal.samples$barcode[normal.samples$project == proj]]
  df.tumor <- df.exp[ ,colnames(df.exp) %in% tumor.samples$barcode[tumor.samples$project == proj]] 
  
  
  # Now we have the genes in our gene.list, we can make a matrix with the expression from the two groups;
  # tumor and normal. We will cut the matrx and cut the for the samples that are T and N.

  # Add columns identifying the pathway and the sample type. This will include a tumor or normal
  # status and the pathway the gene belongs in.
  
  for(pw in pw.list){
    
    # subset the df.normal into the pathways we are comparing
    
    MetabolicPathwayi <- gene.list[gene.list$Pathway ==pw,]
    #Fenton <- gene.list[gene.list$Pathway == "cytosol_Fenton",]
    
    
    # get the proteasome genes 
    
    # get the proton producing genes for each pathway
    
    # make sure no genes overlap between the two groups, drop those genes for comparssion
    
    # continue down to ext step. 
    
    df.normal.X <- df.normal[rownames(df.normal) %in% MetabolicPathwayi$Ensembl,]
    
    df.normal.Y <- df.normal[rownames(df.normal) %in% Fenton$Ensembl,]
    
    M <- cor(cbind(df.normal.X,df.normal.Y))
    p<- corrplot(M, method="ellipse", order="hclust", addrect=2, addCoef.col="black")
   
    ggplot2::ggsave(p, filename = str_glue("~/storage/Metabolic_reprogramming/Figures/Figure_2/{proj}_normal_corplot.png"), device = "png",
                    width = 8, height = 8, units = "in", dpi = "retina", limitsize = F)
    
    cc <- cancor(t(df.tumor.X), t(df.tumor.Y), set.names=c(pw, "Fenton"))
    
    p<- plot(cc)
    
    ggplot2::ggsave(p, filename = str_glue("~/storage/Metabolic_reprogramming/Figures/Figure_2/{proj}_normal_can_cor.png"), device = "png",
                    width = 8, height = 8, units = "in", dpi = "retina", limitsize = F)
    
    
    df.tumor.X <- df.tumor[rownames(df.tumor) %in% MetabolicPathwayi$Ensembl,]
    
    df.tumor.Y <- df.tumor[rownames(df.tumor) %in% Fenton$Ensembl,]
    
    M <- cor(cbind(df.tumor.X,df.tumor.Y))
    p<- corrplot(M, method="ellipse", order="hclust", addrect=2, addCoef.col="black")
    
    ggplot2::ggsave(p, filename = str_glue("~/storage/Metabolic_reprogramming/Figures/Figure_2/{proj}_tumor_corplot.png"), device = "png",
                    width = 8, height = 8, units = "in", dpi = "retina", limitsize = F)
    
    cc <- cancor(X, Y, set.names=c(pw, "Fenton"))
    
    p<- plot(cc)
    
    ggplot2::ggsave(p, filename = str_glue("~/storage/Metabolic_reprogramming/Figures/Figure_2/{proj}_tumor_can_cor.png"), device = "png",
                    width = 8, height = 8, units = "in", dpi = "retina", limitsize = F)
    
    
    
  }
  
  
  
  
}

  
