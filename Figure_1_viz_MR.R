setwd("~/CSBL_shared/RNASeq/TCGA/DEA/tumor_vs_normal")
library(org.Hs.eg.db)
library(clusterProfiler)
library(wesanderson)
library(ggrepel)


#Figure 1A

# purine vs. pyrimidine de novo synthesis. 
X <- c("GART", "PAICS", "PFAS")
X <- as.data.frame(X)
colnames(X) <- "Genes"
Y <- c("CAD", "CTPS1", "DTYMK")
Y <- as.data.frame(Y)
colnames(Y) <- "Genes"
gene.list <- rbind(X,Y)

eg = bitr(gene.list$Genes, fromType="SYMBOL", toType="ENSEMBL", OrgDb="org.Hs.eg.db")

eg$Pathway <- "Purine"
eg$Pathway[5:7] <- "Pyrimidine"

# expression code




# i need the expression FC for each gene for each cancer type. This can be found in the DEA folder
projects <- c("TCGA-BLCA.csv", "TCGA-BRCA.csv", "TCGA-COAD.csv", "TCGA-ESCA.csv", "TCGA-HNSC.csv", "TCGA-KICH.csv", "TCGA-KIRC.csv", "TCGA-KIRP.csv", "TCGA-LIHC.csv", "TCGA-LUAD.csv", "TCGA-LUSC.csv", "TCGA-PRAD.csv", "TCGA-STAD.csv", "TCGA-THCA.csv")

proj <- projects[1]
i <- 1
PathwayMeans <- as.data.frame(matrix(NA, 14, 2))
colnames(PathwayMeans) <- c("Purine", "Pyrimidine")
counter <- 1

for(proj in projects){
  
  dat <- data.table::fread(proj, data.table = FALSE)
  names <- c("ENSEMBL", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
  colnames(dat) <- names
  dat$ENSEMBL <- substr(dat$ENSEMBL, 1, 15)
  df <- dat[dat$ENSEMBL %in% eg$ENSEMBL,]
  PurineMeanFC <- mean(2^(df$log2FoldChange[1:3]))
  PyrimidineMeanFC <- mean(2^(df$log2FoldChange[4:6]))

  
  PathwayMeans$Purine[counter] <- PurineMeanFC
  PathwayMeans$Pyrimidine[counter] <- PyrimidineMeanFC
  
  counter <- (counter+1)
  
}

PathwayMeans$CT <- substring(projects, 6, 9)

# plotting code

Figure_1A <- ggplot2::ggplot(PathwayMeans, aes(Purine, Pyrimidine, label = CT)) +
  geom_point(color ="black") + geom_text_repel(size =3) + geom_abline(intercept = 0, slope = 1) + theme_classic() + xlim(0, 4) + ylim(0, 6) 

#Figure 1B

# purine vs. pyrimidine salvage synthesis. 
X <- c("ATIC", "DGUOK", "HPRT1")
Y <- c("CDA", "DCK", "TK1", "TK2")

X <- as.data.frame(X)
colnames(X) <- "Genes"
Y <- as.data.frame(Y)
colnames(Y) <- "Genes"
gene.list <- rbind(X,Y)

eg = bitr(gene.list$Genes, fromType="SYMBOL", toType="ENSEMBL", OrgDb="org.Hs.eg.db")

eg$Pathway <- "Purine"
eg$Pathway[4:7] <- "Pyrimidine"

# expression code




# i need the expression FC for each gene for each cancer type. This can be found in the DEA folder
projects <- c("TCGA-BLCA.csv", "TCGA-BRCA.csv", "TCGA-COAD.csv", "TCGA-ESCA.csv", "TCGA-HNSC.csv", "TCGA-KICH.csv", "TCGA-KIRC.csv", "TCGA-KIRP.csv", "TCGA-LIHC.csv", "TCGA-LUAD.csv", "TCGA-LUSC.csv", "TCGA-PRAD.csv", "TCGA-STAD.csv", "TCGA-THCA.csv")

# proj <- projects[1]
# i <- 1
PathwayMeans <- as.data.frame(matrix(NA, 14, 2))
colnames(PathwayMeans) <- c("Purine", "Pyrimidine")
counter <- 1

for(proj in projects){
  
  dat <- data.table::fread(proj, data.table = FALSE)
  names <- c("ENSEMBL", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
  colnames(dat) <- names
  dat$ENSEMBL <- substr(dat$ENSEMBL, 1, 15)
  df <- dat[dat$ENSEMBL %in% eg$ENSEMBL,]
  PurineMeanFC <- mean(2^df$log2FoldChange[1:3])
  PyrimidineMeanFC <- mean(2^df$log2FoldChange[4:6])
  
  
  PathwayMeans$Purine[counter] <- PurineMeanFC
  PathwayMeans$Pyrimidine[counter] <- PyrimidineMeanFC
  
  counter <- (counter+1)
  
}

PathwayMeans$CT <- substring(projects, 6, 9)

# plotting code

Figure_1B <- ggplot2::ggplot(PathwayMeans, aes(Purine, Pyrimidine, label = CT)) +
  geom_point(color ="black") + geom_text_repel(size =3) + geom_abline(intercept = 0, slope = 1) + theme_classic() + xlim(0, 6) + ylim(0, 6) 







# expression code

# plotting code


#Figure 1C

# purine de novo vs. purine salvage synthesis. 
X <- c("GART", "PAICS", "PFAS")
Y <- c("ATIC", "DGUOK", "HPRT1")

X <- as.data.frame(X)
colnames(X) <- "Genes"
Y <- as.data.frame(Y)
colnames(Y) <- "Genes"
gene.list <- rbind(X,Y)

eg = bitr(gene.list$Genes, fromType="SYMBOL", toType="ENSEMBL", OrgDb="org.Hs.eg.db")

eg$Pathway <- "De novo"
eg$Pathway[5:7] <- "Salvage"

# expression code

# i need the expression FC for each gene for each cancer type. This can be found in the DEA folder
projects <- c("TCGA-BLCA.csv", "TCGA-BRCA.csv", "TCGA-COAD.csv", "TCGA-ESCA.csv", "TCGA-HNSC.csv", "TCGA-KICH.csv", "TCGA-KIRC.csv", "TCGA-KIRP.csv", "TCGA-LIHC.csv", "TCGA-LUAD.csv", "TCGA-LUSC.csv", "TCGA-PRAD.csv", "TCGA-STAD.csv", "TCGA-THCA.csv")

# proj <- projects[1]
# i <- 1
PathwayMeans <- as.data.frame(matrix(NA, 14, 2))
colnames(PathwayMeans) <- c("De novo", "Salvage")
counter <- 1

for(proj in projects){
  
  dat <- data.table::fread(proj, data.table = FALSE)
  names <- c("ENSEMBL", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
  colnames(dat) <- names
  dat$ENSEMBL <- substr(dat$ENSEMBL, 1, 15)
  df <- dat[dat$ENSEMBL %in% eg$ENSEMBL,]
  DeNovoMeanFC <- mean(2^df$log2FoldChange[1:3])
  SalvageMeanFC <- mean(2^df$log2FoldChange[4:6])
  
  
  PathwayMeans$`De novo`[counter] <- DeNovoMeanFC
  PathwayMeans$Salvage[counter] <- SalvageMeanFC
  
  counter <- (counter+1)
  
}

PathwayMeans$CT <- substring(projects, 6, 9)

# plotting code

Figure_1C <- ggplot2::ggplot(PathwayMeans, aes(`De novo`, Salvage, label = CT)) +
  geom_point(color ="black") + geom_text_repel(size =3) + geom_abline(intercept = 0, slope = 1) + theme_classic() + xlim(0, 3.5) + ylim(0, 6) 




#Figure 1D

# purine vs. pyrimidine salvage synthesis. 
X <- c("CAD", "CTPS1", "DTYMK")
Y <- c("CDA", "DCK", "TK1", "TK2")


X <- as.data.frame(X)
colnames(X) <- "Genes"
Y <- as.data.frame(Y)
colnames(Y) <- "Genes"
gene.list <- rbind(X,Y)

eg = bitr(gene.list$Genes, fromType="SYMBOL", toType="ENSEMBL", OrgDb="org.Hs.eg.db")

eg$Pathway <- "De novo"
eg$Pathway[4:7] <- "Salvage"

# expression code

# i need the expression FC for each gene for each cancer type. This can be found in the DEA folder
projects <- c("TCGA-BLCA.csv", "TCGA-BRCA.csv", "TCGA-COAD.csv", "TCGA-ESCA.csv", "TCGA-HNSC.csv", "TCGA-KICH.csv", "TCGA-KIRC.csv", "TCGA-KIRP.csv", "TCGA-LIHC.csv", "TCGA-LUAD.csv", "TCGA-LUSC.csv", "TCGA-PRAD.csv", "TCGA-STAD.csv", "TCGA-THCA.csv")

# proj <- projects[1]
# i <- 1
PathwayMeans <- as.data.frame(matrix(NA, 14, 2))
colnames(PathwayMeans) <- c("De novo", "Salvage")
counter <- 1

for(proj in projects){
  
  dat <- data.table::fread(proj, data.table = FALSE)
  names <- c("ENSEMBL", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
  colnames(dat) <- names
  dat$ENSEMBL <- substr(dat$ENSEMBL, 1, 15)
  df <- dat[dat$ENSEMBL %in% eg$ENSEMBL,]
  DeNovoMeanFC <- mean(2^df$log2FoldChange[1:3])
  SalvageMeanFC <- mean(2^df$log2FoldChange[4:6])
  
  
  PathwayMeans$`De novo`[counter] <- DeNovoMeanFC
  PathwayMeans$Salvage[counter] <- SalvageMeanFC
  
  counter <- (counter+1)
  
}

PathwayMeans$CT <- substring(projects, 6, 9)

# plotting code

Figure_1D <- ggplot2::ggplot(PathwayMeans, aes(`De novo`, Salvage, label = CT)) +
  geom_point(color ="black") + geom_text_repel() + geom_abline(intercept = 0, slope = 1) + theme_classic() + xlim(0, 6) + ylim(0, 6)




