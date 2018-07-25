# PCA Analysis: a Quick Tutorial
# David Chen
# July 2018

# Data set: 21 CCLE breast cancer cell lines 
# Input variables: log2 RPKM of expression of most variable 217 genes
# Outcome variable is BRCAness status predicted by a support vector machine (SVM) copy number classifier

## Load data:
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
topKsites <- read.csv("data/example_cell_line_gene_expression.csv", header=TRUE, row.names=1)
topKsites <- as.data.frame(t(topKsites)) #dimension required for applying prcomp function

## Pull out the numeric expression columns:
mat <- topKsites[ , ! colnames(topKsites) %in% c("Sample","BRCAness")]
mat <- as.matrix(mat)
mode(mat) <- "numeric"

## Calculate first 2 principal components: 
topKsites.pca <- prcomp(mat)
topKsites$Sample <- rownames(topKsites) #for merging later

## Set aside ALL principal components extracted:
pcs <- topKsites.pca$x

## Merge in cell line information:
pcs <- merge(pcs, topKsites[,colnames(topKsites) %in% c("Sample","BRCAness")], by="row.names")

## Now you are ready to visualize:
library(ggplot2)
ggplot(pcs, aes(x=PC1, y=PC2, color=BRCAness)) +
  geom_point(size=3) +
  theme_classic() +
  theme(axis.title=element_text(size=20), axis.text=element_text(size=20))

## Here is convenient ggplot2-based package for PCA: 
# library(devtools)
# install_github("vqv/ggbiplot")
library(ggbiplot)
ggbiplot(topKsites.pca, var.scale=1, varname.size=0, var.axes=FALSE,
         obs.scale=1, groups=topKsites$BRCAness, ellipse=TRUE, circle=TRUE) +
  theme_classic() +
  theme(legend.direction="horizontal", legend.position="top")
