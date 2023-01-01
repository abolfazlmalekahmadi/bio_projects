setwd("~/Desktop/bioprojects/")

library(Biobase)
library(GEOquery)
library(limma)
library(pheatmap)
library(gplots)
library(ggplot2)
library(reshape2)
library(plyr)
library(stringr)
library(Rtsne)


seriesName <- "GSE48558"
platformName <- "GPL6244"
#Download data in txt format :
tgset <- getGEO(seriesName, GSEMatrix =TRUE, AnnotGPL=TRUE , destdir = "Data")


if (length(tgset) > 1) idx <- grep("GPL6244", attr(tgset, "names")) else idx <- 1
ttgset <- tgset[[idx]]  
#Select samples with "AML Patient" Source Name or "Normal" Phenotype :
gset<- ttgset[,which(ttgset$source_name_ch1 == "AML Patient" | ttgset$`phenotype:ch1` == "Normal")]    

ex <- exprs(gset)
gr <- c(rep('AML_Patient',13), rep('Granulocytes', 2), 'B_Cells', 'T_Cells', rep('Granulocytes', 2),
        rep('Monocytes', 2), 'B_Cells', rep('T_Cells',5), 'B_Cells', 'T_Cells', 'B_Cells', 'T_Cells',
        rep('CD34',3), rep('Granulocytes', 7), rep('AML_Patient',2), 'T_Cells', rep('AML_Patient',3),
        rep('B_Cells',7), 'T_Cells', rep('Monocytes', 4), 'Granulocytes', rep('T_Cells',7))



#Differentiation Analysis
gr2 <- gr
gr2 <- factor(gr2)
gset$description <- gr2
#one-hot encoding matrix
onehot <- model.matrix(~ description + 0, gset)
colnames(onehot) <- levels(gr2)


#linear model using limma, then we use that linear model to fit a bayesian model. This bayesian model assumes that about 1% of genes are significant
fit <- lmFit(gset , onehot)
contrast <-
  makeContrasts(AML_Patient - CD34 , levels = onehot)
fit2 <- contrasts.fit(fit , contrast)
fit2 <- eBayes(fit2 , 0.01)
tT <- topTable(fit2,
               adjust = "fdr",
               sort.by = "B" , number = Inf)
# select important columns
tT2 <- subset(tT , select = c("Gene.symbol" , "Gene.title", "adj.P.Val"  , "logFC"))
aml.up <- subset(tT2, logFC > 1 & adj.P.Val < 0.05)
aml.up.genes <- unique(tT2$Gene.symbol)
aml.up.genes <-unique( as.character(strsplit2( (aml.up$Gene.symbol),"///")))
write.table(aml.up.genes, file="Results/AMLUPGENES.txt", quote=F, row.names = F, col.names = F)


tT2 <- subset(tT , select = c("Gene.symbol" , "Gene.title", "adj.P.Val"  , "logFC"))
aml.down <- subset(tT2, logFC < -1 & adj.P.Val < 0.05)
aml.down.genes <- unique(aml.down$Gene.symbol)
aml.down.genes <- unique(as.character(strsplit2( (aml.down$Gene.symbol),"///")))
write.table(aml.down.genes, file="Results/AMLDOWNGENES.txt", quote=F, row.names = F, col.names = F)

