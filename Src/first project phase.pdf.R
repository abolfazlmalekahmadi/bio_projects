setwd("~/Desktop/bioprojects/")
library(GEOquery)
library(limma)
library(umap)
library(pheatmap)
library(Rtsne)
library(plyr)
library(dplyr)
library(magrittr)
library(ggpubr)



series <- "GSE48558"
gset <- getGEO(series, GSEMatrix =TRUE, AnnotGPL=TRUE , destdir = "Data/")
if (length(gset) > 1) idx <- grep("GPL6244", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
fvarLabels(gset) <- make.names(fvarLabels(gset))
gsms <- paste0("AAAAAAAAAAAAAXXXXXXXXXXXXXXXXXXXXXXXXXXXGXXXGXXXXX",
               "XXXXXXXXXXXXXXXXXXBXTXXXGXGMMBXTXXTTXXTTXBXTXBXTXC",
               "XXXCXXXCXXXXXXXXXXXXXXXXXXXXXXXXXXXXXGGGGGGGAATAAA",
               "BBBBBBBTMMMMGTTTTTTT")

sml <- strsplit(gsms, split="")[[1]]
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]


# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

ex <- exprs(gset)
pdf("Results/Boxplot.pdf" , width  = 64)
boxplot(ex)
dev.off()
#temp <- normalizeQuantiles(ex)
#pdf("Results/Boxplot_norm.pdf" , width = 64)
#boxplot(temp)
#dev.off()
pdf("Results/pheatmap.pdf" , width = 15 , height = 15)
pheatmap(cor(ex), labels_row = sml , labels_col = sml ,)
dev.off()
####
pc <-prcomp(ex)
pdf("Results/Pc.pdf")
plot(pc$x[, 1:2])
plot(pc)
dev.off()
ex.scale <- t(scale(t(ex),scale = F))
pc <-prcomp(ex.scale)
pdf("Results/Pc_scaled.pdf")
plot(pc)
plot(pc$x[, 1:2])
dev.off()

pcr <- data.frame(pc$rotation[,1:3] , Group = sml)
pdf("Results/PCA_samples.pdf")
ggplot(pcr,aes(PC2 ,PC3 , color = Group))+geom_point(size=3)+ theme_bw()
dev.off()
ex_transposed = t(ex)
pers = 3:20
tsne <- Rtsne(ex_transposed, perplexity = 6, check_duplicates = FALSE)
par(mfrow=c(1,2))
labels <- unique(sml)
colors <- rainbow(length(unique(sml)))
bg_colors <- mapvalues(sml, labels, colors)
pdf("Results/tsne_result.pdf")
plot(tsne$Y, col = "blue", pch = 19, cex = 1.5)
plot(tsne$Y, col = "black", bg= bg_colors, pch = 21, cex = 1.5)
legend(x = "topleft", legend = labels, fill = colors)
dev.off()

par(mfrow = c(1,1))

# Cmpute MDS
mds <- ex_transposed %>%
  dist() %>%          
  cmdscale() %>%
  as_tibble()
colnames(mds) <- c("Dim.1", "Dim.2")
# Plot MDS
plot(mds, col = "black", bg= bg_colors, pch = 21, cex = 1.5)
legend(x="top",legend = labels, fill = colors)



