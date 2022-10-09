# Libraries
library(GEOquery)
library(limma)
library(umap)
library(stringr)
library(Rtsne)
library(pheatmap)

# Load Data
curD <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(sub(paste0("/",sub("(.+)/","",curD)),"",curD))


gset <- getGEO("GSE48558", GSEMatrix =TRUE, AnnotGPL=TRUE, destdir='data/')
if (length(gset) > 1) idx <- grep("GPL6244", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

func <- function(x) {
  if (gset$`phenotype:ch1`[x] == "Normal") {
    return("N")
  } else if (gset$source_name_ch1[x] == "AML Patient") {
    return("T")
  }
  else{
    return("E")
  }
}

group <- sapply(1:length(gset$`phenotype:ch1`) , func)

# Delete Extra Samples
temp <- which(group != 'E')
group <- factor(group[temp])
gset <- gset[, temp]

# Draw Boxplot
ex <- exprs(gset)
ord <- order(group)
colors <- c("#1B9E77", "#7570B3")
png("results/boxplot.png" , width = 4000, height = 700)
boxplot(ex[,ord], main = "Boxplot", boxwex=0.8, notch=T, las=3, col=colors[group[ord]], outline = FALSE, cex=0.5)
legend("topright", legend=levels(group), pch=20, col=colors, title="Group", pt.cex=1.5)
dev.off()

# Dimension Reduction

# PCA
scaled <- t(scale(t(ex) , scale = F))

pca_result <- prcomp(scaled)
png("results/pca_genes.png" , width = 1000, height = 400)
plot(pca_result$x[,1:2], main = "gene PCA (genes)", col=group, pch=20, cex=1.5)
legend("topleft", legend=levels(group), pch=20, col=1:nlevels(group), title="Group", pt.cex=1.5)
dev.off()

pca_result <- prcomp(t(scaled))
png("results/pca_samples.png" , width = 1000, height = 400)
plot(pca_result$x[,1:2], main = "gene PCA (samples)", col=group, pch=20, cex=1.5)
legend("topleft", legend=levels(group), pch=20, col=1:nlevels(group), title="Group", pt.cex=1.5)
dev.off()

# tSNE
tsne_result <- Rtsne(ex, dims = 2, perplexity=30, verbose=TRUE, max_iter = 1000)
png("results/tSNE_genes.png" , width = 1000, height = 400)
plot(tsne_result$Y, main="gene tSNE (genes)", xlab="dim1", ylab="f=dim2", col=group, pch=20, cex=1.5)
legend("topright", legend=levels(group), pch=20, col=1:nlevels(group), title="Group", pt.cex=1.5)
dev.off()

tsne_result <- Rtsne(t(ex), dims = 2, perplexity=20, verbose=TRUE, max_iter = 1000)
png("results/tSNE_samples.png" , width = 1000, height = 400)
plot(tsne_result$Y, main="gene tSNE (samples)", xlab="dim1", ylab="f=dim2",  col=group, pch=20, cex=1.5)
legend("topright", legend=levels(group), pch=20, col=1:nlevels(group), title="Group", pt.cex=1.5)
dev.off()

# UMAP
umap_result <- umap(ex, n_neighbors = 15, random_state = 123)
png("results/umap_genes.png" , width = 1000, height = 400)
plot(umap_result$layout, main="gene umap (genes)", xlab="dim1", ylab="dim2", col=group, pch=20, cex=1.5)
legend("topright", legend=levels(group), pch=20, col=1:nlevels(group), title="Group", pt.cex=1.5)
dev.off()

umap_result <- umap(t(ex), n_neighbors = 15, random_state = 123)
png("results/umap_samples.png" , width = 1000, height = 400)
plot(umap_result$layout, main="gene umap (samples)", xlab="dim1", ylab="dim2", col=group, pch=20, cex=1.5)
legend("topright", legend=levels(group), pch=20, col=1:nlevels(group), title="Group", pt.cex=1.5)
dev.off()

# Correlation
png("results/heatmap_reducted.png" , width = 2000, height = 400)
pheatmap(cor(t(tsne_result$Y)), labels_row = factor(gset$source_name_ch1) , labels_col = factor(gset$source_name_ch1), fontsize_col = 5, fontsize_row = 5)    
dev.off()

png("results/heatmap.png" , width = 2000, height = 400)
pheatmap(cor(ex), labels_row = factor(gset$source_name_ch1) , labels_col = factor(gset$source_name_ch1), fontsize_col = 5, fontsize_row = 5)    
dev.off()

# Differential Expression Analysis
gset$group <- group
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(group)

fit <- lmFit(gset, design)  # fit linear model

cont.matrix <- makeContrasts(contrasts="T-N", levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
tT <- subset(tT, select=c("Gene.symbol", "Gene.ID", "adj.P.Val", "logFC"))

tT.up <- subset(tT, adj.P.Val < 0.05 & logFC > 1)
tT.down <- subset(tT, adj.P.Val < 0.05 & logFC < -1)

up_genes <- unique(as.character(strsplit2((tT.up$Gene.symbol),"///")))
down_genes <- unique(as.character(strsplit2((tT.down$Gene.symbol),"///")))

write.table(up_genes, 'results/upgenes.txt', row.names = F, col.names = F, quote = F)
write.table(down_genes, 'results/downgenes.txt', row.names = F, col.names = F, quote = F)

