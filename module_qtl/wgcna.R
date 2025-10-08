library(WGCNA)
library(spqn)
library(sva)
library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
pcgs <- getBM(attributes = c("ensembl_gene_id_version", "hgnc_symbol", "transcript_biotype", "chromosome_name"), 
          filters = c("transcript_biotype", "chromosome_name"), 
          values = list(c("protein_coding"), c(1:22)), mart = ensembl)
allowWGCNAThreads(20)

# load data
lung_logtpms <- read.table("../lung_tpms.tsv", header = TRUE, row.names = 1, sep = "\t")
tum_logtpms <- read.table("../tumor_tpms.tsv", header = TRUE, row.names = 1, sep = "\t")

# genes used in cis-eqtls
expressed_genes_lung <- tail(colnames(read.csv("../qtltools/input_data/Expression_lung_tpm01_20perc_reads6_20perc.csv", head=TRUE, nrows=1)), -1)
expressed_genes_tumor <- tail(colnames(read.csv("../qtltools/input_data/Expression_tumor_tpm01_20perc_reads6_20perc.csv", head=TRUE, nrows=1)), -1)

# keep expressed pcgs
lung_logtpms <- lung_logtpms[,names(lung_logtpms)[(names(lung_logtpms) %in% pcgs$ensembl_gene_id_version) & (names(lung_logtpms) %in% expressed_genes_lung)]]
tum_logtpms <- tum_logtpms[,names(tum_logtpms)[(names(tum_logtpms) %in% pcgs$ensembl_gene_id_version) & (names(tum_logtpms) %in% expressed_genes_tumor)]]

# log-transform expression
lung_logtpms <- log(lung_logtpms + 0.5)
tum_logtpms <- log(tum_logtpms + 0.5)

# estimate number of pcs to remove and regress out pcs
num_pcs_lung <- sva::num.sv(t(lung_logtpms), model.matrix(~ 1, data = lung_logtpms), method="be")
num_pcs_tum <- sva::num.sv(t(tum_logtpms), model.matrix(~ 1, data = tum_logtpms), method="be")

lung_corr_logtpms <- sva_network(lung_logtpms, num_pcs_lung)
tum_corr_logtpms <- sva_network(tum_logtpms, num_pcs_tum)

# correlation matrix of gene expression
lung_cor <- bicor(lung_corr_logtpms)
tum_cor <- bicor(tum_corr_logtpms)
write.table(rownames(lung_cor), "out_wgcna_v6/genes_names_in_order_lung.txt")
write.table(rownames(tum_cor), "out_wgcna_v6/genes_names_in_order_tumor.txt")

# normalize correlations
lung_avg <- apply(lung_logtpms, 2, mean)
tum_avg <- apply(tum_logtpms, 2, mean)

pdf("out_wgcna_v6/cor_before_norm_lung.pdf")
plot_signal_condition_exp(lung_cor, lung_avg, signal=0.001)
dev.off()

pdf("out_wgcna_v6/cor_before_norm_tum.pdf")
plot_signal_condition_exp(tum_cor, tum_avg, signal=0.001)
dev.off()

lung_cor_norm <- normalize_correlation(lung_cor, lung_avg, 10, 1000, 9)
tum_cor_norm <- normalize_correlation(tum_cor, tum_avg, 10, 1000, 9)

pdf("out_wgcna_v6/cor_after_norm_lung.pdf")
plot_signal_condition_exp(lung_cor_norm, lung_avg, signal=0.001)
dev.off()

pdf("out_wgcna_v6/cor_after_norm_tum.pdf")
plot_signal_condition_exp(tum_cor_norm, tum_avg, signal=0.001)
dev.off()

### WGCNA

# 1. select power

sft_lung <- pickSoftThreshold(lung_cor_norm, dataIsExpr = FALSE, networkType = "unsigned", powerVector = seq(from = 1, to = 41, by = 2), verbose = 0)

powers <- seq(from = 1, to = 41, by = 2)
pdf("out_wgcna_v6/pickSoftTreshold_lung.pdf")
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft_lung$fitIndices[,1], -sign(sft_lung$fitIndices[,3])*sft_lung$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft_lung$fitIndices[,1], -sign(sft_lung$fitIndices[,3])*sft_lung$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red") 
# Mean connectivity as a function of the soft-thresholding power
plot(sft_lung$fitIndices[,1], sft_lung$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) 
text(sft_lung$fitIndices[,1], sft_lung$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

sft_tum <- pickSoftThreshold(tum_cor_norm, dataIsExpr = FALSE, networkType = "unsigned", powerVector = seq(from = 2, to = 20, by = 2), verbose = 0)

powers <- seq(from = 2, to = 20, by = 2)
pdf("out_wgcna_v6/pickSoftTreshold_tum.pdf")
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft_tum$fitIndices[,1], -sign(sft_tum$fitIndices[,3])*sft_tum$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft_tum$fitIndices[,1], -sign(sft_tum$fitIndices[,3])*sft_tum$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red") 
# Mean connectivity as a function of the soft-thresholding power
plot(sft_tum$fitIndices[,1], sft_tum$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) 
text(sft_tum$fitIndices[,1], sft_tum$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

# 2. make weighted adjacency matrix

adj_lung <- adjacency.fromSimilarity(lung_cor_norm, power=sft_lung$powerEstimate, type = "unsigned")
adj_tum <- adjacency.fromSimilarity(tum_cor_norm, power=sft_tum$powerEstimate, type = "unsigned")

# 3. make weighted graph TOM similarity matrix

tom_lung <- TOMsimilarity(adj_lung)
tom_tum <- TOMsimilarity(adj_tum)
row.names(tom_lung) <- row.names(lung_cor)
row.names(tom_tum) <- row.names(tum_cor)

saveRDS(tom_lung, "out_wgcna_v6/TOM_network_lung.rds")
saveRDS(tom_tum, "out_wgcna_v6/TOM_network_tum.rds")

# 4. hierachical clustering of topological similarity

dist_lung <- as.dist(1 - tom_lung)
dist_tum <- as.dist(1 - tom_tum)

dendro_lung <- hclust(dist_lung)
dendro_tum <- hclust(dist_tum)

cutd_lung <- cutreeDynamic(dendro_lung, distM = as.matrix(dist_lung), minClusterSize = 10)
cutd_tum <- cutreeDynamic(dendro_tum, distM = as.matrix(dist_tum), minClusterSize = 10)

modules_lung <- mergeCloseModules(lung_corr_logtpms, cutd_lung, cutHeight = 0.1)
modules_tumor <- mergeCloseModules(tum_corr_logtpms, cutd_tum, cutHeight = 0.1)

saveRDS(modules_lung, "out_wgcna_v6/merged_modules_lung.rds")
saveRDS(modules_tumor, "out_wgcna_v6/merged_modules_tum.rds")
write.table(data.frame(gene=names(lung_corr_logtpms), mods=modules_lung$colors), "out_wgcna_v6/modules_lung.tsv", sep="\t", quote=F, row.names=F)
write.table(data.frame(gene=names(tum_corr_logtpms), mods=modules_tumor$colors), "out_wgcna_v6/modules_tumor.tsv", sep="\t", quote=F, row.names=F)
