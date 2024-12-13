library(WGCNA)
library(tidyverse)
library(DESeq2)

d = readRDS('../../data/processed_datasets/master_allData_batchCorrected.RDS')
m = d$subject_specimen


counts = as.data.frame(d$pbmc_gene_expression$raw_count$raw_data)
counts[1:5, 1:5]
counts = counts[, which(colnames(counts) %in% m$specimen_id)]

m = m %>%
    filter(specimen_id %in% colnames(counts))

# vst transform
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = m ,
                              design = ~ as.factor(infancy_vac))
dds
dds = estimateSizeFactors(dds)
vsd <- vst(dds, blind=TRUE)
dim(assay(vsd))

# HVF
topVarGenes <- head(order(rowVars(assay(vsd)),decreasing=TRUE),5000)
topVarGenes <- rownames(vsd)[topVarGenes]

dds <- dds[topVarGenes,]
dds
dds = estimateSizeFactors(dds)
vsd <- vst(dds, blind=TRUE)
dim(assay(vsd))

bef_plusHVF=plotPCA(vsd, "dataset") # PCA before

assay(vsd) <- limma::removeBatchEffect(assay(vsd), batch=as.vector(dds$dataset))


vsd$subject_id = as.factor(vsd$subject_id)

bef_noHVF1 = plotPCA(vsd, "dataset") # PCA before
bef_noHVF2 = plotPCA(vsd, "subject_id") # PCA before

assay(vsd) <- limma::removeBatchEffect(assay(vsd), batch=as.vector(dds$dataset))
aft_noHVF1 = plotPCA(vsd, "dataset")
aft_noHVF2 = plotPCA(vsd, "subject_id")

aft_noHVF1
#### still some pretty bad batch effects - will have to try and remove these

# filter for inout into WGCNA
ds=c('2020_dataset', '2021_dataset')

# meta-data, keep only 2020 and 2021
m = m %>%
  filter(dataset %in% ds, timepoint %in% c(0,1))


# transpose for input into WGCNA
expr = data.frame(t(assays(vsd)[[1]]))
expr = expr %>%
  filter(rownames(expr) %in% m$specimen_id)
expr[1:5, 1:5]

dim(expr)

# Group data in a dendogram to check outliers
sampleTree = hclust(dist(expr), method = "average")

par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

#Plot a line showing the cut-off
abline(h = 65, col = "red") 
abline(h = 63, col = "red") 
abline(h =50, col = "red") 

#Determine clusters below the line
#help("cutreeStatic")
clust = cutreeStatic(sampleTree, cutHeight = 50, minSize = 10)
table(clust)
clust

#Cluster 1 contains the samples we want to keep.
keepSamples = (clust==1)
expr = expr[keepSamples, ]
#dim(expression0)
nGenes = ncol(expr)
nSamples = nrow(expr)

#Regrouping samples
sampleTree2 = hclust(dist(expr), method = "average")

m = m %>%
    filter(specimen_id %in% rownames(expr))

m$infancy_vac_num = 1
m[m$infancy_vac == 'wP', 'infancy_vac_num'] = 2
table(m$infancy_vac_num)

m$dataset_num = 2020
m[m$dataset == "2021_dataset", 'dataset_num'] = 2021
table(m$dataset_num)

traitColors = numbers2colors(m[,c('timepoint', 'infancy_vac_num', 'subject_id', 'dataset_num')], signed = FALSE)
colnames(traitColors) = c('timepoint', 'infancy_vac_num', 'subject_id', 'dataset_num')

#Plot a sample dendogram with the colors below
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels =colnames(traitColors), 
                    main = "Sample dendrogram and trait heatmap")


## CREATE NETWORK

# Choose a set of soft threshold parameters
powers = c(c(1:20), seq(from = 22, to=30, by=2))

sft = pickSoftThreshold(expr, powerVector = powers, verbose = 5) 

par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red") 
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) 
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# Turn data expression into topological overlap matrix
power=sft$powerEstimate
power

adjacency = adjacency(expr, power = power, type='signed')
TOM = TOMsimilarity(adjacency); # Turn adjacency into topological overlap
dissTOM = 1-TOM
dim(dissTOM)


## CONSTRUCT MODELS

# Plot gene tree
geneTree = hclust(as.dist(dissTOM), method = "average");
#pdf(file = "3-gene_cluster.pdf", width = 12, height = 9);
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)

# We like large modules, so we set the minimum module size relatively high:
# minModuleSize = 30;
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,deepSplit = 2, 
                            pamRespectsDendro = FALSE,minClusterSize = 30)
table(dynamicMods)
length(table(dynamicMods)) 
# Convert numeric labels into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath

plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",dendroLabels = FALSE,
                    hang = 0.03,addGuide = TRUE, guideHang = 0.05,main = "Gene dendrogram and module colors")

## MERGE MODULES

# Calculate eigengenes
MEList = moduleEigengenes(expr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

# Merge close modules
MEDissThres=c(0.6, 0.4, 0.3)
abline(h=MEDissThres, col = "red")

merge = mergeCloseModules(expr, dynamicColors, cutHeight = MEDissThres[[1]], verbose = 3) 
mergedColors = merge$colors
unique(mergedColors)

merge2 = mergeCloseModules(expr, dynamicColors, cutHeight = MEDissThres[[2]], verbose = 3) 
mergedColors2 = merge2$colors
unique(mergedColors2)

merge3 = mergeCloseModules(expr, dynamicColors, cutHeight = MEDissThres[[3]], verbose = 3) 
mergedColors3 = merge3$colors
unique(mergedColors3)


# Plot merged module tree
#pdf(file = "5-merged_Module_Tree.pdf", width = 12, height = 9)  
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors, mergedColors2, mergedColors3), main='2020, 2021 : TP 0 and 1, HVG 5000',
                    c("Dynamic Tree Cut", "Merged dynamic 0.6", "Merged dynamic 0.4", "Merged dynamic 0.3"), dendroLabels = FALSE, 
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05)  



mergedMEs = merge3$newMEs  
dim(merge3$newMEs)

library(pheatmap)
pheatmap(merge$oldMEs,cluster_col=T,cluster_row=T,show_rownames=F,show_colnames=T,fontsize=6)
pheatmap(merge$newMEs,cluster_col=T,cluster_row=T,show_rownames=F,show_colnames=T,fontsize=6)

# export gene lists
mod = list()
for (i in 1:length(merge3$newMEs)){
    modules = c(substring(names(merge3$newMEs)[i], 3));
    genes = colnames(expr)
    inModule = is.finite(match(dynamicColors,modules))
    modGenes = genes[inModule]
    mod[[i]] = modGenes
}
names(mod) = colnames(merge3$newMEs)

lapply(mod, length)


# save


MEList2 =  moduleEigengenes(expr, colors = mergedColors2)
MEList2$eigengenes[1:5, 1:5]
# export gene lists
mod = list()
for (i in 1:length(MEList2$eigengenes)){
    modules = c(substring(names(MEList2$eigengenes)[i], 3));
    genes = colnames(expr)
    inModule = is.finite(match(MEList2$validColors,modules))
    modGenes = genes[inModule]
    mod[[i]] = modGenes
}
names(mod) = colnames(MEList2$eigengenes)
MEList2$modules = mod

MEList1 =  moduleEigengenes(expr, colors = mergedColors)
MEList1$eigengenes[1:5, 1:5]
# export gene lists
mod = list()
for (i in 1:length(MEList1$eigengenes)){
    modules = c(substring(names(MEList1$eigengenes)[i], 3));
    genes = colnames(expr)
    inModule = is.finite(match(MEList1$validColors,modules))
    modGenes = genes[inModule]
    mod[[i]] = modGenes
}
names(mod) = colnames(MEList1$eigengenes)
MEList1$modules = mod

save(MEList1, MEList2, file = '../data/WGCNA/first_pass.RData')

expr_all = t(assays(vsd)[[1]])
expr_all[1:5, 1:5]
dim(expr_all)
expr_eigengenes_1 = moduleEigengenes(expr = expr_all, colors = MEList1$validColors)
expr_eigengenes_2 = moduleEigengenes(expr = expr_all, colors = MEList2$validColors)

save(expr_eigengenes_1, expr_eigengenes_2, file = '../data/WGCNA/first_pass_eigengenes-expr.RData')


MEList1$eigengenes[1:5, 1:5]
expr_eigengenes_1$eigengenes['77', 'MEivory']
