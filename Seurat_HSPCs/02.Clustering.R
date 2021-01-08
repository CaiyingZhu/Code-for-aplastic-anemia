rm(list=ls())
pwd <- getwd()
library(reticulate)
library(Seurat)
library(ggplot2)
library(cowplot)

## load data
load(paste0(pwd,"/int/01.seurat.data.Ctrl.Rdata"))
load(paste0(pwd,"/int/01.seurat.data.AA.Rdata"))

## Normalization and scale data
ctrl.clean.obj <- NormalizeData(object = ctrl.clean.obj, normalization.method = "LogNormalize", scale.factor = 10000)
ctrl.clean.obj <- ScaleData(ctrl.clean.obj, display.progress = F, vars.to.regress = c("percent.mito","nUMI"))

AA.clean.obj <- NormalizeData(object = AA.clean.obj, normalization.method = "LogNormalize", scale.factor = 10000)
AA.clean.obj <- ScaleData(AA.clean.obj, display.progress = F, vars.to.regress= c("percent.mito","nUMI"))

#### Perform a canonical correlation analysis (CCA)
pcaGenes <- function(objs){
	 unique(unlist(lapply(objs, function(obj){
		this.pca <- prcomp(t(as.matrix(obj@data)))
		pc.dim <- 1:3
		pc.gene.num <- 200
		this.pca.genes <- unique(unlist(lapply(pc.dim, function(x){
			head(row.names(this.pca$rotation[order(abs(this.pca$rotation[,x]), decreasing=T),]), n=pc.gene.num)
		})))
		return(this.pca.genes)
	     })))
     }

## PCA genes
pc.genes <- pcaGenes(c(ctrl.clean.obj, AA.clean.obj))
length(pc.genes) # 644
save(pc.genes,file = paste0(pwd,"/int/pc.gene.644.Rdata"))

# group label
ctrl.clean.obj@meta.data$source <- "Ctrl"
AA.clean.obj@meta.data$source <- "AA"

# combine
hspc.obj <- RunCCA(ctrl.clean.obj, AA.clean.obj, genes.use = pc.genes, num.cc = 30, group.by = "source")
dim(hspc.obj@raw.data) # 16571  4881
dim(hspc.obj@data)     # 16571  4881

# visualize results of CCA plot CC1 versus CC2 and look at a violin plot
p1 <- DimPlot(object = hspc.obj, reduction.use = "cca", group.by = "source", pt.size = 0.5, do.return = TRUE)
p2 <- VlnPlot(object = hspc.obj, features.plot = "CC1", group.by = "source", do.return = TRUE)
plot_grid(p1, p2)
ggsave(paste(pwd, "/output/merge/combine.CC1_and_CC2.pdf", sep=""), height=5, width=12)

# measure of correlation strength for each CC
p3 <- MetageneBicorPlot(hspc.obj, grouping.var = "source", dims.eval = 1:30, display.progress = FALSE)
plot(p3)
ggsave(paste(pwd, "/output/merge/combine.CC_CorStrength.pdf", sep=""), height=5, width=7)

# heatmap the 30 CC
pdf(paste(pwd, "/output/merge/combine.CC_heatmap.pdf", sep=""), height=12, width=7)
DimHeatmap(object = hspc.obj, reduction.type = "cca", cells.use = 500, dim.use = 1:30, do.balanced = TRUE)
dev.off()

# choose CC23 for dimension cutoff
# Align CCA subspaces
dimSelected <- 23
hspc.aligned.obj <- AlignSubspace(hspc.obj, reduction.type = "cca", grouping.var = "source", dims.align = 1:dimSelected)
# visualize aligned CCA
p1 <- VlnPlot(object = hspc.aligned.obj, features.plot = "ACC1", group.by = "source", do.return = TRUE)
p2 <- VlnPlot(object = hspc.aligned.obj, features.plot = "ACC2", group.by = "source", do.return = TRUE)
pdf(paste(pwd, "/output/merge/aligned.CCA1_and_CCA2.pdf", sep=""), height=5, width=7)
plot_grid(p1, p2)
dev.off()

#### Dimension reduction and clustering
hspc.aligned.obj <- RunTSNE(hspc.aligned.obj, reduction.use = "cca.aligned", dims.use = 1:dimSelected, do.fast = T, dim.embed = 2)
hspc.aligned.obj <- FindClusters(hspc.aligned.obj, reduction.type = "cca.aligned",resolution = c(0.6,0.8,1,1.2), dims.use = 1:dimSelected, force.recalc =T)
# Visualization
p1 <- TSNEPlot(hspc.aligned.obj, do.return = T, pt.size = 0.5, group.by = "source") + xlab("source")
p2 <- TSNEPlot(hspc.aligned.obj, do.return = T, pt.size = 0.5, group.by = "SampleID") + xlab("SampleID")
p3 <- TSNEPlot(hspc.aligned.obj, do.label = T, do.return = T, pt.size = 0.5, group.by ="res.0.6") + xlab("res_0.6")
p4 <- TSNEPlot(hspc.aligned.obj, do.label = T, do.return = T, pt.size = 0.5, group.by ="res.0.8") + xlab("res_0.8")
p5 <- TSNEPlot(hspc.aligned.obj, do.label = T, do.return = T, pt.size = 0.5, group.by ="res.1") + xlab("res_1")
p6 <- TSNEPlot(hspc.aligned.obj, do.label = T, do.return = T, pt.size = 0.5, group.by ="res.1.2") + xlab("res_1.2")
plot_grid(p1, p2, p3, p4, p5, p6, nrow=2)
ggsave(paste(pwd, "/output/aligned.tSNE_clustering.dimUse", dimSelected, ".pdf",sep=""), height=8, width=14, device="pdf")

# use resolution 0.8 to identify clusters
hspc.aligned.obj <- FindClusters(hspc.aligned.obj, reduction.type = "cca.aligned", resolution = 0.8, dims.use = 1:dimSelected, force.recalc =T)

## Conserved cell type markers
# - Find markers in cluster in Ctrl cells, denoted as A
# - Find markers in cluster in AA cells, denoted as B
# - Conserved markers: intersect(A, B)

conserved.markers <- c()
for(i in levels(hspc.aligned.obj@ident)){
	 print(paste0("Processing cluster c", i, " ..."))
	 this.marker <- FindConservedMarkers(
		hspc.aligned.obj,
		ident.1 = i,
		grouping.var = "source",
		print.bar = FALSE, only.pos = TRUE
	 )
	 this.marker <- this.marker[order(this.marker[,2] + this.marker[,7], decreasing =T), ]
	 this.marker$gene <- row.names(this.marker)
	 this.marker$cluster <- i
	 conserved.markers <- rbind(conserved.markers, this.marker)
     }

# stat numbers
table(conserved.markers$cluster)
#  0   1   2   3   4   5   6   7   8
# 618 508 274 183 642  45 481 411 486

num.topMarkers <- 6
for(i in levels(hspc.aligned.obj@ident)){
	FeaturePlot(object = hspc.aligned.obj, features.plot = head(conserved.markers$gene[conserved.markers$cluster ==i], n=num.topMarkers), min.cutoff = "q8", cols.use = c("lightgrey","blue"), pt.size = 0.5)
	ggsave(paste(pwd, "/output/hspc.topFeaturePlot.c", i, ".ctrl1-4.AA.pdf", sep=""), height=7, width=7, device="pdf")
}

## rename 
oldCelltype <- c("0","1","2","3","4","5","6","7","8")
newCelltype <- c("Neu2","MLP","HSC/MPP","LMPP","MEP","Neu1","MD1","MD2","EBM")

for(i in 0:(length(oldCelltype)-1)){
	 hspc.aligned.obj <- RenameIdent(object = hspc.aligned.obj, old.ident.name = i, new.ident.name = newCelltype[i+1])
     }
	 
new.ident.order <- c("HSC/MPP","LMPP","MEP","MLP","EBM","Neu1","Neu2","MD1","MD2")
hspc.aligned.obj@ident <- factor(hspc.aligned.obj@ident,levels=new.ident.order)

hspc.aligned.obj@meta.data <- hspc.aligned.obj@meta.data[,1:16]
hspc.aligned.obj@meta.data$cellGroup <- hspc.aligned.obj@ident

## find all markers
hspc.allMarkers <- FindAllMarkers(hspc.aligned.obj, only.pos = TRUE)
table(hspc.allMarkers$cluster)
## HSC/MPP    LMPP     MEP     MLP     EBM    Neu1    Neu2     MD1     MD2
##    430     300     872     543     799     123     877     673     636

## Dendrogram of clusters
data.byCluster <- matrix(rep(1, length(pc.genes)*length(levels(hspc.aligned.obj@ident))), ncol=length(levels(hspc.aligned.obj@ident)))
colnames(data.byCluster) <- new.ident.order
for(i in new.ident.order){
	 data.byCluster[,i] <- Matrix::rowMeans(hspc.aligned.obj@data[pc.genes,hspc.aligned.obj@ident == i])
     }	 
group.clust <- hclust(as.dist(1-cor(data.byCluster,method="spearman")), method = "ward.D2")
library(dendextend)
pdf(paste0(pwd, "/output/hspc.dendrogram_Clusters.pdf"),width=6,height=3)
as.dendrogram(group.clust) %>% 
	plot(horiz=FALSE, axes=FALSE, main="Clustering of cell groups")
dev.off()
	
## save data	 
save(hspc.aligned.obj,pc.genes,new.ident.order, file = paste0(pwd,"/int/02.Seurat.AACtrl.aligned.Rdata"))

