rm(list=ls())
library(reticulate)
use_python("/public/home/zhucy/anaconda3/bin/python")
py_config()
py_available(umap)
library(Seurat)
library(ggplot2)
library(cowplot)
library(uwot)
setwd("/public/home/zhucy/SingleCell/AA/Tcell/Seurat")
pwd <- getwd()

#cell.label <- "CD4T"
cell.label <- "CD8T"
## load Ctrl and AA objects
load(paste0(pwd,"/int/01.seurat.data.v2.AA.Rdata"))
load(paste0(pwd,"/int/01.seurat.data.v2.ctrl.Rdata"))
load(paste0(pwd,"/int/02.celltype.",cell.label,".meta.Rdata"))

## Ctrl obj
Ctrl.obj <- NormalizeData(object = Ctrl.obj, normalization.method = "LogNormalize", scale.factor = 10000)
Ctrl.obj <- ScaleData(Ctrl.obj, display.progress = F, vars.to.regress = c("percent.mito","nUMI"))

## AA obj
AA.obj <- NormalizeData(object = AA.obj, normalization.method = "LogNormalize", scale.factor = 10000)
AA.obj <- ScaleData(AA.obj, display.progress = F, vars.to.regress = c("percent.mito","nUMI"))

## Perform a canonical correlation analysis (CCA)
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
pc.genes <- pcaGenes(c(Ctrl.obj, AA.obj))
length(pc.genes)

## combine
combine.obj <- RunCCA(Ctrl.obj, AA.obj, genes.use = pc.genes, num.cc = 30, group.by = "source")
dim(combine.obj@raw.data) # 14641  2057 CD8T: 14641  2024
dim(combine.obj@data)     # 14641  2057 CD8T: 14641  2024

## visualize results of CCA plot CC1 versus CC2 and look at a violin plot
p1 <- DimPlot(object = combine.obj, reduction.use = "cca", group.by = "source", pt.size = 0.5, do.return = TRUE)
p2 <- VlnPlot(object = combine.obj, features.plot = "CC1", group.by = "source", do.return = TRUE)
plot_grid(p1, p2)
ggsave(paste(pwd, "/output/combine.CC1_and_CC2.",cell.label,".pdf", sep=""), height=5, width=12)
dev.off()

## measure of correlation strength for each CC
p3 <- MetageneBicorPlot(combine.obj, grouping.var = "source", dims.eval = 1:30, display.progress = FALSE)
plot(p3)
ggsave(paste(pwd, "/output/combine.CC_CorStrength.",cell.label,".pdf", sep=""), height=5, width=7)
dev.off()

## heatmap the 30 CC
pdf(paste(pwd, "/output/combine.CC_heatmap.",cell.label,".pdf", sep=""), height=12, width=7)
DimHeatmap(object = combine.obj, reduction.type = "cca", cells.use = 500, dim.use = 1:30, do.balanced = TRUE)
dev.off()

#### Dimension reduction 
combine.obj<- FindVariableGenes(object = combine.obj, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.5, x.high.cutoff = 3, y.cutoff = 0.5)
combine.obj <- RunPCA(combine.obj, pcs.compute = 30, pcs.print=NULL)
p1 <- DimPlot(object = combine.obj, reduction.use = "pca", group.by = "source", pt.size = 0.5, do.return = TRUE, dim.1=1, dim.2=2)
p2 <- DimPlot(object = combine.obj, reduction.use = "pca", group.by = "source", pt.size = 0.5, do.return = TRUE, dim.1=2, dim.2=3)
p3 <- VlnPlot(object = combine.obj, features.plot = "PC1", group.by = "source", do.return = TRUE)
p4 <- VlnPlot(object = combine.obj, features.plot = "PC2", group.by = "source", do.return = TRUE)
plot_grid(p1, p2, p3, p4)
ggsave(paste(pwd, "/output/PC1_and_PC2.",cell.label,".pdf", sep=""), height=5, width=8)

## choose CC23 for dimension cutoff
## check for dimension choose
## PCElbowPlot
PCElbowPlot(object = combine.obj)
ggsave(paste0(pwd, "/output/", cell.label, ".PCElbowPlot.pdf"))

## Align CCA subspaces
dimSelected <- 6
aligned.obj <- AlignSubspace(combine.obj, reduction.type = "cca", grouping.var = "source", dims.align = 1:dimSelected)

# visualize aligned CCA
p1 <- VlnPlot(object = aligned.obj, features.plot = "ACC1", group.by = "source", do.return = TRUE)
p2 <- VlnPlot(object = aligned.obj, features.plot = "ACC2", group.by = "source", do.return = TRUE)
pdf(paste(pwd, "/output/aligned.CCA1_and_CCA2.",cell.label,".pdf", sep=""), height=5, width=7)
plot_grid(p1, p2)
dev.off()

## Run uwot and visualization
aligned.obj <- RunUMAP(aligned.obj, reduction.use = "cca.aligned", dims.use = 1:dimSelected, min_dist=0.3)
p1 <- DimPlot(aligned.obj, reduction.use="umap", do.return = T, pt.size = 1, group.by = "cellGroup") + xlab("cellGroup")
p2 <- DimPlot(aligned.obj, reduction.use="umap", do.return = T, pt.size = 1, group.by = "Disease_types") + xlab("Disease_types")
plot_grid(p1, p2)
ggsave(paste(pwd, "/output/aligned.UMAP.dimUse", dimSelected,".",cell.label,".pdf", sep=""), height=4, width=12, device="pdf")

## prepared data for uwot
new.ident.order <- c("Naive","Memory","Effector")		 
data <- as.matrix(aligned.obj@dr$cca@cell.embeddings)

set.seed(2)
umap <- umap(data[,1:dimSelected],n_neighbors = 30, n_components = 2, metric = "euclidean", spread = 0.8, min_dist = 0.2)
umap <- as.data.frame(umap)
colnames(umap) <- c("UMAP1","UMAP2")
rownames(umap) <- rownames(data)
meta <- aligned.obj@meta.data
umap.plot <- cbind(umap,meta)
umap.plot$cellGroup <- factor(umap.plot$cellGroup,levels=new.ident.order)
umap.plot$source <- factor(umap.plot$source,levels=c("Ctrl","AA"))
umap.plot$Tissue <- factor(umap.plot$Tissue,levels=c("PB","BM"))
colnames(umap.plot)[11] <- "Disease"
umap.plot$Disease <- factor(umap.plot$Disease,levels=c("Ctrl","non-SAA","SAA"))

ggplot(umap.plot,aes(x=UMAP1,y=UMAP2))+theme_classic()+
     geom_point(aes(x=UMAP1, y=UMAP2,colour=cellGroup,shape=Tissue),alpha=1,size=0.5)+
     theme(panel.background = element_blank(),panel.border=element_rect(fill='transparent', color='black'))
ggsave(paste(pwd,"/output/UMAP.",cell.label,".cellGroup.uwot.pdf", sep=""), height=2, width=3, device="pdf")
			
ggplot(umap.plot,aes(x=UMAP1,y=UMAP2))+theme_classic()+
     geom_point(aes(x=UMAP1, y=UMAP2,colour=Disease,shape=Tissue),alpha=1,size=0.5)+
     theme(panel.background = element_blank(),panel.border=element_rect(fill='transparent', color='black'))	 
ggsave(paste(pwd,"/output/UMAP.",cell.label,".Disease_types.uwot.pdf", sep=""), height=2, width=2.8, device="pdf")

## Fraction of clusters by patients
suppressPackageStartupMessages(library(dplyr))
T.counts <- sapply(unique(aligned.obj@meta.data$SampleID), function(x){
	 table(as.factor(aligned.obj@meta.data$cellGroup)[ aligned.obj@meta.data$SampleID == x])
     })
T.counts <- as.data.frame(t(prop.table(T.counts, 2)), stringsAsFactors=F)
T.counts.melt <- reshape2::melt(T.counts %>% tibble::rownames_to_column( var = "group" ))
names(T.counts.melt) <- c("Sample", "Cluster", "Percentage")
T.counts.melt$Cluster <- factor(T.counts.melt$Cluster, levels = new.ident.order, ordered=T)
T.counts.melt <- T.counts.melt[order(T.counts.melt$Cluster),]
oneCluster <- "Naive"
T.oneCluster <- T.counts.melt[T.counts.melt$Cluster == oneCluster & !T.counts.melt$Sample %in% c("Ctrl2","Ctrl4"), ]
T.sample.order <- c("Ctrl2","Ctrl4", as.character(T.oneCluster$Sample)[ order(T.oneCluster$Percentage,decreasing=T) ])

ggplot(T.counts.melt, aes(x=factor(Sample, levels=T.sample.order), y=Percentage)) +
	geom_bar(aes(fill=factor(Cluster, levels=new.ident.order)), stat="identity") +
	geom_hline(yintercept=cumsum(rev(T.counts.melt$Percentage[ T.counts.melt$Sample == "Ctrl2"])), linetype="dashed") +
	guides(fill=guide_legend("")) +
	xlab("") + theme(axis.text.x=element_text(angle=90))
ggsave(paste0(pwd, "/output/",cell.label,".clusterFrequency_byPatient.aligned.pdf"))
dev.off()

##  Find markers
##  Conserved cell type markers
# - Find markers in cluster in Ctrl cells, denoted as A
# - Find markers in cluster in AA cells, denoted as B
# - Conserved markers: intersect(A, B)
table(names(aligned.obj@ident)==rownames(aligned.obj@meta.data))
aligned.obj@ident <- factor(aligned.obj@meta.data$cellGroup,levels=new.ident.order)
names(aligned.obj@ident) <- rownames(aligned.obj@meta.data)

conserved.markers <- c()
for(i in levels(aligned.obj@ident)){
	 print(paste0("Processing cluster c", i, " ..."))
	 this.marker <- FindConservedMarkers(
		aligned.obj,
		ident.1 = i,
		grouping.var = "source",
		print.bar = FALSE, only.pos = TRUE
	 )
	 this.marker <- this.marker[order(this.marker[,2] + this.marker[,7], decreasing =T), ]
	 this.marker$gene <- row.names(this.marker)
	 this.marker$cluster <- i
	 conserved.markers <- rbind(conserved.markers, this.marker)
     }

## find all markers
allMarkers <- FindAllMarkers(aligned.obj, only.pos = TRUE)

save(aligned.obj,allMarkers,pc.genes,conserved.markers, new.ident.order, file = paste0(pwd,"/int/04.Seurat.AACtrl.aligned.",cell.label,".Rdata"))

## output marker genes
write.table(allMarkers, file = paste0(pwd,"/output/allMarkers.",cell.label,".xls"),row.names=F,sep="\t",quote=F)

## heatmap of top marker genes
## function to select genes
markerSelect <- function(
     DEGs,
     obj,
     adjP = 0.05,
     fc = log(2),
     only.pos = T,
     gene.num = NULL,
     conserved = F,
     order = F
     ){
     if(isTRUE(conserved)){
         p_val_adj <- 5
         avg_logFC <- 2
         } else {
             p_val_adj <- "p_val_adj"
             avg_logFC <- "avg_logFC"
             }
     topMarkers <- unique(unlist(lapply(unique(DEGs$cluster), function(x){
     this.data <- DEGs[DEGs$cluster == x & DEGs[,p_val_adj] <= adjP & abs(DEGs[,avg_logFC]) >= fc,]
     if (isTRUE(only.pos)){
         this.data <- DEGs[DEGs$cluster == x & DEGs[,p_val_adj] <= adjP & DEGs[,avg_logFC] >= fc,]
         }
     this.data <- this.data[ order(this.data[,avg_logFC], decreasing=T),]
     if(is.null(gene.num)){
         return(this.data$gene)
         }
     return(head(this.data$gene, n=gene.num))
     })))
     if(isTRUE(order)){
         topMarkers <- topMarkers[hclust(as.dist(1-cor(t(obj@scale.data[topMarkers,]))), method="ward.D2")$order]
         }
     return(topMarkers)
     }

## heatmap function
library(pheatmap)
geneHeatmap <- function(
     obj,
     genes,
     label,
     annotate,
	 order,
     minScale = -2.5,
     maxScale = 2.5,
     cluster_rows = T,
     cluster_cols = F,
     cutree_rows = NA,
     cutree_cols = NA,
     clustering_distance_rows = "euclidean",
     clustering_method = "ward.D2",
     show_rownames = T,
     show_colnames= F,
     fontsize_row = 5,
	 width=NA,
	 height=NA
     ){
     sub.data <- as.matrix(obj@scale.data[genes, order(obj@ident)])
     if(length(genes) == 1){
         sub.data <- t(sub.data)
         row.names(sub.data) <- genes
         cluster_rows <- F
         cutree_rows <- NA
         }
     sub.obj <- ScaleData(obj)
     sub.data[ sub.data < minScale] <- minScale
     sub.data[ sub.data > maxScale] <- maxScale
     annotate.df <- sapply(annotate, function(x){
      if( x == "ident"){
         as.character(slot(sub.obj, x))
         } else {
         as.character(sub.obj@meta.data[[x]])
         }
     })
     annotate.df <- as.data.frame(annotate.df, row.names=row.names(sub.obj@meta.data),stringsAsFactors=F)[order(sub.obj@ident),]
     annotate.df[[1]] <- factor(annotate.df[[1]],levels=levels(sub.obj@ident))
	 annotate.df[[2]] <- factor(annotate.df[[2]],levels=order)
	 annotate.df <- annotate.df[order(annotate.df[[1]],annotate.df[[2]]),]
	 sub.data <- sub.data[,rownames(annotate.df)]
     gaps_col <- match(unique(annotate.df[[1]]), annotate.df[[1]]) - 1
     gaps_col <- gaps_col[-1]
     pdf(paste0(pwd,"/output/heatmap.", label, ".pdf"),width=width,height=height)
     pheatmap(
         sub.data,
         color=gplots::colorpanel(n=100, low="#FF00FF",mid="#000000", high="#FFFF00"),
         cluster_rows=cluster_rows, cluster_cols=cluster_cols,
         cutree_rows = cutree_rows, cutree_cols = cutree_cols,
         clustering_distance_rows = clustering_distance_rows,
         clustering_method = clustering_method,
         show_rownames = show_rownames, show_colnames = show_colnames,
         fontsize_row = fontsize_row,
         annotation_col = annotate.df,
         gaps_col = gaps_col
         )
     dev.off()
     }
	
T.topMarkers <- markerSelect(allMarkers,fc = log(1.5), aligned.obj, gene.num=20)
geneHeatmap(obj = aligned.obj,order=c("BM","PB"), genes = T.topMarkers,label=paste0(cell.label,".allMarkers.Tissue"), 
     annotate = list("Cluster"="ident","Tissue"="Tissue"), cutree_rows = 3, cluster_rows=F, fontsize_row = 4,
	 width=6,height=4)

## identified differentially expressed genes (DEGs)
## DEGs between non-SAA and Ctrl
## We exclude Ctrl4 cells from this analysis 
## and cell-cell interaction analysis since
## we detected abnormal activation signaling
## in T cells of this sample.
## reset ident
aligned.obj@meta.data$cellGroup.Disease <- paste0(aligned.obj@meta.data$cellGroup,"_",aligned.obj@meta.data$Disease_types)
aligned.obj <- SetAllIdent(aligned.obj,id = "cellGroup.Disease")
cell.use <- row.names(aligned.obj@meta.data[aligned.obj@meta.data$SampleID != "Ctrl4",])
sub.obj <- SubsetData(aligned.obj, cells.use = cell.use, subset.raw = T)

nonSAA_Ctrl.DEGs <- c()
for(i in unique(sub.obj@meta.data$cellGroup)){
	 print(paste0(i, "_non-SAA vs ", i, "_Ctrl"))
	 this.DEGs <- FindMarkers(sub.obj,ident.1 = paste0(i, "_non-SAA"), ident.2 = paste0(i, "_Ctrl"))
	 this.DEGs$gene <- row.names(this.DEGs)
	 this.DEGs$cluster <- i
	 nonSAA_Ctrl.DEGs <- rbind(nonSAA_Ctrl.DEGs, this.DEGs)
     }

## set thresholds
fc <- 0.25
adjP <- 0.05
nonSAA_Ctrl_DEGsig <- nonSAA_Ctrl.DEGs[ abs(nonSAA_Ctrl.DEGs$avg_logFC) >= fc & nonSAA_Ctrl.DEGs$p_val_adj <= adjP,]
write.table(nonSAA_Ctrl_DEGsig, file = paste0(pwd,"/output/nonSAA_Ctrl.",cell.label,".DEGs.aligned.xls"),row.names=F,sep="\t",quote=F)

## DEGs 
SAA_nonSAA.DEGs <- c()
for(i in unique(aligned.obj@meta.data$cellGroup)){
	 print(paste0(i, "_SAA vs ", i, "_non-SAA"))
	 this.DEGs <- FindMarkers(aligned.obj,ident.1 = paste0(i, "_SAA"), ident.2 = paste0(i, "_non-SAA"))
	 this.DEGs$gene <- row.names(this.DEGs)
	 this.DEGs$cluster <- i
	 SAA_nonSAA.DEGs <- rbind(SAA_nonSAA.DEGs, this.DEGs)
     }

SAA_nonSAA_DEGsig <- SAA_nonSAA.DEGs[ abs(SAA_nonSAA.DEGs$avg_logFC) >= fc & SAA_nonSAA.DEGs$p_val_adj <= adjP,]
write.table(SAA_nonSAA_DEGsig, file = paste0(pwd,"/output/SAA_nonSAA.",cell.label,".DEGs.aligned.xls"),row.names=F,sep="\t",quote=F)
