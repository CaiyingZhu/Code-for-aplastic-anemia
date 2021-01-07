### subroutines
## select high variable genes
pcaGenes <- function(objs){
	 unique(unlist(lapply(objs, function(obj){
		this.pca <- prcomp(t(as.matrix(obj@data)))
		pc.dim <- 1:10
		pc.gene.num <- 200
		this.pca.genes <- unique(unlist(lapply(pc.dim, function(x){
			head(row.names(this.pca$rotation[order(abs(this.pca$rotation[,x]), decreasing=T),]), n=pc.gene.num)
		})))
		return(this.pca.genes)
	     })))
     }

## select marker genes	 
library(pheatmap)
markerSelect <- function(DEGs, obj, pval=0.05, fc=log(2), only.pos=F, gene.num = NULL, conserved=F, order=F){
	 if(isTRUE(conserved)){
		p_val_adj <- 5
		avg_logFC <- 2
	 } else {
		p_val_adj <- "p_val_adj"
		avg_logFC <- "avg_logFC"
	 }
	 topMarkers <- unique(unlist(lapply(unique(DEGs$cluster), function(x){
		this.data <- DEGs[DEGs$cluster == x & DEGs[,p_val_adj] <= pval & abs(DEGs[,avg_logFC]) >= fc,]
		if (isTRUE(only.pos)){
			this.data <- DEGs[DEGs$cluster == x & DEGs[,p_val_adj] <= pval & DEGs[,avg_logFC] >= fc,]
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

## caculate Roe value
Roe <- function(data, condition, cellType, samples, ctrl){
  cellType.class <- unique(data[[cellType]][ data[[samples]] %in% ctrl])
  samples.class <- unique(data[[samples]])
  condition.class <- unique(data[[condition]])
  this.df <- c()
  # O/E
  for(cdt in condition.class){
     cdt.data <- data[ data[[condition]] == cdt & data[[cellType]] %in% cellType.class,]
    ctrl.total <- nrow(cdt.data[ cdt.data[[samples]] %in% ctrl,])
    for(ctyp in cellType.class){
      ctrl.ctyp <- nrow(cdt.data[ cdt.data[[samples]] %in% ctrl & cdt.data[[cellType]] == ctyp, ])
      ratio.expected <- ctrl.ctyp / ctrl.total
      for(ptt in samples.class){
        if(!ptt %in% ctrl){
          ptt.total <- nrow(cdt.data[ cdt.data[[samples]] == ptt, ])
          ptt.ctyp <- nrow(cdt.data[ cdt.data[[samples]] == ptt & cdt.data[[cellType]] == ctyp, ])
          ptt.O2E <- ptt.ctyp / (ptt.total * ratio.expected)
          ppt.chisq.pvalue <- chisq.test(
            matrix(
              c(ctrl.ctyp,
                ctrl.total - ctrl.ctyp,
                ptt.ctyp,
                ptt.total - ptt.ctyp),
              nrow=2
            ))$p.value
          this.df <- rbind(this.df, c(cdt, ctyp, ptt, ptt.O2E, ppt.chisq.pvalue))
        }
      }
    }
  }
  this.df <- as.data.frame(this.df, stringsAsFactors=F)
  colnames(this.df) <- c("condition", "cellType", "samples", "O2E", "chisq.p")
  this.df$O2E <- as.numeric(this.df$O2E)
  this.df$chisq.p <- as.numeric(this.df$chisq.p)
  this.estimate <- c()
  for(cdt in condition.class){
    for(ctyp in cellType.class){
      this.O2E <- as.numeric(this.df$O2E[this.df$condition==cdt & this.df$cellType==ctyp])
      this.O2E <- this.O2E[!is.na(this.O2E)]
      this.mean <- mean(this.O2E)
      this.sem <- sd(this.O2E, na.rm=T)/sqrt(length(this.O2E))
      this.estimate <- rbind(this.estimate, c(cdt, ctyp, this.mean, this.sem))
    }
  }
  this.estimate <- as.data.frame(this.estimate, stringsAsFactors=F)
  colnames(this.estimate) <- c("condition", "cellType", "meanO2E", "sem")
  this.estimate$meanO2E <- as.numeric(this.estimate$meanO2E)
  this.estimate$sem <- as.numeric(this.estimate$sem)
  return(list("samples"=this.df, "group"=this.estimate))
}

# set working directory
setwd("/public/home/zhucy/SingleCell/AA/SCENIC")
pwd <- getwd()
library(reticulate)
use_python("/public/home/zhucy/anaconda3/bin/python")
py_config()
py_available(umap)

# load libraries
suppressPackageStartupMessages({
     library(SCENIC)
     library(AUCell)
     library(RcisTarget)
     library(SingleCellExperiment)
     library(Seurat)
     })  
scenicOptions <- readRDS("int/scenicOptions.Rds")
# high confident AUC 
regulonAUC <- getAUC(loadInt(scenicOptions, "aucell_regulonAUC"))
names(dimnames(regulonAUC)) <- NULL
dim(regulonAUC) # 664 6259
# cell  info
cellInfo <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "cellInfo"), ifNotExists="null")
dim(cellInfo)   # 6259   19

# binary nonDup
regulonBinary <- loadInt(scenicOptions, "aucell_binary_nonDupl")
dim(regulonBinary) # 444 6259
#### Perform a canonical correlation analysis (CCA)
# create Seurat object
tf.obj <- CreateSeuratObject(
 	 raw.data = regulonBinary,
	 meta.data = cellInfo,
	 min.cells = 0, min.features = 0
     )
	 
tf.obj@data <- tf.obj@raw.data
tf.obj@scale.data <- tf.obj@raw.data

tf.activity.obj <- CreateSeuratObject(
     raw.data = regulonAUC[which(row.names(regulonAUC) %in% row.names(regulonBinary)), ],
	 meta.data = cellInfo,
	 min.cells = 0, min.features = 0
     )
tf.activity.obj@data <- tf.activity.obj@raw.data
tf.activity.obj@scale.data <- tf.activity.obj@raw.data

## PCA genes
pc.genes <- pcaGenes(c(tf.obj))
length(pc.genes) # 305

# combine
tf.obj <- RunPCA(tf.obj, pc.genes = pc.genes, pcs.compute = 30, pcs.print=NULL)
p1 <- DimPlot(object = tf.obj, reduction.use = "pca", group.by = "source", pt.size = 0.5, do.return = TRUE, dim.1=1, dim.2=2)
p2 <- DimPlot(object = tf.obj, reduction.use = "pca", group.by = "source", pt.size = 0.5, do.return = TRUE, dim.1=2, dim.2=3)
p3 <- VlnPlot(object = tf.obj, features.plot = "PC1", group.by = "source", do.return = TRUE)
p4 <- VlnPlot(object = tf.obj, features.plot = "PC2", group.by = "source", do.return = TRUE)
plot_grid(p1, p2, p3, p4)
ggsave(paste(pwd, "/output/hspc.PC1_and_PC2.pdf", sep=""), height=5, width=8)

## check for dimension selection
# PCElbowPlot
PCElbowPlot(object = tf.obj)
ggsave(paste(pwd, "/output/hspc.DR.PCElbowPlot.pdf", sep=""))
# heatmap the 20 CC
pdf(paste(pwd, "/output/hspc.PCA_heatmap.pdf", sep=""), height=12, width=7)
DimHeatmap(object = tf.obj, reduction.type = "pca", cells.use = 500, dim.use = 1:20, do.balanced = TRUE)
dev.off()
save(tf.obj,tf.activity.obj,file=paste0(pwd,"/int/tf.activity.obj.Rdata"))

dimSelected <- 10
# umap and Clustering
tf.obj <- RunUMAP(tf.obj, reduction.use = "pca", dims.use = 1:dimSelected, min_dist=0.3)
tf.obj <- FindClusters(tf.obj, reduction.type = "pca", resolution = c(0.2, 0.3, 0.35, 0.4), dims.use = 1:dimSelected, force.recalc =T)

# Visualization
p1 <- DimPlot(tf.obj, reduction.use="umap", do.return = T, pt.size = 0.5, group.by = "source") + xlab("source")
p2 <- DimPlot(tf.obj, reduction.use="umap", do.return = T, pt.size = 0.5, group.by = "SampleID") + xlab("SampleID")
p3 <- DimPlot(tf.obj, reduction.use="umap", do.label = T, do.return = T, pt.size = 0.5, group.by ="res.0.2") + xlab("res_0.2")
p4 <- DimPlot(tf.obj, reduction.use="umap", do.label = T, do.return = T, pt.size = 0.5, group.by ="res.0.3") + xlab("res_0.3")
p5 <- DimPlot(tf.obj, reduction.use="umap", do.label = T, do.return = T, pt.size = 0.5, group.by ="res.0.35") + xlab("res_0.35")
p6 <- DimPlot(tf.obj, reduction.use="umap", do.label = T, do.return = T, pt.size = 0.5, group.by ="res.0.4") + xlab("res_0.4")
plot_grid(p1, p2, p3, p4, p5, p6, nrow=1)
ggsave(paste(pwd, "/output/hspc.umap_clustering.dimUse", dimSelected, ".pdf", sep=""), height=4, width=21, device="pdf")
p7 <- DimPlot(tf.obj, reduction.use="umap", do.label = T, do.return = T, pt.size = 0.5, group.by = "cellGroup", no.legend=T)
plot_grid(p1, p4, p7, nrow=1)
ggsave(paste(pwd, "/output/hspc.umap_clustering.projectionBySeurat.dimUse", dimSelected, ".pdf", sep=""), height=4, width=10, device="pdf")

tf.obj@meta.data$Disease_types <- factor(tf.obj@meta.data$Disease_types,levels=c("Ctrl","non-SAA","SAA"))
DimPlot(tf.obj, reduction.use="umap", do.return = T, pt.size = 0.1, group.by = "Disease_types") + xlab("Disease_types")
ggsave(paste(pwd, "/output/hspc.umap_clustering.projectionBySeurat.dimUse", dimSelected, ".Disease_types.pdf", sep=""), height=2.5, width=3.2, device="pdf")

# regulon defined clusters
tf.obj <- FindClusters(tf.obj, reduction.type = "pca", resolution = 0.3, dims.use = 1:dimSelected, force.recalc =T)

library(scales)
cell_counts.byCluster <- data.frame()
for( i in levels(tf.obj@ident)){
	 this.data <- table(tf.obj@meta.data$cellGroup[tf.obj@ident %in% i])
	 this.data <- data.frame("cellType"=names(this.data), "counts"=as.numeric(this.data), "freq"=as.numeric(this.data)/sum(this.data), "tf.cellType"=i, stringsAsFactors=F)
	 cell_counts.byCluster <- rbind(cell_counts.byCluster, this.data)
     }

# display in pie chart
ident.order <- c(3, 0, 4, 2, 1, 5, 6)
new.ident.order <- levels(tf.obj@meta.data$cellGroup)

cols <- c("#D39200","#93AA00","#00BA38","#FF61C3","#F8766D","#619CFF","#DB72FB","#00C19F","#00B9E3")
cell_counts.byCluster$tf.cellType <- factor(cell_counts.byCluster$tf.cellType, levels=ident.order)
cell_counts.byCluster$cellType <- factor(cell_counts.byCluster$cellType, levels=new.ident.order)
ggplot(cell_counts.byCluster, aes(x="", y=freq, fill=factor(cellType))) + 
	 geom_bar(width = 1, stat = "identity", alpha=0.5) +
	 coord_polar("y", start=0) +
	 xlab("") + ylab("") + 
	 theme(axis.text.x=element_blank()) +
	 facet_wrap(~tf.cellType) +
	 ggtitle("Distribution of cell types\nin tf-network based clusters")
     ggsave(paste0(pwd, "/output/hspc.cell_types_pieChart.newOrder.pdf"))

#### Identify cell type markers
# find all markers
tf.allMarkers <- FindAllMarkers(tf.obj)
table(tf.allMarkers$cluster)
# 0  1  2  3  4  5  6
# 32 55 44 38 42 30 40

## Feature plot of cluster specific markers genes
gene.use <- c("MEIS1 (75g)","HOXA9 (33g)","IRF1 (332g)","POLE3 (562g)","TAL1 (642g)","GATA1 (1052g)","KLF1_extended (1671g)",
     "EBF1 (17g)","TCF3 (218g)","CEBPB (1418g)","SPI1 (666g)","IRF8 (590g)","SPIB (1419g)","MAFB_extended (54g)","CEBPE (48g)")
for(i in gene.use){
     FeaturePlot(object = tf.obj, features.plot = i, cols.use = c("lightgrey","blue"), pt.size = 0.001, reduction.use="umap") +
     theme(axis.text=element_text(size=5,color="black"),axis.text.x=element_text(size=5,color="black"),axis.text.y = element_text(size = 5,color="black"),
     axis.title.y = element_text(size = 5),axis.title.x = element_text(size = 5,color="black"))
     ggsave(paste(pwd, "/output/Umap_clustering.dim.tf.top.markers.",i,".NoLgend.pdf",sep=""), height=2, width=1.8, device="pdf")
     }

## UMAP of regulon clusters
DimPlot(tf.obj, reduction.use="umap", pt.size = 0.1, 
	cols.use=rev(RColorBrewer::brewer.pal(n=8, "YlOrRd")[-1]), no.legend=T, plot.title="regulon cluster") + 
    theme(axis.text=element_text(size=8,color="black"),axis.text.x=element_text(size=8,color="black"),axis.text.y = element_text(size = 8,color="black"),
	axis.title.y = element_text(size = 8),axis.title.x = element_text(size = 8,color="black"))
ggsave(paste0(pwd, "/output/umap_by_regolon_cluster.new.order.pdf"),height=2.4,width=2.2)
 
## View regulon activity
rm(list=ls())
pwd <- getwd()
library(Seurat)
load(paste0(pwd, "/int/05.tf_clustered.Rdata"))

## R1=cluster3; R2=cluster0; R3=cluster4;
## R4=cluster2; R5=cluster1; R6=cluster5;
## R7=cluster6
ident.order <- c(3, 0, 4, 2, 1, 5, 6)
tf.obj@ident <- factor(tf.obj@ident, levels=ident.order)

## UMAP colored by regulon clusters colored and
## disease classification
DimPlot(tf.obj, reduction.use="umap", pt.size = 0.1, label.size = 8,
	 cols.use=rev(RColorBrewer::brewer.pal(n=8, "YlOrRd")[-1]), no.legend=F, plot.title="regulon cluster") +
	 theme(axis.text=element_text(size=8,color="black"),axis.text.x=element_text(size=8,color="black"),axis.text.y = element_text(size = 8,color="black"),
	 axis.title.y = element_text(size = 8),axis.title.x = element_text(size = 8,color="black"))
ggsave(paste0(pwd, "/output/umap_by_regolon_cluster.new.legend.pdf"),height=2.4, width=3)

DimPlot(tf.obj, reduction.use="umap", pt.size = 0.1, group.by = "Disease_types", no.legend=T, plot.title="Disease_types")
ggsave(paste0(pwd, "/output/umap_by_Disease_types.pdf"))


## Differential regulons
tf.obj@meta.data$celltype.source.tf <- paste0(tf.obj@ident, "_", tf.obj@meta.data$Disease_types)
tf.obj <- StashIdent(tf.obj, save.name = "Celltype.clsuter.tf")
tf.obj <- SetAllIdent(tf.obj, id = "celltype.source.tf")

AA_Ctrl.DEGs <- c()
for(t in c("non-SAA","SAA")){
	 for(i in sort(unique(tf.obj@meta.data$Celltype.clsuter.tf))){
	 	 print(paste0(i,"_",t," vs ", i, "_Ctrl"))
		  this.DEGs <- FindMarkers(tf.obj, ident.1 = paste0(i, "_",t), ident.2 = paste0(i, "_Ctrl"), min.pct=0)
		  this.DEGs$gene <- row.names(this.DEGs)
		  this.DEGs$cluster <- i
		  this.DEGs$disease <- t
		  AA_Ctrl.DEGs <- rbind(AA_Ctrl.DEGs, this.DEGs)
         }
	 }
	 
table(AA_Ctrl.DEGs$cluster)
# 0  1  2  3  4  5  6
# 11 23 19 16 13 25 29

## output DEGs
write.table(AA_Ctrl.DEGs, paste(pwd, "/output/tf.cluster_DEG.xls",sep=""), quote=F, sep="\t", row.names=T, col.names=NA)
## save data
save(AA_Ctrl.DEGs, tf.obj, file = paste0(pwd, "/int/06.tf.DEG.Rdata"))

## filter differential regulons
load(paste0(pwd, "/int/06.tf.DEG.Rdata"))
adjP <- 0.05
fc <- 0
AA_Ctrl.DEGs.filtered <- AA_Ctrl.DEGs[ AA_Ctrl.DEGs$p_val_adj <= adjP & abs(AA_Ctrl.DEGs$avg_logFC) >= fc , ]
table(AA_Ctrl.DEGs.filtered$cluster)
# 0  1  2  3  4  5  6
# 11 19 17 14  8 11 10

##  heatmap and summary of DEG frequency in each cell types
for(t in disease){
	 AA_Ctrl.DEGs.this.filtered <- AA_Ctrl.DEGs.filtered[AA_Ctrl.DEGs.filtered$disease == t ,]
	 AA_Ctrl.DEGs.this <- AA_Ctrl.DEGs[AA_Ctrl.DEGs$disease == t,]    
     tf.DEG.freq <- list()
     for( i in unique(AA_Ctrl.DEGs.this.filtered$gene)){
	      this.value <- list()
	      for(j in sort(unique(tf.obj@meta.data$Celltype.clsuter.tf))){
		     this.data <- AA_Ctrl.DEGs.this[ AA_Ctrl.DEGs.this$cluster %in% j, ] 
		     if(i %in% this.data$gene){
			     this.value$log2FC <- c(this.value$log2FC, this.data$avg_logFC[this.data$gene == i])
			     this.value$adj.pvalue <- c(this.value$adj.pvalue, this.data$p_val_adj[this.data$gene == i])
		         } else {
			         this.value$log2FC <- c(this.value$log2FC, 0)
			         this.value$adj.pvalue <- c(this.value$adj.pvalue, 1)
		             }
	             }
	     tf.DEG.freq$log2FC <- rbind(tf.DEG.freq$log2FC, this.value$log2FC)
	     tf.DEG.freq$adj.pvalue <- rbind(tf.DEG.freq$adj.pvalue, this.value$adj.pvalue)
         }
         colnames(tf.DEG.freq$log2FC) <- sort(unique(tf.obj@meta.data$Celltype.clsuter.tf))
         row.names(tf.DEG.freq$log2FC) <- unique(AA_Ctrl.DEGs.this.filtered$gene)
         colnames(tf.DEG.freq$adj.pvalue) <- sort(unique(tf.obj@meta.data$Celltype.clsuter.tf))
         row.names(tf.DEG.freq$adj.pvalue) <- unique(AA_Ctrl.DEGs.this.filtered$gene)
         # corplot
         library(pheatmap)
         data2plot <- tf.DEG.freq$log2FC
         sig.matrix <- matrix(ifelse(tf.DEG.freq$adj.pvalue <= adjP & abs(tf.DEG.freq$log2FC) >= fc, "*", ""), nrow(tf.DEG.freq$log2FC))
         pdf(paste(pwd, "/output/tf.DEG_byCluster.",t,".pdf", sep=""), width=4, height=5)
         pheatmap(data2plot,
         	 main = paste0("FC of ",t," to Ctrl"), 
	         cluster_rows = F, 
	         cluster_cols = F,
	         display_numbers = sig.matrix,
	         number_color = "red",
             border_color = NA,
	         fontsize_row = 4.5,
	         show_rownames = T,
	         show_colnames = T,
	         rot = 0,
	         labels_col = paste(colnames(data2plot), " (", apply(sig.matrix, 2, function(x){sum(x!="")}), ")", sep="")
              )
         dev.off()
         ## save data
         save(AA_Ctrl.DEGs.this.filtered, tf.DEG.freq, file=paste0(pwd, "/int/07.tf.DEG.filtered.",t,".RData"))        	 
         gene.list <- unique(unlist(sapply(ident.order, function(x){
             AA_Ctrl.DEGs.this.filtered$gene[AA_Ctrl.DEGs.this.filtered$cluster %in% x]
             })))
         library(pheatmap)
         data2plot <- data2plot[gene.list, match(ident.order, colnames(data2plot))]
         sig.matrix <- matrix(ifelse(tf.DEG.freq$adj.pvalue <= adjP & abs(tf.DEG.freq$log2FC) >= fc, "*", ""), nrow(tf.DEG.freq$log2FC))
         data2plot[ data2plot < -0.4] <- -0.4
         data2plot[ data2plot > 0.4] <- 0.4
         row.names(data2plot) <- gsub("_extended", "", row.names(data2plot))
         row.names(data2plot) <- gsub(" \\(.*", "", row.names(data2plot))		 
         pdf(paste(pwd, "/output/tf.DEG_byCluster.new.noStar.noNumber.",t,".pdf", sep=""), width=2, height=6)
         pheatmap(data2plot,
         main = paste0("FC of ",t," to Ctrl"), 
         cluster_rows = F, 
         cluster_cols = F,
         number_color = "black",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         border_color = "grey",
         fontsize_row = 8,
         show_rownames = T,
         show_colnames = T,
         rot = 0,
         labels_col = paste(colnames(data2plot), " (", apply(sig.matrix, 2, function(x){sum(x!="")}), ")", sep="")
         )
         dev.off()		 
	  }

## UMAP of cluster 3 and cluster 4
tf.list <- c("GATA2 (680g)", )
for(c in c(3,4)){
     if(c==3){
	     i <- "RUNX3 (185g)"
		 } else {
		 i <- "GATA2 (680g)"
		 }
     sub.obj <- SubsetData(tf.obj, cells.use = rownames(tf.obj@meta.data[tf.obj@meta.data$res.0.3 == c,]),subset.raw=T)
	 FeaturePlot(object = sub.obj,features.plot = i, cols.use = c("lightgrey","blue"), pt.size = 0.1, reduction.use="umap") +
	 theme(axis.text=element_text(size=8,color="black"),axis.text.x=element_text(size=8,color="black"),axis.text.y = element_text(size = 8,color="black"),axis.title.y = element_text(size = 8),axis.title.x = element_text(size = 8,color="black"))
     ggsave(paste0(pwd, "/output/",c,"_", i, "_activity.pdf"),height=2,width=2)
     DimPlot(sub.obj, reduction.use="umap", pt.size = 0.1, group.by = "Disease_types", cols.use=c("#C0C0C0","#7F00FF", "#F44016"), no.legend=T, plot.title="disease") + 
	 theme(axis.text=element_text(size=8,color="black"),axis.text.x=element_text(size=8,color="black"),axis.text.y = element_text(size = 8,color="black"),axis.title.y = element_text(size = 8),axis.title.x = element_text(size = 8,color="black"))
	 ggsave(paste0(pwd, "/output/umap_by_disease.",c,"_",i,".pdf"),height=2, width=2)
	 }

## caculate Roe value	 
sub.obj <- SubsetData(tf.obj, cells.use = rownames(tf.obj@meta.data[!tf.obj@meta.data$SampleID %in% c("Ctrl5","Ctrl6") & tf.obj@meta.data$Disease_types != "SAA",]),subset.raw=T)
hspc.tf.Roe <- Roe(sub.obj@meta.data, condition="Tissue", cellType="res.0.3", samples="SampleID", ctrl=c("Ctrl1","Ctrl4"))
roe.rescale <- hspc.tf.Roe$samples
roe.rescale$O2E[roe.rescale$O2E > 6] <- 6
ident.order <- c(3, 0, 4, 2, 1, 5, 6)
hspc.tf.Roe$group$cellType <- factor(hspc.tf.Roe$group$cellType, levels=ident.order,order=T)
ggplot(data=hspc.tf.Roe$group, aes(x=cellType, y=meanO2E)) +
	geom_bar(stat="identity", alpha=0.8, color="white", fill="grey") +
	geom_point(data=roe.rescale, aes(x=cellType, y=O2E, color=samples, size=-log10(chisq.p)), alpha=0.5) +
	geom_errorbar(aes(ymin=meanO2E-sem, ymax=meanO2E+sem), width=.2,position=position_dodge(.9)) +
	geom_hline(yintercept = 1, color="grey", linetype="dashed") +
	ggtitle("hspc cell type distribution\nin CAA vs. Ctrl (TF based)") +
	guides(fill=guide_legend(title=""), color=guide_legend(title="")) +
	xlab("") + ylab("Ro/e") + coord_flip()
ggsave(paste0(pwd, "/output/hspc.tf.Roe.pdf"), width=8, height=5)
