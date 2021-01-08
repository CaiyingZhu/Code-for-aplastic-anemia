## integrate cells from two healthy control from our recent published 
## study(Xie XW, et al. National Science Review. 2020) and cells from
## this study.
rm(list=ls())
pwd <- getwd()
library(reticulate)
use_python("/public/home/zhucy/anaconda3/bin/python")
py_config()
py_available(umap)
library(Seurat)
library(ggplot2)
library(cowplot)
setwd("/public/home/zhucy/SingleCell/AA/Seurat")

##load data
load(paste0(pwd,"/int/02.Celllabel.transfer.XieXW.Rdata"))
load(paste0(pwd,"/int/02.Seurat.Ctrl.AA.aligned.Rdata"))

## update object
reference.obj <- UpdateSeuratObject(hspc.aligned.obj)

## splite AA and two CD34+ control
## add label and split object by dataset
reference.obj@meta.data$dataset <- ""
reference.obj@meta.data$dataset[reference.obj@meta.data$SampleID %in% c("Ctrl1","Ctrl4")] <- "dataset1"
reference.obj@meta.data$dataset[!reference.obj@meta.data$SampleID %in% c("Ctrl1","Ctrl4")] <- "dataset2"
ABC.clean.obj@meta.data$dataset <- "dataset3"
reference.list <- SplitObject(reference.obj, split.by = "dataset")

## integrate data
merge.obj <- merge(reference.obj,ABC.clean.obj)
AACtrl.list <- list("dataset1" = reference.list[[1]],"dataset2" = reference.list[[2]],"dataset" = ABC.clean.obj)

## normalization and find variance of each object
for(i in 1:length(x = AACtrl.list)) {
     AACtrl.list[[i]] <- NormalizeData(object = AACtrl.list[[i]], verbose = FALSE)
     AACtrl.list[[i]] <- FindVariableFeatures(object = AACtrl.list[[i]],selection.method = "vst", nfeatures = 2000, verbose = FALSE)
     }
AACtrl.interAnchor <- FindIntegrationAnchors(object.list = AACtrl.list,scale = TRUE, l2.norm = TRUE,dims = 1:30, k.anchor = 5, k.filter = 200, k.score = 30)
AACtrl.interObj <- IntegrateData(AACtrl.interAnchor, new.assay.name = "integrated", features.to.integrate = rownames(merge.obj@assays$RNA@data),dims = 1:30, k.weight = 100, sd.weight = 1)
DefaultAssay(object = AACtrl.interObj) <- "integrated"
AACtrl.interObj <- FindVariableFeatures(object = AACtrl.interObj, selection.method = "vst", nfeatures = 2000)
AACtrl.interObj <- ScaleData(object = AACtrl.interObj,verbose = FALSE)
AACtrl.interObj <- RunPCA(object = AACtrl.interObj,npcs = 30,features=AACtrl.interObj@assays$integrated@var.features[1:1000], verbose = FALSE)

## estimate and select dimensions for UMAP
plot <- ElbowPlot(AACtrl.interObj,ndims=30)
plot_grid(plot)
ggsave(paste(pwd, "/output/ElbowPlot.integrated.pdf", sep=""))

new.ident.order <- c("HSC/MPP","LMPP","MEP","MLP","EBM","Neu1","Neu2","MD1","MD2")
AACtrl.interObj@active.ident <- factor(AACtrl.interObj@meta.data$cellGroup,levels=new.ident.order)
names(AACtrl.interObj@active.ident) <- rownames(AACtrl.interObj@meta.data)
AACtrl.interObj@meta.data$cellGroup <- factor(AACtrl.interObj@meta.data$cellGroup,levels = new.ident.order)

## select top 10 dims
dim.Use <- 10
cols <- c("#D39200","#93AA00","#00BA38","#FF61C3","#F8766D","#619CFF","#DB72FB","#00C19F","#00B9E3")
AACtrl.interObj <- RunUMAP(object = AACtrl.interObj, reduction = "pca", dims = 1:dim.Use, min.dist = 0.3,spread = 1)
p1 <- DimPlot(object = AACtrl.interObj, cols = c("grey","red","blue"), reduction = "umap", group.by = "dataset")
p2 <- DimPlot(object = AACtrl.interObj, cols = cols,reduction = "umap", group.by = "cellGroup", label = TRUE, repel = TRUE) + NoLegend()
plot_grid(p1, p2)
ggsave(paste(pwd, "/output/Umap_clustering.dim.Use",dim.Use,".integrated.pdf",sep=""), height=4, width=9, device="pdf")

AACtrl.interObj@meta.data$Disease_types <- factor(AACtrl.interObj@meta.data$Disease_types,levels=c("Ctrl","non-SAA","SAA"),order=T)
DimPlot(object = AACtrl.interObj, cols = c("#7F7472","#FFCC99","#D94520"), pt.size = 0.1,reduction = "umap", group.by = "Disease_types")
ggsave(paste(pwd, "/output/Umap_clustering.dim.Use",dim.Use,".integrated.Disease_types.pdf",sep=""), height=2.2, width=3.2, device="pdf")

DimPlot(object = AACtrl.interObj, cols = cols,reduction = "umap", pt.size = 0.1,group.by = "cellGroup", label = TRUE, repel = TRUE)
ggsave(paste(pwd, "/output/Umap_clustering.dim.Use",dim.Use,".integrated.cellGroup.pdf",sep=""), height=2.2, width=3.6, device="pdf")

## only visualized 4 normal control
sub.obj <- subset(AACtrl.interObj,cells=rownames(AACtrl.interObj@meta.data[AACtrl.interObj@meta.data$source=="Ctrl",]))
p1 <- DimPlot(object = sub.obj, cols = cols,reduction = "umap", group.by = "cellGroup", label = TRUE, repel = TRUE) + NoLegend()
p2 <- DimPlot(object = sub.obj, cols = c("#BBBBD2","#D26D15","#7575EB","#912ED2"), reduction = "umap", group.by = "SampleID")
plot_grid(p1, p2)
ggsave(paste(pwd, "/output/Umap_clustering.dim.Use",dim.Use,".integrated.Ctrl.by.SampleID.pdf",sep=""), height=4, width=9, device="pdf")

DefaultAssay(object = AACtrl.interObj) <- "RNA"
Idents(AACtrl.interObj) <- factor(AACtrl.interObj@meta.data$cellGroup)
## save data
save(AACtrl.interObj, file = paste0(pwd,"/int/03.Seurat.integrated.AACtrl.Rdata"))

## identified differentially expressed genes (DEGs)
## Conserved cell type markers
# - Find markers in cluster in Ctrl cells, denoted as A
# - Find markers in cluster in AA cells, denoted as B
# - Conserved markers: intersect(A, B)

conserved.markers <- c()
for(i in levels(AACtrl.interObj@active.ident)){
	 print(paste0("Processing cluster ", i, " ..."))
	 this.marker <- FindConservedMarkers(
		AACtrl.interObj,
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
conserved.markers$cluster <- factor(conserved.markers$cluster,levels=levels(AACtrl.interObj@active.ident),order=T)
table(conserved.markers$cluster)
# HSC/MPP    LMPP     MEP     MLP     EBM    Neu1    Neu2     MD1     MD2
#    261     156     632     534     556      40     697     503     426

library(dplyr)
top12 <- conserved.markers %>% group_by(cluster) %>% top_n(12, Ctrl_avg_logFC)
unigene=unique(top12$gene)
unigene=as.character(unigene)
length(unigene) # 102

num.topMarkers <- c(0,4,8)
AACtrl.interObj.plot <- AACtrl.interObj
data <- as.matrix(AACtrl.interObj.plot@assays$RNA@data)
tab <- data[unigene,]
melt <- reshape2::melt(tab)

## visualization distribution of expression value
ggplot(melt,aes(x=value)) + geom_line(stat="density",colour="black")
ggsave(filename =paste0(pwd,"/output/expr_distribution.expr.pdf"),height=3,width=6,onefile=F)
     
data[data >= 3] <- 3
dgdata <- as(as.matrix(data), "dgCMatrix")
AACtrl.interObj.plot@assays$RNA@data <- dgdata

for(i in levels(AACtrl.interObj@active.ident)){
    for(j in num.topMarkers){
	     FeaturePlot(object = AACtrl.interObj.plot, features = head(conserved.markers$gene[conserved.markers$cluster ==i],12)[(j+1):(j+4)], min.cutoff = "q8", cols = c("lightgrey","blue"), pt.size = 0.25)
	     ggsave(paste(pwd, "/output/hspc.topFeaturePlot.",sub("/",".",i),".subset.", j, ".pdf", sep=""), height=6, width=7, device="pdf")
         }
     }

## surface markers used for sorting
surface.markers <- c("CD34","CD38","FLT3","PTPRC","MME","ITGA6","THY1")
FeaturePlot(object = AACtrl.interObj.plot, features = surface.markers, ncol = 4,label.size=2, min.cutoff = "q8", cols = c("lightgrey","blue"), pt.size = 0.05)
ggsave(paste(pwd, "/output/hspc.topFeaturePlot.surface.markers.pdf", sep=""), height=5, width=11, device="pdf")

## top cluster markers genes    
select.markers <- c("CD52","SPINK2","GATA1","VPREB1","CLC","MPO","S100A9","IRF8","TGFBI")
FeaturePlot(object = AACtrl.interObj.plot, features = select.markers, ncol = 3 , min.cutoff = "q8", cols = c("lightgrey","blue"), pt.size = 0.1) + NoLegend()
ggsave(paste(pwd, "/output/hspc.topFeaturePlot.selected.markers.top10.pdf", sep=""), height=8, width=9, device="pdf")

## identified DEGs in each cluster
hspc.allMarkers <- FindAllMarkers(AACtrl.interObj, only.pos = TRUE)

## subroutines
## select top markers
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
## heatmap 
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
	 obj <- ScaleData(obj)
     sub.data <- as.matrix(obj@assays$RNA@scale.data[genes, order(obj@active.ident)])
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
      if( x == "active.ident"){
         as.character(slot(sub.obj, x))
         } else {
         as.character(sub.obj@meta.data[[x]])
         }
     })
     annotate.df <- as.data.frame(annotate.df, row.names=row.names(sub.obj@meta.data),stringsAsFactors=F)[order(sub.obj@active.ident),]
     annotate.df[[1]] <- factor(annotate.df[[1]],levels=levels(sub.obj@active.ident))
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

library(pheatmap)	 
topMarkers <- markerSelect(hspc.allMarkers, AACtrl.interObj, gene.num=10)
geneHeatmap(obj = AACtrl.interObj, order = c("Ctrl","AA"), genes = topMarkers,label = "hspc.topMarkers", 
      annotate = list("Cluster" = "active.ident","Source" = "source"), cutree_rows = 3, cluster_rows = F,
	  width = 4,height = 2.4,fontsize_row = 3)

## filter out not significant genes
hspc.allMarkers_filter <- hspc.allMarkers[hspc.allMarkers$p_val_adj <= 0.05,]
table(hspc.allMarkers_filter$cluster)
# HSC/MPP    LMPP     MEP     MLP     EBM    Neu1    Neu2     MD1     MD2
#    493     243     819     369     588     108    1018     721     745

## write output
write.table(hspc.allMarkers_filter, file = paste(pwd, "/output/hspc.allMarkers.xls", sep=""), quote=F, sep="\t", row.names=F)

## identify DEGs within each cluster
## reset idents
AACtrl.interObj@meta.data$cellGroup.Disease <- paste0(AACtrl.interObj@meta.data$cellGroup, "_", AACtrl.interObj@meta.data$Disease_types)
Idents(AACtrl.interObj) <- factor(AACtrl.interObj@meta.data$cellGroup.Disease)

## set thresholds
fc <- log(2)
adjP <- 0.05

## DEGs between SAA/non-SAA and Ctrl
for(j in c("non-SAA","SAA")){
     AA_Ctrl.DEGs <- c()
     for(i in new.ident.order){
		 print(paste0(i, "_",j," vs ", i, "_Ctrl"))
		 this.DEGs <- FindMarkers(AACtrl.interObj,assay="RNA",ident.1 = paste0(i, "_",j), ident.2 = paste0(i, "_Ctrl"))
		 this.DEGs$Gene <- row.names(this.DEGs)
		 this.DEGs$cluster <- i
		 AA_Ctrl.DEGs <- rbind(AA_Ctrl.DEGs, this.DEGs)
  	 	 }
	 AA_Ctrl.DEGs$Disease <- j
	 AA_Ctrl.DEGs$cluster <- factor(AA_Ctrl.DEGs$cluster,levels = new.ident.order)
	 AA_Ctrl.DEGs <- AA_Ctrl.DEGs[order(AA_Ctrl.DEGs$cluster,AA_Ctrl.DEGs$avg_logFC),]
	 ## save DEGs
	 save(AA_Ctrl.DEGs, file = paste0(pwd,"/output/04.hspc.",j,"_vs_Ctrl_DEG.Rdata"))
	 ## output DEGs
     write.table(AA_Ctrl.DEGs, file = paste(pwd, "/output/hspc.cluster_",j,"_vs_Ctrl_DEG.xls", sep=""), quote=F, sep="\t", row.names=T, col.names=NA)
     AA_Ctrl.DEGs <- AA_Ctrl.DEGs[ abs(AA_Ctrl.DEGs$avg_logFC) >= fc & AA_Ctrl.DEGs$p_val_adj <= adjP,]
     geneList <- unique(AA_Ctrl.DEGs$Gene)
     DEGedges <- data.frame("Gene"=geneList, stringsAsFactors=F)
     for(i in new.ident.order){
	     genes.up <- AA_Ctrl.DEGs$Gene[ AA_Ctrl.DEGs$cluster == i & AA_Ctrl.DEGs$avg_logFC > 0 ]
	     genes.down <- AA_Ctrl.DEGs$Gene[ AA_Ctrl.DEGs$cluster == i & AA_Ctrl.DEGs$avg_logFC < 0 ]
	     DEGedges[, i] <- 0
	     DEGedges[ DEGedges$Gene %in% genes.up, i] <- 1
	     DEGedges[ DEGedges$Gene %in% genes.down, i] <- -1
         }
     DEGedges <- DEGedges[order(DEGedges[,2], DEGedges[,3],DEGedges[,4], DEGedges[,5], DEGedges[,6], DEGedges[,7], DEGedges[,8], DEGedges[,9], DEGedges[,10], decreasing=T),]
     row.names(DEGedges) <- DEGedges$Gene
     DEGedges$Gene <- NULL
	 ## heatmap of DEGs
     library(pheatmap)
     pdf(paste(pwd, "/output/hspc.DEGs_overlap_in_cell_types.",j,".vs.Ctrl.pdf", sep=""), width=3, height=4.5)
     pheatmap(DEGedges, cluster_rows = F, cluster_cols = F,
     	 color=colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(3),
     	 border_color = NA, show_rownames = F,
         main = "DEGs across HSPCs")
         dev.off()
	 DEGs.num <- sapply(new.ident.order, function(x){
	 genes.up <- AA_Ctrl.DEGs$Gene[ AA_Ctrl.DEGs$cluster == x & AA_Ctrl.DEGs$avg_logFC > 0]
	 genes.down <- AA_Ctrl.DEGs$Gene[ AA_Ctrl.DEGs$cluster == x & AA_Ctrl.DEGs$avg_logFC < 0]
	 c(length(genes.up), length(genes.down))
         })
	 DEGs.num.origin <- as.data.frame(t(DEGs.num))
	 row.names(DEGs.num) <- c(j, "Ctrl")
	 DEGs.num <- as.data.frame(t(DEGs.num))
	 ## heatmap of DEG counts
	 pdf(paste(pwd, "/output/hspc.DEGs_num_",j,"_vs_Ctrl.pdf", sep=""), width=3, height=5)
	 pheatmap(DEGs.num, cluster_rows = F, cluster_cols = F,
	 	 display_numbers = DEGs.num.origin, number_format = "%.0f", fontsize_number = 13,
		 legend = F, border_color = "white",
		 fontsize = 8,
		 main = "Number of DEGs")
	     dev.off()
         }
dev.off()

## DEGs frequency across nine cell clusters
pwd <- getdw()
load(paste0(pwd,"/output/04.hspc.non-SAA_vs_Ctrl_DEG.Rdata"))
nonSAA_DEGs <- AA_Ctrl.DEGs

load(paste0(pwd,"/output/04.hspc.SAA_vs_Ctrl_DEG.Rdata"))
SAA_DEGs <- AA_Ctrl.DEGs
DEGs <- rbind(nonSAA_DEGs,SAA_DEGs)
fc <-log(2)
adjP <- 0.05

AA_Ctrl.DEGs <- DEGs[ abs(DEGs$avg_logFC) >= fc & DEGs$p_val_adj <= adjP,]
geneList <- unique(AA_Ctrl.DEGs$Gene)
DEGedges <- data.frame("Gene"=geneList, stringsAsFactors=F)
for(i in new.ident.order){
     genes.up <- AA_Ctrl.DEGs$Gene[ AA_Ctrl.DEGs$cluster == i & AA_Ctrl.DEGs$avg_logFC > 0 ]
     genes.down <- AA_Ctrl.DEGs$Gene[ AA_Ctrl.DEGs$cluster == i & AA_Ctrl.DEGs$avg_logFC < 0 ]
     DEGedges[, i] <- 0
     DEGedges[ DEGedges$Gene %in% genes.up, i] <- 1
     DEGedges[ DEGedges$Gene %in% genes.down, i] <- -1
     }
DEGedges <- DEGedges[order(DEGedges[,2], DEGedges[,3],DEGedges[,4], DEGedges[,5], DEGedges[,6], DEGedges[,7], DEGedges[,8], DEGedges[,9], DEGedges[,10], decreasing=T),]
row.names(DEGedges) <- DEGedges$Gene
DEGedges$Gene <- NULL
library(pheatmap)
pdf(paste(pwd, "/output/hspc.DEGs_overlap_in_cell_types.AA.vs.Ctrl.pdf", sep=""), width=3, height=4.5)
pheatmap(DEGedges, cluster_rows = F, cluster_cols = F,
 	 color=colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(3),
 	 border_color = NA, show_rownames = F,
     main = "DEGs across HSPCs")
     dev.off()
DEGedges$Freq <- apply(DEGedges, 1, function(x){sum(x!=0)})
pdf(paste(pwd, "/output/hspc.DEGs_overlap_frequency_distribution.AA_vs_Ctrl.pdf", sep=""),width=4, height=5)
freq.hist <- hist(DEGedges$Freq, col="blue", xlab="Number of cell types",ylab="Number of DEGs", las=1, main="", border="white")
lines(x=1:(max(DEGedges$Freq)-1) + 0.5, y=cumsum(freq.hist$counts) / sum(freq.hist$counts) * max(freq.hist$counts), lwd=2, col="black", type="s")
text(x=1:(max(DEGedges$Freq)-1) + 0.5, y=freq.hist$counts + 3, adj=c(0.5,0), xpd=T,paste0(format(freq.hist$counts/sum(freq.hist$counts) * 100, digits=2), "%"),cex=0.7)
dev.off()
