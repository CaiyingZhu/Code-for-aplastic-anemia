## R version 3.5.0; Seurat version 2.3.4
rm(list=ls())
pwd <- getwd()
library(Seurat)
library(dplyr)
library(Hmisc)
#### import UMI table, 'genes' x 'cells'
## read in ctrl counts data and meta table
ctrl.file <- "../data/Ctrl_umi_Codinggene_Tcell.txt"
ctrl.table <- read.table(ctrl.file, head=T, stringsAsFactors=F)

# remove 1st column, GeneName
ctrl.table <- ctrl.table[,-1]
# rename row names 
row.names(ctrl.table) <- ctrl.table[["GeneSymbol"]]

# remove GeneSymbol
ctrl.table$GeneSymbol <- NULL

# assign row names
row.names(meta.table) <- meta.table[["CellID"]]
ctrl.meta <- meta.table[colnames(ctrl.table),]

## AA cells 
AA.file <- "../data/AA_umi_Codinggene_Tcell.txt"
AA.table <- read.table(AA.file, head=T, stringsAsFactors=F)
# remove 1st column, GeneName
AA.table <- AA.table[,-1]

row.names(AA.table) <- AA.table[["GeneSymbol"]]
# remove GeneSymbol
AA.table$GeneSymbol <- NULL

#### import meta data for cells
meta.file <- "../data/meta_Table_Tcell.txt"
meta.table <- read.table(meta.file, head=T, stringsAsFactors=F)

## remove Respones
treat.cellID <- meta.table$CellID[ meta.table$Treatment != "untreated"]
AA.table <- AA.table[, ! names(AA.table) %in% treat.cellID]
meta.table <- meta.table[ meta.table$Treatment == "untreated",]

#### create Seurat object
ctrl.obj <- CreateSeuratObject(
	 raw.data = as.matrix(ctrl.table), min.cells = 3, min.genes = 200, 
	 meta.data = ctrl.meta, project = "Ctrl"
     )

AA.obj <- CreateSeuratObject(
	 raw.data = AA.table, min.cells = 3, min.genes = 200, 
	 meta.data = meta.table, project = "AA"
     )

#### Qaulity control
## Ctrl cells
# mitochondrial abundance
mito.genes <- grep(pattern = "^MT-", x = rownames(x = ctrl.obj@data), value = TRUE)
percent.mito <- Matrix::colSums(ctrl.obj@raw.data[mito.genes, ])/Matrix::colSums(ctrl.obj@raw.data)
ctrl.obj <- AddMetaData(object = ctrl.obj, metadata = percent.mito, col.name = "percent.mito")

## filter cells by cell type
## separately filter out Ctrl2
## and Ctrl4 T cells.
fold.sd <- 2
# nGene: mean +/- fold.sd * sd, 
# nUMI: mean +/- fold.sd * sd, 
# percent.mito: mean + fold.sd * sd
cells.rm <- unique(unlist(lapply(levels(ctrl.obj@ident), function(cellIdent){
	 # for each cell type
	 this.iterm.data <- ctrl.obj@meta.data[ctrl.obj@ident == cellIdent,]
	 this.iterm.genes <- unlist(lapply(c("nGene", "nUMI", "percent.mito"), function(x){
		# for each iterm
		this.mean <- mean(this.iterm.data[[x]], na.rm=T)
		this.sd <- sd(this.iterm.data[[x]], na.rm=T)
		if(x %in% c("nGene", "nUMI")){
			this.rm.cells <- which(this.iterm.data[[x]] <= this.mean - fold.sd * this.sd | this.iterm.data[[x]] >= this.mean + fold.sd * this.sd)
		}	else {
			this.rm.cells <- which(this.iterm.data[[x]] >= this.mean + fold.sd * this.sd)
		}
		return(this.rm.cells)
	 }))
	 return(this.iterm.genes)
     })))
ctrl.clean.obj <- SubsetData(object = ctrl.obj, cells.use = ctrl.obj@cell.names[-cells.rm],subset.raw = TRUE)

## AA cells
# mitochondrial abundance
mito.genes <- grep(pattern = "^MT-", x = rownames(x = AA.obj@data), value = TRUE)
percent.mito <- Matrix::colSums(AA.obj@raw.data[mito.genes, ])/Matrix::colSums(AA.obj@raw.data)
AA.obj <- AddMetaData(object = AA.obj, metadata = percent.mito, col.name = "percent.mito")
# filter cells by cell type
fold.sd <- 2
# nGene: mean +/- fold.sd * sd, 
# nUMI: mean +/- fold.sd * sd, 
# percent.mito: mean + fold.sd * sd
cells.rm <- unique(unlist(lapply(levels(AA.obj@ident), function(cellIdent){
	# for each cell type
	this.iterm.data <- AA.obj@meta.data[AA.obj@ident == cellIdent,]
	this.iterm.genes <- unlist(lapply(c("nGene", "nUMI", "percent.mito"), function(x){
		# for each iterm
		this.mean <- mean(this.iterm.data[[x]], na.rm=T)
		this.sd <- sd(this.iterm.data[[x]], na.rm=T)
		if(x %in% c("nGene", "nUMI")){
			this.rm.cells <- which(this.iterm.data[[x]] <= this.mean - fold.sd * this.sd | this.iterm.data[[x]] >= this.mean + fold.sd * this.sd)
		}	else {
			this.rm.cells <- which(this.iterm.data[[x]] >= this.mean + fold.sd * this.sd)
		}
		return(this.rm.cells)
	}))
	return(this.iterm.genes)
})))
AA.clean.obj <- SubsetData(object = AA.obj, cells.use = AA.obj@cell.names[-cells.rm],subset.raw = TRUE)

# plot "nGene", "nUMI", "percent.mito"
ctrl.nGene <- VlnPlot(object = ctrl.clean.obj, features.plot = c("nGene"), do.return = T) + ylim(0, 3000) + xlab("Ctrl")
ctrl.nUMI <- VlnPlot(object = ctrl.clean.obj, features.plot = c("nUMI"), do.return = T) + ylim(0, 3*10^4) + xlab("Ctrl")
ctrl.mito <- VlnPlot(object = ctrl.clean.obj, features.plot = c("percent.mito"), do.return=T) + ylim(0,0.2) + xlab("Ctrl")
AA.nGene <- VlnPlot(object = AA.clean.obj, features.plot = c("nGene"), do.return = T) + ylim(0, 3000) + xlab("AA")
AA.nUMI <- VlnPlot(object = AA.clean.obj, features.plot = c("nUMI"), do.return = T) + ylim(0, 3*10^4) + xlab("AA")
AA.mito <- VlnPlot(object = AA.clean.obj, features.plot = c("percent.mito"), do.return = T) + ylim(0,0.2) + xlab("AA")

plot_grid(ctrl.nGene, AA.nGene)
ggsave(paste(pwd, "/output/nGene.pdf", sep=""))
plot_grid(ctrl.nUMI, AA.nUMI)
ggsave(paste(pwd, "/output/nUMI.pdf", sep=""))
plot_grid(ctrl.mito, AA.mito)
ggsave(paste(pwd, "/output/rMito.pdf", sep=""))

#### Remove ribosome and mitochondria genes
ctrl.clean.obj@data <- ctrl.clean.obj@data[grep(pattern="^MT-|^RPL|^RPS", x=row.names(ctrl.clean.obj@data), value=T, invert = T),]
ctrl.clean.obj@raw.data <- ctrl.clean.obj@raw.data[row.names(ctrl.clean.obj@data),]
AA.clean.obj@data <- AA.clean.obj@data[grep(pattern="^MT-|^RPL|^RPS", x=row.names(AA.clean.obj@data), value=T, invert = T),]
AA.clean.obj@raw.data <- AA.clean.obj@raw.data[row.names(AA.clean.obj@data),]

#### Normalization and scale data
ctrl.clean.obj <- NormalizeData(object = ctrl.clean.obj, normalization.method = "LogNormalize", scale.factor = 10000)
ctrl.clean.obj <- ScaleData(ctrl.clean.obj, display.progress = F)
#### Normalization and scale data
AA.clean.obj <- NormalizeData(object = AA.clean.obj, normalization.method = "LogNormalize", scale.factor = 10000)
AA.clean.obj <- ScaleData(AA.clean.obj, display.progress = F)

## save data
save(ctrl.obj, ctrl.clean.obj, AA.obj, AA.clean.obj, file=paste0(pwd,"/int/01.seurat.data.Rdata"))

#### Clustering of Ctrl cells
## load data
rm(list=ls())
pwd <- getwd()
library(Seurat)
load(paste0(pwd,"/int/01.seurat.data.Rdata"))

#cell.label <- "CD4T"
cell.label <- "CD8T"
# pca based gene selection
cell.idx <- which(ctrl.clean.obj@meta.data$Celltype == cell.label)
sub.obj <- SubsetData(
	 ctrl.clean.obj, 
	 cells.use = row.names(ctrl.clean.obj@meta.data)[cell.idx]
     )	 

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
pc.genes <- pcaGenes(c(sub.obj))

## Dimension reduction using PCA
sub.obj <- RunPCA(object = sub.obj, pc.genes = pc.genes, do.print = FALSE, pcs.compute = 60)

## check for dimension choose
# PCElbowPlot
PCElbowPlot(object = sub.obj)
ggsave(paste0(pwd, "/output/", cell.label, ".PCElbowPlot.pdf"))

#-- PCA suggest PC 1 to 5
tsne.perplexity <- 30
pca.pcs <- 5
sub.obj <- RunTSNE(sub.obj, reduction.use = "pca", 
	     dims.use = 1:pca.pcs, do.fast = T, dim.embed = 2, 
	     perplexity = tsne.perplexity)
sub.obj <- FindClusters(sub.obj, reduction.type = "pca", 
#	     resolution = 0.3,  ## CD4T
	     resolution = 0.35,  ## CD8T	
	     dims.use = 1:pca.pcs, force.recalc =T)

# Visualization
p0 <- TSNEPlot(sub.obj, do.label = T, do.return = T, pt.size = 0.5, group.by = "Tissue") + xlab("Tissue")
p1 <- TSNEPlot(sub.obj, do.label = T, do.return = T, pt.size = 0.5, group.by = "LibraryID") + xlab("LibraryID")
p2 <- TSNEPlot(sub.obj, do.label = T, do.return = T, pt.size = 0.5, group.by = "res.0.3") + xlab("res.0.3")
plot_grid(p0, p1, p2, nrow=1)
ggsave(paste(pwd, "/output/", cell.label, ".PCAonly.tSNE_clustering.dimUse", pca.pcs, ".prePlxty_",tsne.perplexity, ".pdf", sep=""), height=4, width=10, device="pdf")
save(sub.obj, pc.genes, file=paste0(pwd, "/int/02.", cell.label, ".PCA_clustering.RData"))

#### Cell assignment from AA to Ctrl
rm(list=ls())
pwd <- getwd()
library(Seurat)
load(paste0(pwd,"/int/01.seurat.data.Rdata"))

#cell.label <- "CD4T"
cell.label <- "CD8T"
load(paste0(pwd, "/int/02.", cell.label, ".PCA_clustering.RData"))

# knn
library(class)
Knum <- 10
train.set <- sub.obj@data[pc.genes,]
train.cl <- sub.obj@ident

aa.obj <- SubsetData(AA.clean.obj, cells.use=row.names(AA.clean.obj@meta.data)[AA.clean.obj@meta.data$Celltype %in% cell.label])
aa.set <- aa.obj@data[pc.genes,]
aa.cl <- knn(train=t(train.set), test=t(aa.set), cl=train.cl, prob=T, k=Knum)
aa.obj@meta.data$assigned <- aa.cl
sub.obj@meta.data$source <- "Ctrl"
sub.obj@meta.data$assigned <- train.cl
aa.obj@meta.data$source <- "AA"
mg.obj <- MergeSeurat(object1=sub.obj, object2=aa.obj)
mg.obj <-  ScaleData(mg.obj, display.progress = F)
save(mg.obj, file=paste0(pwd, "/int/03.", cell.label, ".cell_assignment.Rdata"))
