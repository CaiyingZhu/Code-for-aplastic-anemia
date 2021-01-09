## quality control on single cells and 
## create seurat objects of control and AA
rm(list=ls())
pwd <- getwd()
library(reticulate)
library(Seurat) # version 2
library(ggplot2)
library(cowplot)

#### import UMI table, 'genes' x 'cells'
## Control normal cells
ctrl.file <- "../data/Ctrl_codingGene_NHSPC.txt"
ctrl.table <- read.table(ctrl.file, head=T, stringsAsFactors=F)

# remove 1st column, GeneName, rename row names 
ctrl.table <- ctrl.table[,-1]
row.names(ctrl.table) <- ctrl.table[["GeneSymbol"]]

# remove GeneSymbol
ctrl.table$GeneSymbol <- NULL

## AA cells 
AA.file <- "../data/AA_codingGene_HSPC.txt"
AA.table <- read.table(AA.file, head=T, stringsAsFactors=F)
# remove 1st column, GeneName, rename row names 
AA.table <- AA.table[,-1]
row.names(AA.table) <- AA.table[["GeneSymbol"]]

# remove GeneSymbol
AA.table$GeneSymbol <- NULL

#### import meta data for cells
meta.file <- "../data/metaTable_HSPC.txt"
meta.table <- read.table(meta.file, head=T, stringsAsFactors=F)
# assign row names
row.names(meta.table) <- meta.table[["CellID"]]

## remove cells colected from samples response to treatment
treat.cellID <- meta.table$CellID[ meta.table$Treatment != "untreated"]
AA.table <- AA.table[, ! names(AA.table) %in% treat.cellID]
meta.table <- meta.table[ meta.table$Treatment == "untreated",]

## seperate ctrl and AA meta data
ctrl.meta <- meta.table[meta.table$SampleID == c("Ctrl1","Ctrl4"),]
AA.meta <- meta.table[colnames(AA.table),]

#### create Seurat object
ctrl.obj <- CreateSeuratObject(
	raw.data = as.matrix(ctrl.table), min.cells = 3, min.genes = 200, 
	meta.data = ctrl.meta, project = "Ctrl")
	
AA.obj <- CreateSeuratObject(
	raw.data = as.matrix(AA.table), min.cells = 3, min.genes = 200, 
	meta.data = meta.table, project = "AA")
	
#### Qaulity control
## Ctrl cells
# mitochondrial abundance
mito.genes <- grep(pattern = "^MT-", x = rownames(x = ctrl.obj@data), value = TRUE)
percent.mito <- Matrix::colSums(ctrl.obj@raw.data[mito.genes, ])/Matrix::colSums(ctrl.obj@raw.data)
ctrl.obj <- AddMetaData(object = ctrl.obj, metadata = percent.mito, col.name = "percent.mito")

# filter cells
fold.sd <- 2
# nGene: mean +/- fold.sd * sd, 
# nUMI: mean +/- fold.sd * sd, 
# percent.mito: mean + fold.sd * sd
cells.rm <- unique(unlist(lapply(levels(ctrl.obj@ident), function(cellIdent){
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
ctrl.clean.obj <- SubsetData(object = ctrl.obj, cells.use = ctrl.obj@cell.names[-cells.rm],subset.raw =T)

## AA cells
# mitochondrial abundance
mito.genes <- grep(pattern = "^MT-", x = rownames(x = AA.obj@data), value = TRUE)
percent.mito <- Matrix::colSums(AA.obj@raw.data[mito.genes, ])/Matrix::colSums(AA.obj@raw.data)
AA.obj <- AddMetaData(object = AA.obj, metadata = percent.mito, col.name = "percent.mito")
# filter cells
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
AA.clean.obj <- SubsetData(object = AA.obj, cells.use = AA.obj@cell.names[-cells.rm])

# plot "nGene", "nUMI", "percent.mito"
ctrl.nGene <- VlnPlot(object = ctrl.clean.obj, features = c("nGene"), group.by="SampleID") + ylim(0, 3000) + xlab("Ctrl")
ctrl.nUMI <- VlnPlot(object = ctrl.clean.obj, features = c("nUMI"), group.by="SampleID") + ylim(0, 3*10^4) + xlab("Ctrl")
ctrl.mito <- VlnPlot(object = ctrl.clean.obj, features = c("percent.mito"), group.by="SampleID") + ylim(0,0.2) + xlab("Ctrl")
AA.nGene <- VlnPlot(object = AA.clean.obj, features = c("nGene"), group.by="SampleID") + ylim(0, 3000) + xlab("AA")
AA.nUMI <- VlnPlot(object = AA.clean.obj, features = c("nUMI"), group.by="SampleID") + ylim(0, 3*10^4) + xlab("AA")
AA.mito <- VlnPlot(object = AA.clean.obj, features = c("percent.mito"), group.by="SampleID") + ylim(0,0.2) + xlab("AA")

plot_grid(ctrl.nGene, AA.nGene)
ggsave(paste(pwd, "/output/nGene.pdf", sep=""),width=8,height=4)
plot_grid(ctrl.nUMI, AA.nUMI)
ggsave(paste(pwd, "/output/nUMI.pdf", sep=""),width=8,height=4)
plot_grid(ctrl.mito, AA.mito)
ggsave(paste(pwd, "/output/rMito.pdf", sep=""),width=8,height=4)

## remove samples with less than 20 cells
## we also discard the sample with less than 20 cells of HSPCs and without T cells
AA.clean.obj <- SubsetData(object =AA.clean.obj, cells.use = AA.clean.obj@cell.names[! AA.clean.obj@meta.data$SampleID %in% c("P16","P0")], subset.raw =T)

#### Remove ribosome and mitochondria genes
ctrl.clean.obj@data <- ctrl.clean.obj@data[grep(pattern="^MT-|^RPL|^RPS", x=row.names(ctrl.clean.obj@data), value=T, invert = T),]
ctrl.clean.obj@raw.data <- ctrl.clean.obj@raw.data[row.names(ctrl.clean.obj@data),]
AA.clean.obj@data <- AA.clean.obj@data[grep(pattern="^MT-|^RPL|^RPS", x=row.names(AA.clean.obj@data), value=T, invert = T),]
AA.clean.obj@raw.data <- AA.clean.obj@raw.data[row.names(AA.clean.obj@data),]

## save data
save(ctrl.obj, ctrl.clean.obj, file = paste0(pwd,"/int/01.seurat.data.Ctrl.Rdata"))
save(AA.obj, AA.clean.obj, file = paste0(pwd,"/int/01.seurat.data.AA.Rdata"))

