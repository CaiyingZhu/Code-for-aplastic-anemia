## create seurat object of cells from two healthy control from our recent 
## published study(Xie XW, et al. National Science Review. 2020)
rm(list=ls())
pwd <- getwd()
## set directory containing Seuratv3 package
.libPaths("/public/home/zhucy/software/R-3.5.0/UpdateLibrary")
library(reticulate)
use_python("/public/home/zhucy/anaconda3/bin/python")
py_config()
py_available(umap)
library(Seurat)
library(ggplot2)
library(cowplot)

#### load two normal control expression data
load("/public/home/zhucy/SingleCell/NHSPC/Res/XieXW/Ctrl_XieXW.Rdata")

#### create Seurat object
XieXW.obj <- CreateSeuratObject(counts, project = "XieXW", min.cells = 3, min.features = 200,
	  names.field = 1, names.delim = "_", meta.data = meta)
	
dim(XieXW.obj@assays$RNA@data)   # 16494  1491
dim(XieXW.obj@assays$RNA@counts) # 16494  1491

#### Qaulity control
# mitochondrial abundance
mito.genes <- grep(pattern = "^MT-", x = rownames(x = XieXW.obj@assays$RNA@data), value = TRUE)
percent.mito <- Matrix::colSums(XieXW.obj@assays$RNA@counts[mito.genes, ])/Matrix::colSums(XieXW.obj@assays$RNA@counts)
XieXW.obj <- AddMetaData(XieXW.obj, metadata = percent.mito, col.name = "percent.mito")

# filter cells
fold.sd <- 2
# nFeature_RNA: mean +/- fold.sd * sd, 
# nCount_RNA: mean +/- fold.sd * sd, 
# percent.mito: mean + fold.sd * sd
cells.rm <- unique(unlist(lapply(unique(XieXW.obj@meta.data$Celltype), function(cellIdent){
	this.iterm.data <- XieXW.obj@meta.data[XieXW.obj@meta.data$Celltype == cellIdent,]
	this.iterm.genes <- unlist(lapply(c("nFeature_RNA", "nCount_RNA", "percent.mito"), function(x){
		# for each iterm
		this.mean <- mean(this.iterm.data[[x]], na.rm=T)
		this.sd <- sd(this.iterm.data[[x]], na.rm=T)
		if(x %in% c("nFeature_RNA", "nCount_RNA")){
			this.rm.cells <- which(this.iterm.data[[x]] <= this.mean - fold.sd * this.sd | this.iterm.data[[x]] >= this.mean + fold.sd * this.sd)
		}	else {
			this.rm.cells <- which(this.iterm.data[[x]] >= this.mean + fold.sd * this.sd)
		}
		return(this.rm.cells)
	}))
	return(this.iterm.genes)
})))
length(cells.rm) # 113
XieXW.clean.obj <- subset(XieXW.obj, cells = rownames(XieXW.obj@meta.data)[-cells.rm])

## Cells per sample
table(XieXW.clean.obj@meta.data$SampleID)
# Ctrl5 Ctrl6
#   648   730

## data dimension
dim(XieXW.clean.obj@assays$RNA@data)   # 16494  1378
dim(XieXW.clean.obj@assays$RNA@counts) # 16494  1378

#### Remove ribosome and mitochondria genes
XieXW.clean.obj@assays$RNA@data <- XieXW.clean.obj@assays$RNA@data[grep(pattern="^MT-|^RPL|^RPS", x=row.names(XieXW.clean.obj@assays$RNA@data), value=T, invert = T),]
XieXW.clean.obj@assays$RNA@counts <- XieXW.clean.obj@assays$RNA@counts[row.names(XieXW.clean.obj@assays$RNA@data),]

dim(XieXW.clean.obj@assays$RNA@data)   # 16385  1378
dim(XieXW.clean.obj@assays$RNA@counts) # 16385  1378

#### Normalization and scale data
XieXW.clean.obj <- NormalizeData(object = XieXW.clean.obj, normalization.method = "LogNormalize", scale.factor = 10000)
XieXW.clean.obj <- ScaleData(XieXW.clean.obj,vars.to.regress = c("percent.mito","nCount_RNA"))

## save data
save(XieXW.obj, XieXW.clean.obj,file = paste0(pwd,"/int/01.seurat.data.v3.Ctrl5-6.Rdata"))
