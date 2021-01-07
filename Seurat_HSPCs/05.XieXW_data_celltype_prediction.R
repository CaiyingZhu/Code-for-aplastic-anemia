## cell type prediction
rm(list=ls())
pwd <- getwd()
library(reticulate)
use_python("/public/home/zhucy/anaconda3/bin/python")
py_config()
py_available(umap)
library(Seurat)
library(dplyr)
library(Matrix)
library(uwot)
library(ggplot2)
library(sctransform)

## cell clusters in AA and two normal control were used as reference
## load reference Object
load(paste0(pwd,"/int/02.Seurat.Ctrl.AA.aligned.Rdata"))

## load query Seurat Object
load(paste0(pwd,"/int/01.seurat.data.v3.Ctrl5-6.Rdata"))
## Update object
reference.obj <- UpdateSeuratObject(hspc.aligned.obj)

## add label
XieXW.clean.obj@meta.data$source <- "Ctrl"
XieXW.clean.obj@meta.data$dataset <- "XieXW"
reference.obj@meta.data$dataset <- "reference"

XieXW.clean.obj@meta.data$Sex[XieXW.clean.obj@meta.data$SampleID=="Ctrl5"] <- "Male"
XieXW.clean.obj@meta.data$Sex[XieXW.clean.obj@meta.data$SampleID=="Ctrl6"] <- "Female"
XieXW.clean.obj@meta.data$Age[XieXW.clean.obj@meta.data$SampleID=="Ctrl5"] <- 33
XieXW.clean.obj@meta.data$Age[XieXW.clean.obj@meta.data$SampleID=="Ctrl6"] <- 29

inter.list <- list("reference" = reference.obj, "XieXW" = XieXW.clean.obj)
for(i in 1:length(inter.list)) {
     inter.list[[i]] <- NormalizeData(inter.list[[i]], verbose = FALSE)
     inter.list[[i]] <- ScaleData(inter.list[[i]], vars.to.regress = c("percent.mito","nCount_RNA"))
     inter.list[[i]] <- FindVariableFeatures(inter.list[[i]],selection.method = "vst" ,nfeatures = 3000)
	 }

## label transfer
anchors <- FindTransferAnchors(reference = inter.list[[1]], query = inter.list[[2]], 
	       k.anchor = 10, max.features = 1000, dims=1:30)
predictions <- TransferData(anchorset = anchors, refdata = inter.list[[1]]@active.ident, dims=1:30)
result <- data.frame(predictions[c("predicted.id","prediction.score.max")])
result <- cbind(inter.list[[2]]@meta.data,result)
result$predicted.id <- factor(result$predicted.id,levels = levels(reference.obj@active.ident))
inter.list[[2]]@meta.data$cellGroup <- result$predicted.id
inter.list[[2]]@active.ident <- factor(inter.list[[2]]@meta.data$cellGroup, levels = levels(reference.obj@active.ident))

XieXW.clean.obj <- inter.list[[2]]

## save data
save(XieXW.clean.obj, inter.list ,file = paste0(pwd,"/int/02.Celllabel.transfer.XieXW.Rdata"))
