## transcription factors (TFs) and genes co-expression networks 
## were inferred based on a subset of AA and control scRNA-seq
## expression data                                             
rm(list=ls())
pwd <- getwd()
## load libraries
.libPaths("/public/home/zhucy/software/R-3.5.0/UpdateLibrary")
suppressPackageStartupMessages({
	library(SCENIC)
	library(AUCell)
	library(RcisTarget)
	library(SingleCellExperiment)
     })

## load a subset of expression matrix
load("../int/hspc.aligned.data.Rdata")
# cell info
cellInfo <- hspc.aligned.obj@meta.data
# expression matrix
exprMat <- as.matrix(hspc.aligned.obj@data)

## Initialize settings
scenicOptions <- initializeScenic(org="hgnc", dbDir="resource", nCores=20)
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

## Co-expression network
runCorrelation(exprMat, scenicOptions)
runGenie3(exprMat, scenicOptions)

## Build Gene Regulatory Networks(GRNs)
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions)

## load all HSPCs expression matrix
load(../int/All.hspc.data.Rdata)
# cell info
cellInfo <- hspc.integrated.obj@meta.data
saveRDS(cellInfo, file="int/cellInfo.Rds")

# expression matrix
exprMat <- as.matrix(hspc.integrated.obj@data)

## Score regulon activity of each HSPC
runSCENIC_3_scoreCells(scenicOptions, exprMat)

### Binarize activity of regulons
runSCENIC_4_aucell_binarize(scenicOptions)