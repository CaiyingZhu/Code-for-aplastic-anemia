rm(list=ls())
## set output directory of cellphoneDB as working directory
setwd("/public/home/zhucy/SingleCell/AA/cellphoneDB/")
library(limma)
library(ggplot2)
library(ggpubr)
pwd <- getwd()
dir <- "Tissue"

## read in expr and pvalue matrix
exprMat <- read.table(paste0(dir,"/means.txt"),row.names=2,sep="\t",head=T)
pMat <- read.table(paste0(dir,"/pvalues.txt"),row.names=2,sep="\t",head=T)
colnames(exprMat) <- gsub("HSC.MPP","HSCMPP",colnames(exprMat))
colnames(exprMat) <- gsub(".BM.",".BM_",fixed = TRUE,colnames(exprMat))
colnames(exprMat) <- gsub(".PB.",".PB_",fixed = TRUE,colnames(exprMat))

colnames(pMat) <- gsub("HSC.MPP","HSCMPP",colnames(pMat))
colnames(pMat) <- gsub(".BM.",".BM_",fixed = TRUE,colnames(pMat))
colnames(pMat) <- gsub(".PB.",".PB_",fixed = TRUE,colnames(pMat))

## prepared meta table
meta <- data.frame("group"=colnames(exprMat)[-c(1:10)],row.names=colnames(exprMat)[-c(1:10)])

## split group pair
tab <- data.frame(strsplit2(meta$group,'_',fixed=TRUE))
colnames(tab) <- c("source","target")
tab$source <- as.character(tab$source)
tab$target <- as.character(tab$target)

meta <- cbind(meta,tab)
for(i in colnames(tab)){
     tab <- data.frame(strsplit2(meta[[i]],'.',fixed=TRUE))
	 colnames(tab) <- paste0(i,c(".phenotype",".cellGroup",".Disease_type",".Tissue"))
     meta <- cbind(meta,tab)
	 }	 
meta$source.phenotype <- as.character(meta$source.phenotype)
meta$source.cellGroup <- as.character(meta$source.cellGroup)
meta$source.Disease_type <- as.character(meta$source.Disease_type)
meta$source.Tissue <- as.character(meta$source.Tissue)
meta$target.phenotype <- as.character(meta$target.phenotype)
meta$target.cellGroup <- as.character(meta$target.cellGroup)
meta$target.Disease_type <- as.character(meta$target.Disease_type)
meta$target.Tissue <- as.character(meta$target.Tissue)
	 
#count significant pair number in cluster-cluster
{
CountSig <- function(pMat){
     data <- pMat
	 data <- data[,-c(1:10)]
	 CountMat <- matrix(NA,nrow=ncol(data),ncol=1)
     row.names(CountMat) <- colnames(data)
	 colnames(CountMat) <- "count"
	 for(i in colnames(data)){
		 k=as.data.frame(table(data[,i]<=0.05))
		 CountMat[i,1]=k[2,2]
		     }
		 CountMat[is.na(CountMat)] <- 0			 
		 return(CountMat)
	 }
}

## counting significant pair number in cluster-cluster
Sigcount <- CountSig(pMat)
table(rownames(Sigcount)==rownames(meta))
meta <- cbind(meta,Sigcount)

## remove not significant molecular pairs
pMat <- pMat[,rownames(meta)]

pMat_conver <- pMat
pMat_conver[pMat_conver > 0.05] <- 1
pMat_conver <- pMat_conver[rowSums(pMat_conver) < dim(pMat_conver)[2],]

exprMat <- exprMat[rownames(pMat_conver),colnames(pMat_conver)]
pMat <- pMat[rownames(pMat_conver),colnames(pMat_conver)]

dim(pMat)        # 374 3844
dim(exprMat)     # 374 3844
dim(pMat_conver) # 374 3844

## retain T-HSPCs interaction pairs
Thspc <- meta[meta$source.phenotype %in% c("CD4T","CD8T") & meta$target.phenotype == "hspc",]
hspcT <- meta[meta$source.phenotype == "hspc"  & meta$target.phenotype %in% c("CD4T","CD8T"),]

pMat_filter.Thspc <- pMat_conver[,rownames(Thspc)]
pMat_filter.Thspc <- pMat_filter.Thspc[rowSums(pMat_filter.Thspc) < dim(pMat_filter.Thspc)[2],]

pMat_filter.hspcT <- pMat_conver[,rownames(hspcT)]
pMat_filter.hspcT <- pMat_filter.hspcT[rowSums(pMat_filter.hspcT) < dim(pMat_filter.hspcT)[2],]

dim(pMat_filter.Thspc) #  155 945
dim(pMat_filter.hspcT) #  197 945

pMat_filter.Thspc <- pMat[rownames(pMat_filter.Thspc),colnames(pMat_filter.Thspc)]
pMat_filter.hspcT <- pMat[rownames(pMat_filter.hspcT),colnames(pMat_filter.hspcT)]

expMat_filter.Thspc <- exprMat[rownames(pMat_filter.Thspc),colnames(pMat_filter.Thspc)]
expMat_filter.hspcT <- exprMat[rownames(pMat_filter.hspcT),colnames(pMat_filter.hspcT)]

## replace "_" as ":" to link molecules in a complex
## pvalue matrix
rownames(pMat_filter.Thspc) <- gsub("ACVR_","ACVR:",fixed=TRUE,rownames(pMat_filter.Thspc))
rownames(pMat_filter.Thspc) <- gsub("ACVR1_","ACVR1:",fixed=TRUE,rownames(pMat_filter.Thspc))
rownames(pMat_filter.Thspc) <- gsub("BMPR1B_","BMPR1B:",fixed=TRUE,rownames(pMat_filter.Thspc))
rownames(pMat_filter.Thspc) <- gsub("BMR1B_","BMR1B:",fixed=TRUE,rownames(pMat_filter.Thspc))
rownames(pMat_filter.hspcT) <- gsub("_TGFR",":TGFR",fixed=TRUE,rownames(pMat_filter.hspcT))
rownames(pMat_filter.hspcT) <- gsub("PlexinA1_","PlexinA1:",fixed=TRUE,rownames(pMat_filter.hspcT))
rownames(pMat_filter.hspcT) <- gsub("PlexinA3_","PlexinA3:",fixed=TRUE,rownames(pMat_filter.hspcT))
rownames(pMat_filter.hspcT) <- gsub("PlexinA4_","PlexinA4:",fixed=TRUE,rownames(pMat_filter.hspcT))

## expression matrix
rownames(expMat_filter.Thspc) <- gsub("ACVR_","ACVR:",fixed=TRUE,rownames(expMat_filter.Thspc))
rownames(expMat_filter.Thspc) <- gsub("ACVR1_","ACVR1:",fixed=TRUE,rownames(expMat_filter.Thspc))
rownames(expMat_filter.Thspc) <- gsub("BMPR1B_","BMPR1B:",fixed=TRUE,rownames(expMat_filter.Thspc))
rownames(expMat_filter.Thspc) <- gsub("BMR1B_","BMR1B:",fixed=TRUE,rownames(expMat_filter.Thspc))
rownames(expMat_filter.hspcT) <- gsub("_TGFR",":TGFR",fixed=TRUE,rownames(expMat_filter.hspcT))
rownames(expMat_filter.hspcT) <- gsub("PlexinA1_","PlexinA1:",fixed=TRUE,rownames(expMat_filter.hspcT))
rownames(expMat_filter.hspcT) <- gsub("PlexinA3_","PlexinA3:",fixed=TRUE,rownames(expMat_filter.hspcT))
rownames(expMat_filter.hspcT) <- gsub("PlexinA4_","PlexinA4:",fixed=TRUE,rownames(expMat_filter.hspcT))

## add molecular pair to pair column
pMat_filter.Thspc$pair <- rownames(pMat_filter.Thspc)
expMat_filter.Thspc$pair <- rownames(expMat_filter.Thspc)

## change colnames HSPCs-T cells to T-HSPCs
tab <- data.frame(strsplit2(colnames(pMat_filter.hspcT),'_',fixed=TRUE))
colnames(pMat_filter.hspcT) <- paste(tab$X2,tab$X1,sep="_")

tab <- data.frame(strsplit2(colnames(expMat_filter.hspcT),'_',fixed=TRUE))
colnames(expMat_filter.hspcT) <- paste(tab$X2,tab$X1,sep="_")

## change molecular pair name order from HSPC-T to T-HSPCs interaction, add pair column
tab <- data.frame(strsplit2(rownames(pMat_filter.hspcT),'_',fixed=TRUE))
pMat_filter.hspcT$pair <- paste(tab$X2,tab$X1,sep="_")
pMat_filter.hspcT <- pMat_filter.hspcT[,colnames(pMat_filter.Thspc)]

tab <- data.frame(strsplit2(rownames(expMat_filter.hspcT),'_',fixed=TRUE))
expMat_filter.hspcT$pair <- paste(tab$X2,tab$X1,sep="_")
expMat_filter.hspcT <- expMat_filter.hspcT[,colnames(expMat_filter.Thspc)]

## meta table order
source <- hspcT[,c("target","target.phenotype","target.cellGroup","target.Disease_type","target.Tissue")]
target <- hspcT[,c("source","source.phenotype","source.cellGroup","source.Disease_type","source.Tissue")]

hspcT[,c("target","target.phenotype","target.cellGroup","target.Disease_type","target.Tissue")] <- target
hspcT[,c("source","source.phenotype","source.cellGroup","source.Disease_type","source.Tissue")] <- source

table(rownames(hspcT)==hspcT$group)
tab <- data.frame(strsplit2(rownames(hspcT),'_',fixed=TRUE))
hspcT$group <- paste(tab$X2,tab$X1,sep="_")
rownames(hspcT) <- paste(tab$X2,tab$X1,sep="_")
hspcT <- hspcT[rownames(Thspc),]

## combine T cells-HSPCs and HSPCs-T cells interaction pairs
table(colnames(pMat_filter.hspcT)==colnames(pMat_filter.Thspc))
table(colnames(expMat_filter.hspcT)==colnames(expMat_filter.Thspc))

pMat_filter <- rbind(pMat_filter.Thspc, pMat_filter.hspcT)
expMat_filter <- rbind(expMat_filter.Thspc, expMat_filter.hspcT)

## remove duplicated pairs
pMat_filter <- pMat_filter[!duplicated(pMat_filter$pair),]
expMat_filter <- expMat_filter[!duplicated(expMat_filter$pair),]

Thspc$interaction <- hspcT$count + Thspc$count
Thspc$count <- NULL

## when separete bone marrow and peripheral blood cells, there were some cell clusters 
## with less than 10 cells in SAA or non-SAA, we consider to remove those clusters.
## remove cellGroup with less than 10 cells (CD4T Effector, Neu1 and EBM)
Thspc$phenotype.cellGroup <- paste(Thspc$source.phenotype,Thspc$source.cellGroup,sep=".")
Thspc_rm <- Thspc[ Thspc$phenotype.cellGroup != "CD4T.Effector",]

source.cellGroup.order <- c("Naive","Memory","Effector")
target.cellGroup.order <- c("HSCMPP","LMPP","MEP","MLP","EBM","Neu1","Neu2","MD1","MD2")
Thspc_rm$source.cellGroup <- factor(Thspc_rm$source.cellGroup,levels = source.cellGroup.order,order=T)
Thspc_rm$target.cellGroup <- factor(Thspc_rm$target.cellGroup,levels = target.cellGroup.order,order=T)
Thspc_rm$group <- as.character(Thspc_rm$group)

## remain the interaction within each disease type
type <- c("HD","CAA","SAA")
meta_filter <- c()
for(i in type){
     this.meta <- Thspc_rm[Thspc_rm$source.Disease_type == i & Thspc_rm$target.Disease_type == i,]
     meta_filter <- rbind(meta_filter,this.meta)
	 }
meta_filter$source.Disease_type <- factor(meta_filter$source.Disease_type,levels=type,order=T)

## interaction distribution
ggplot(meta_filter,aes(x=interaction)) + geom_line(stat="density",colour="black")
ggsave(filename =paste0(dir,"/output/interaction_distribution.pdf"),height=3,width=6,onefile=F)

## remained non-SAA and control
nonSAACtrl.meta <- meta_filter[meta_filter$source.Disease_type != "SAA",]
nonSAACtrl.meta$Tissue.Disease_type <- paste(nonSAACtrl.meta$source.Tissue,nonSAACtrl.meta$source.Disease_type,sep=".")
nonSAACtrl.meta$Tissue.Disease_type <- factor(nonSAACtrl.meta$Tissue.Disease_type,levels=c("BM.CAA","BM.HD","PB.CAA","PB.HD"))
nonSAACtrl.meta$target.Disease_type <- factor(nonSAACtrl.meta$target.Disease_type,levels=c("CAA","HD"))
my_comparisons <- list(c("BM.CAA", "BM.HD"), c("PB.CAA", "PB.HD"), c("BM.CAA", "PB.CAA"))
## boxplot 
p <- ggboxplot(nonSAACtrl.meta, x = "Tissue.Disease_type", 
     y = "interaction",
	 fill="target.Disease_type",
     color = "black",
	 outlier.shape = NA,
	 palette = c("FireBrick","DarkGray"),
     facet.by = "source.phenotype",
	 font.label = list(size = 8, color = "black"),
	 title="Number of Interaction Pairs in non-SAA and Ctrl",
     xlab="",
	 order=c("BM.HD","BM.CAA","PB.HD","PB.CAA"),
	 ylab="Number of Interaction Pairs",
	 short.panel.labs = TRUE)
     # Use only p.format as label. Remove method name.
     p + stat_compare_means(comparisons = my_comparisons,label.y = c(37, 37, 40), method="t.test",aes(label = paste0("p = ", ..p.format..))) +
	 theme(axis.text=element_text(size=6,color="black"),axis.text.x=element_text(angle=45, hjust=1,vjust=1))
     ggsave(filename =paste0(dir,"/output/Interaction number boxplot in BM vs PB non-SAA.unpaired.pdf"),height=4.2,width=4,onefile=F)

## unique interaction in each disease type
pMat_filter_saa <- pMat_filter[,c(rownames(meta_filter[meta_filter$source.Disease_type == "SAA",]),"pair")]
pMat_filter_caa <- pMat_filter[,c(rownames(meta_filter[meta_filter$source.Disease_type == "CAA",]),"pair")]
pMat_filter_ctrl <- pMat_filter[,c(rownames(meta_filter[meta_filter$source.Disease_type == "HD",]),"pair")]

## retain significant pairs
## non-SAA
pMat_filter_caa_conver <- pMat_filter_caa
pMat_filter_caa_conver[1:90][pMat_filter_caa_conver[1:90] > 0.05] <- 1
pMat_filter_caa_conver <- pMat_filter_caa_conver[rowSums(pMat_filter_caa_conver[1:90]) < dim(pMat_filter_caa_conver[1:90])[2],]
pMat_filter_caa <- pMat_filter_caa[rownames(pMat_filter_caa_conver),]
expMat_filter_caa <- expMat_filter[rownames(pMat_filter_caa_conver),colnames(pMat_filter_caa)]
dim(pMat_filter_caa)   # 102  91
dim(expMat_filter_caa) # 102  91

## Ctrl
pMat_filter_ctrl_conver <- pMat_filter_ctrl
pMat_filter_ctrl_conver[1:90][pMat_filter_ctrl_conver[1:90] > 0.05] <- 1
pMat_filter_ctrl_conver <- pMat_filter_ctrl_conver[rowSums(pMat_filter_ctrl_conver[1:90]) < dim(pMat_filter_ctrl_conver[1:90])[2],]
pMat_filter_ctrl <- pMat_filter_ctrl[rownames(pMat_filter_ctrl_conver),]
expMat_filter_ctrl <- expMat_filter[rownames(pMat_filter_ctrl_conver),colnames(pMat_filter_ctrl)]
dim(pMat_filter_ctrl)   # 67 91
dim(expMat_filter_ctrl) # 67 91

## melt convert pvalue and expression matrix 
library(data.table)
melt.conver <- function(pMat,expMat){
     pMelt <- reshape2::melt(pMat)
	 expMelt <- reshape2::melt(expMat)
	 colnames(pMelt) <- c("pair","group","pvalue")
	 colnames(expMelt) <- c("pair","group","expr")
	 Melt_pval_exp <- cbind(pMelt,expr=expMelt[["expr"]])
     tab <- data.frame(strsplit2(Melt_pval_exp$group,'_',fixed=TRUE))
	 colnames(tab) <- c("source","target")
     tab$source <- as.character(tab$source)
	 tab$target <- as.character(tab$target)
     Melt_pval_exp <- cbind(Melt_pval_exp,tab)
	 for(i in colnames(tab)){
	     tab_split <- data.frame(strsplit2(Melt_pval_exp[[i]],'.',fixed=TRUE))
		 colnames(tab_split) <- paste0(i,c(".phenotype",".cellGroup",".Disease_type",".Tissue"))
		 Melt_pval_exp <- cbind(Melt_pval_exp,tab_split)
		 }
	 Melt_pval_exp$phenotype.cellGroup <- paste(Melt_pval_exp$source.phenotype,Melt_pval_exp$source.cellGroup,sep=".")
	 phenotype.order <- c("CD4T.Naive","CD4T.Memory","CD8T.Naive","CD8T.Memory","CD8T.Effector")
     Melt_pval_exp$phenotype.cellGroup <- factor(Melt_pval_exp$phenotype.cellGroup,levels=phenotype.order,order=T)
	 Melt_pval_exp$source.cellGroup <- factor(Melt_pval_exp$source.cellGroup,levels=source.cellGroup.order,order=T)
	 Melt_pval_exp$target.cellGroup <- factor(Melt_pval_exp$target.cellGroup,levels=target.cellGroup.order,order=T)
	 Melt_pval_exp$source.phenotype <- factor(Melt_pval_exp$source.phenotype,levels=c("CD4T","CD8T"),order=T)
	 return(Melt_pval_exp)
	 }

Melt_caa_sig <- melt.conver(pMat=pMat_filter_caa, expMat=expMat_filter_caa)
Melt_ctrl_sig <- melt.conver(pMat=pMat_filter_ctrl, expMat=expMat_filter_ctrl)

## cell type specific pairs
library(plyr)
library(ggplot2)
phenotype.order <- c("CD4T.Naive","CD4T.Memory","CD8T.Naive","CD8T.Memory","CD8T.Effector")
caa_uniq <- c()
for(i in phenotype.order){
     caa_bm_sig <- Melt_caa_sig[Melt_caa_sig$phenotype.cellGroup==i & Melt_caa_sig$source.Tissue=="BM" ,]
	 caa_pb_sig <- Melt_caa_sig[Melt_caa_sig$phenotype.cellGroup==i & Melt_caa_sig$source.Tissue=="PB" ,]
     caa_sig_rbind <- rbind(caa_bm_sig,caa_pb_sig)
	 ctrl_bm_sig <- Melt_ctrl_sig[Melt_ctrl_sig$phenotype.cellGroup==i & Melt_ctrl_sig$source.Tissue=="BM",]
	 ctrl_pb_sig <- Melt_ctrl_sig[Melt_ctrl_sig$phenotype.cellGroup==i & Melt_ctrl_sig$source.Tissue=="PB",]
     caa_bm_sig_pval <- caa_bm_sig[caa_bm_sig$pvalue <= 0.05,]
	 caa_pb_sig_pval <- caa_pb_sig[caa_pb_sig$pvalue <= 0.05,]
     ctrl_bm_sig_pval <- ctrl_bm_sig[ctrl_bm_sig$pvalue <= 0.05,]
	 ctrl_pb_sig_pval <- ctrl_pb_sig[ctrl_pb_sig$pvalue <= 0.05,]
	 caa_bm_sig_pval_uniq <- caa_bm_sig_pval[!caa_bm_sig_pval$pair %in% c(ctrl_bm_sig_pval$pair,ctrl_pb_sig_pval$pair),]
     caa_pb_sig_pval_uniq <- caa_pb_sig_pval[!caa_pb_sig_pval$pair %in% c(ctrl_bm_sig_pval$pair,ctrl_pb_sig_pval$pair),]	
     caa_sig_pval_uniq <- rbind(caa_bm_sig_pval_uniq,caa_pb_sig_pval_uniq) 
	 caa_uniq <- rbind(caa_uniq,caa_sig_pval_uniq)
	 count.data <- ddply(caa_sig_pval_uniq, .(pair), summarize, counts= length(pair))
	 count.data <- count.data[order(count.data$counts),]
     pairs <- unique(count.data$pair)
	 caa_sig <- caa_sig_rbind[caa_sig_rbind$pair %in% pairs,]
	 caa_sig$pair <- factor(caa_sig$pair,levels=pairs,order=T)
	 caa_sig$pvalue[caa_sig$pvalue > 0.05] <- 0.1
	 caa_sig$pvalue[caa_sig$pvalue <= 0.001] <- 0.001
	 caa_sig$expr[caa_sig$expr >= 2.0] <- 2.0
	 if(T){
	 ggplot(caa_sig,aes(x=factor(target.cellGroup,levels=target.cellGroup.order),y=pair))+
         geom_point(aes(color=expr,size=-log10(pvalue)),alpha=0.8) +
		 scale_size_continuous(range = c(0.2,1.2)) + 
      	 scale_colour_gradient(low="#0000CD",high="#FF0000",limits=c(0.0,2.0)) +
         xlab("") + ylab("") +
         theme_bw()+theme(axis.text=element_text(size=7,color="black"),axis.text.x=element_text(angle=45, hjust=1,vjust=1),
     	 axis.text.y = element_text(size = 7,color="black"),axis.title.y = element_text(size = 7),
         panel.grid.minor = element_blank())+labs(title=i) +
		 facet_grid(~source.Tissue, scales = "free",space="free")
         ggsave(filename =paste0(dir,"/output/dotplot_significant.pair.",i,".non-SAA.pdf"),height=5.6,width=4,onefile=F)
	     }
		 if(T){
	     ggplot(caa_sig,aes(x=factor(target.cellGroup,levels=target.cellGroup.order),y=pair))+
       	     geom_point(aes(color=expr,size=-log10(pvalue)),alpha=0.8) +
			 scale_size_continuous(range = c(0.2,1.2)) + 
      	     scale_colour_gradient(low="#0000CD",high="#FF0000",limits=c(0.0,2.0)) +
         	 xlab("") + ylab("") +
             theme_bw()+ 
			 facet_grid(~source.Tissue, scales = "free",space="free") +
			 theme(axis.text=element_text(size=7,color="black"),axis.text.x=element_text(angle=45, hjust=1,vjust=1),
     	     axis.text.y = element_blank(),axis.title.y = element_blank(),
             panel.grid.minor = element_blank())+labs(title=i)+ theme(legend.position='none')			 
             ggsave(filename =paste0(dir,"/output/dotplot_significant.pair.",i,".non-SAA.NoLegend.pdf"),height=5.6,width=1.4,onefile=F)
	         }
	     }

## plot in one figure (unique pairs in BM and PB)
count.data <- ddply(caa_uniq, .(source.Tissue,pair), summarize, counts= length(pair))
count.data <- count.data[order(count.data$counts),]
pairs <- unique(count.data$pair) 
table(count.data$source.Tissue)
#BM PB
#45 51
CAA_sig_bm <- Melt_caa_sig[Melt_caa_sig$source.Tissue=="BM" ,]
CAA_sig_pb <- Melt_caa_sig[Melt_caa_sig$source.Tissue=="PB" ,]
CAA_sig_bm <- CAA_sig_bm[CAA_sig_bm$pvalue <= 0.05,]
CAA_sig_pb <- CAA_sig_pb[CAA_sig_pb$pvalue <= 0.05,]

over_pair <- intersect(as.vector(count.data[count.data$source.Tissue=="BM",]$pair),as.vector(count.data[count.data$source.Tissue=="PB",]$pair))
over_pair <- pairs[pairs %in% over_pair]
bm_uniq_pair <- setdiff(as.vector(count.data[count.data$source.Tissue=="BM",]$pair),as.vector(CAA_sig_pb$pair))
bm_uniq_pair <- pairs[pairs %in% bm_uniq_pair]
pb_uniq_pair <- setdiff(as.vector(count.data[count.data$source.Tissue=="PB",]$pair),as.vector(CAA_sig_bm$pair))
pb_uniq_pair <- pairs[pairs %in% pb_uniq_pair]
pairs_order <- c(over_pair,bm_uniq_pair,pb_uniq_pair)
Melt_sig_pair <- Melt_caa_sig[Melt_caa_sig$pair %in% pairs_order,]
Melt_sig_pair$pvalue[Melt_sig_pair$pvalue > 0.05] <- 0.1
Melt_sig_pair$pvalue[Melt_sig_pair$pvalue <= 0.001] <- 0.001
Melt_sig_pair$tag <- ""
Melt_sig_pair$tag[Melt_sig_pair$pair %in% over_pair] <- "overlap"
Melt_sig_pair$tag[Melt_sig_pair$pair %in% bm_uniq_pair] <- "BM"
Melt_sig_pair$tag[Melt_sig_pair$pair %in% pb_uniq_pair] <- "PB"
Melt_sig_pair$tag <- factor(Melt_sig_pair$tag,levels=c("overlap","BM","PB"))
Melt_sig_pair$source <- gsub(".CAA.",".",fixed=TRUE,Melt_sig_pair$source)
Melt_sig_pair$source <- factor(Melt_sig_pair$source,levels=c("CD4T.Naive.PB","CD4T.Memory.PB","CD8T.Naive.PB","CD8T.Memory.PB","CD8T.Effector.PB",
                         "CD4T.Naive.BM","CD4T.Memory.BM","CD8T.Naive.BM","CD8T.Memory.BM","CD8T.Effector.BM"))
Melt_sig_pair$expr[Melt_sig_pair$expr >= 2.0] <- 2.0
ggplot(Melt_sig_pair,aes(x=factor(target.cellGroup,levels=target.cellGroup.order),y=pair))+
     geom_point(aes(color=expr,size=-log10(pvalue)),alpha=0.8) +
	 scale_size_continuous(range = c(0.2,1.2)) + 
     scale_colour_gradient(low="#0000CD",high="#FF0000",limits=c(0.0,2.0)) +
     xlab("") + ylab("") +
     theme_bw()+theme(axis.text=element_text(size=7,color="black"),axis.text.x=element_text(angle=90, hjust=1,vjust=1),
     axis.text.y = element_text(size = 7,color="black"),axis.title.y = element_text(size = 7),
     panel.grid.minor = element_blank()) + 
	 facet_grid(tag~source, scales = "free",space="free")
     ggsave(filename =paste0(dir,"/output/dotplot_significant.pair.total.pair.BMvsPB.non-SAA.pdf"),height=6.2,width=10,onefile=F)

ggplot(Melt_sig_pair,aes(x=factor(target.cellGroup,levels=target.cellGroup.order),y=pair))+
     geom_point(aes(color=expr,size=-log10(pvalue)),alpha=0.8) +
	 scale_size_continuous(range = c(0.2,1.2)) + 
     scale_colour_gradient(low="#0000CD",high="#FF0000",limits=c(0.0,2.0)) +
  	 xlab("") + ylab("") +
     theme_bw()+ 
	 facet_grid(tag~source, scales = "free",space="free") +
	 theme(axis.text=element_text(size=7,color="black"),axis.text.x=element_text(angle=45, hjust=1,vjust=1),
     axis.text.y = element_blank(),axis.title.y = element_blank(),
     panel.grid.minor = element_blank()) + theme(legend.position='none')			 
     ggsave(filename =paste0(dir,"/output/dotplot_significant.pair.total.pair.BMvsPB.non-SAA.NoLegend.pdf"),height=6.2,width=7,onefile=F)
	    
## VennDiagram
library(VennDiagram)
Melt_caa_sig_pval <- Melt_caa_sig[Melt_caa_sig$pvalue <= 0.05,]
Melt_ctrl_sig_pval <- Melt_ctrl_sig[Melt_ctrl_sig$pvalue <= 0.05,]
for(i in c("CD4T","CD8T")){
     pdf(file=paste0(dir,"/output/CAAvsCtrl_venn diagram.BAvsPB.",i,".pdf"),4,4,onefile=F)
     p <- venn.diagram(list(BMAA=Melt_caa_sig_pval[Melt_caa_sig_pval$source.phenotype==i & Melt_caa_sig_pval$source.Tissue=="BM",]$pair,
	     BMCtrl=Melt_ctrl_sig_pval[Melt_ctrl_sig_pval$source.phenotype==i & Melt_ctrl_sig_pval$source.Tissue=="BM",]$pair,
		 PBAA=Melt_caa_sig_pval[Melt_caa_sig_pval$source.phenotype==i & Melt_caa_sig_pval$source.Tissue=="PB",]$pair,
		 PBCtrl=Melt_ctrl_sig_pval[Melt_ctrl_sig_pval$source.phenotype==i & Melt_ctrl_sig_pval$source.Tissue=="PB",]$pair),
         filename= NULL,
         lwd=1,
         lty=3,
         main = i,
         main.cex =1.2,
         cex = 2,
         cat.cex = 1,
         col=c('steelblue1','darkorange1','#98FB98','#9370DB'),
         fill=c('steelblue1','darkorange1','#98FB98','#9370DB'),
         cat.col="black")
         grid.draw(p)
         dev.off()
	 }

## SAA vs Ctrl
## remain SAA and Ctrl
## remove HSPCs cell type with less than 10 cells
meta_filter_rm <- meta_filter[!meta_filter$target.cellGroup %in% c("EBM","Neu1"),]
meta_filter_rm_caa <- meta_filter_rm[meta_filter_rm$source.Disease_type != "CAA",]

meta_filter_rm_caa$Tissue.Disease_type <- paste(meta_filter_rm_caa$source.Tissue,meta_filter_rm_caa$source.Disease_type,sep=".")
meta_filter_rm_caa$Tissue.Disease_type <- factor(meta_filter_rm_caa$Tissue.Disease_type,levels=c("BM.SAA","BM.HD","PB.SAA","PB.HD"))
meta_filter_rm_caa$target.Disease_type <- factor(meta_filter_rm_caa$target.Disease_type,levels=c("SAA","HD"))
my_comparisons <- list( c("BM.SAA", "BM.HD"), c("PB.SAA", "PB.HD"), c("BM.SAA", "PB.SAA"))

## boxplot 
p <- ggboxplot(meta_filter_rm_caa, x = "Tissue.Disease_type", 
     y = "interaction",
	 fill="target.Disease_type",
     color = "black",
	 outlier.shape=NA,
	 palette = c("FireBrick","DarkGray"),
     facet.by = "source.phenotype",
	 font.label = list(size = 8, color = "black"),
	 title="Number of Interaction Pairs in SAA and Ctrl",
     xlab="",
	 order=c("BM.HD","BM.SAA","PB.HD","PB.SAA"),
	 ylab="Number of Interaction Pairs",
	 short.panel.labs = TRUE)
     # Use only p.format as label. Remove method name.
     p + stat_compare_means(comparisons = my_comparisons,label.y = c(40, 40, 44), method="t.test",aes(label = paste0("p = ", ..p.format..))) +
	 theme(axis.text=element_text(size=8,color="black"),axis.text.x=element_text(angle=45, hjust=1,vjust=1))
     ggsave(filename =paste0(dir,"/output/Interaction number boxplot in BM vs PB SAA.pdf"),height=4.4,width=4,onefile=F)

pMat_filter_saa <- pMat_filter[,c(rownames(meta_filter_rm_caa[meta_filter_rm_caa$source.Disease_type == "SAA",]),"pair")]
pMat_filter_ctrl <- pMat_filter[,c(rownames(meta_filter_rm_caa[meta_filter_rm_caa$source.Disease_type == "HD",]),"pair")]

## SAA
pMat_filter_saa_conver <- pMat_filter_saa
pMat_filter_saa_conver[1:70][pMat_filter_saa_conver[1:70] > 0.05] <- 1
pMat_filter_saa_conver <- pMat_filter_saa_conver[rowSums(pMat_filter_saa_conver[1:70]) < dim(pMat_filter_saa_conver[1:70])[2],]
pMat_filter_saa <- pMat_filter_saa[rownames(pMat_filter_saa_conver),]
expMat_filter_saa <- expMat_filter[rownames(pMat_filter_saa_conver),colnames(pMat_filter_saa)]
dim(pMat_filter_saa)   # 127  71
dim(expMat_filter_saa) # 127  71

## Ctrl
pMat_filter_ctrl_conver <- pMat_filter_ctrl
pMat_filter_ctrl_conver[1:70][pMat_filter_ctrl_conver[1:70] > 0.05] <- 1
pMat_filter_ctrl_conver <- pMat_filter_ctrl_conver[rowSums(pMat_filter_ctrl_conver[1:70]) < dim(pMat_filter_ctrl_conver[1:70])[2],]
pMat_filter_ctrl <- pMat_filter_ctrl[rownames(pMat_filter_ctrl_conver),]
expMat_filter_ctrl <- expMat_filter[rownames(pMat_filter_ctrl_conver),colnames(pMat_filter_ctrl)]
dim(pMat_filter_ctrl)   # 62 71
dim(expMat_filter_ctrl) # 62 71

## melt convert
Melt_saa_sig <- melt.conver(pMat=pMat_filter_saa, expMat=expMat_filter_saa)
Melt_ctrl_sig <- melt.conver(pMat=pMat_filter_ctrl, expMat=expMat_filter_ctrl)

## cell type specific pairs
library(plyr)
library(ggplot2)
saa_uniq <- c()
for(i in phenotype.order){
     saa_bm_sig <- Melt_saa_sig[Melt_saa_sig$phenotype.cellGroup==i & Melt_saa_sig$source.Tissue=="BM" ,]
	 saa_pb_sig <- Melt_saa_sig[Melt_saa_sig$phenotype.cellGroup==i & Melt_saa_sig$source.Tissue=="PB" ,]
     saa_sig_rbind <- rbind(saa_bm_sig,saa_pb_sig)
	 ctrl_bm_sig <- Melt_ctrl_sig[Melt_ctrl_sig$phenotype.cellGroup==i & Melt_ctrl_sig$source.Tissue=="BM",]
	 ctrl_pb_sig <- Melt_ctrl_sig[Melt_ctrl_sig$phenotype.cellGroup==i & Melt_ctrl_sig$source.Tissue=="PB",]
     saa_bm_sig_pval <- saa_bm_sig[saa_bm_sig$pvalue <= 0.05,]
	 saa_pb_sig_pval <- saa_pb_sig[saa_pb_sig$pvalue <= 0.05,]
     ctrl_bm_sig_pval <- ctrl_bm_sig[ctrl_bm_sig$pvalue <= 0.05,]
	 ctrl_pb_sig_pval <- ctrl_pb_sig[ctrl_pb_sig$pvalue <= 0.05,]
	 saa_bm_sig_pval_uniq <- saa_bm_sig_pval[!saa_bm_sig_pval$pair %in% c(ctrl_bm_sig_pval$pair,ctrl_pb_sig_pval$pair),]
     saa_pb_sig_pval_uniq <- saa_pb_sig_pval[!saa_pb_sig_pval$pair %in% c(ctrl_pb_sig_pval$pair,ctrl_bm_sig_pval$pair),]	
     saa_sig_pval_uniq <- rbind(saa_bm_sig_pval_uniq,saa_pb_sig_pval_uniq) 
	 saa_uniq <- rbind(saa_uniq,saa_sig_pval_uniq)
	 count.data <- ddply(saa_sig_pval_uniq, .(pair), summarize, counts= length(pair))
	 count.data <- count.data[order(count.data$counts),]
     pairs <- unique(count.data$pair)
	 saa_sig <- saa_sig_rbind[saa_sig_rbind$pair %in% pairs,]
	 saa_sig$pair <- factor(saa_sig$pair,levels=pairs,order=T)
	 saa_sig$pvalue[saa_sig$pvalue > 0.05] <- 0.1
	 saa_sig$pvalue[saa_sig$pvalue <= 0.001] <- 0.001
	 saa_sig$expr[saa_sig$expr >= 2.0] <- 2.0
	 if(T){
	 ggplot(saa_sig,aes(x=factor(target.cellGroup,levels=target.cellGroup.order),y=pair))+
         geom_point(aes(color=expr,size=-log10(pvalue)),alpha=0.8) +
		 scale_size_continuous(range = c(0.2,1.2)) + 
      	 scale_colour_gradient(low="#0000CD",high="#FF0000") +
         xlab("") + ylab("") +
         theme_bw()+theme(axis.text=element_text(size=5,color="black"),axis.text.x=element_text(angle=45, hjust=1,vjust=1),
     	 axis.text.y = element_text(size = 5,color="black"),axis.title.y = element_text(size = 5),
         panel.grid.minor = element_blank())+labs(title=i) +
		 facet_grid(~source.Tissue, scales = "free",space="free")
         ggsave(filename =paste0(dir,"/output/dotplot_significant.pair.",i,".SAA.pdf"),height=6.4,width=6,onefile=F)
	     }
		 if(T){
	     ggplot(saa_sig,aes(x=factor(target.cellGroup,levels=target.cellGroup.order),y=pair))+
       	     geom_point(aes(color=expr,size=-log10(pvalue)),alpha=0.8) +
			 scale_size_continuous(range = c(0.2,1.2)) + 
      	     scale_colour_gradient(low="#0000CD",high="#FF0000",limits=c(0.0,2.0)) +
         	 xlab("") + ylab("") +
             theme_bw()+ 
			 facet_grid(~source.Tissue, scales = "free",space="free") +
			 theme(axis.text=element_text(size=5,color="black"),axis.text.x=element_text(angle=45, hjust=1,vjust=1),
     	     axis.text.y = element_blank(),axis.title.y = element_blank(),
             panel.grid.minor = element_blank())+labs(title=i)+ theme(legend.position='none')			 
             ggsave(filename =paste0(dir,"/output/dotplot_significant.pair.",i,".SAA.NoLegend.pdf"),height=6.4,width=1.2,onefile=F)
	         }
         }
 
## plot in one figure (unique pairs in BM and PB)
count.data <- ddply(saa_uniq, .(source.Tissue,pair), summarize, counts= length(pair))
count.data <- count.data[order(count.data$counts),]
pairs <- unique(count.data$pair) 
table(count.data$source.Tissue)
#BM PB
#68 71

Melt_saa_sig_pval <- Melt_saa_sig[Melt_saa_sig$pvalue <= 0.05,]
Melt_saa_sig_pval_bm <- Melt_saa_sig_pval[Melt_saa_sig_pval$source.Tissue == "BM",]
Melt_saa_sig_pval_pb <- Melt_saa_sig_pval[Melt_saa_sig_pval$source.Tissue == "PB",]
Melt_ctrl_sig_pval <- Melt_ctrl_sig[Melt_ctrl_sig$pvalue <= 0.05,]

over_pair <- intersect(as.vector(count.data[count.data$source.Tissue=="BM",]$pair),as.vector(count.data[count.data$source.Tissue=="PB",]$pair))
over_pair <- pairs[pairs %in% over_pair]
bm_uniq_pair <- setdiff(as.vector(count.data[count.data$source.Tissue=="BM",]$pair),as.vector(Melt_saa_sig_pval_pb$pair))
bm_uniq_pair <- pairs[pairs %in% bm_uniq_pair]
pb_uniq_pair <- setdiff(as.vector(count.data[count.data$source.Tissue=="PB",]$pair),as.vector(Melt_saa_sig_pval_bm$pair))
pb_uniq_pair <- pairs[pairs %in% pb_uniq_pair]
pairs_order <- c(over_pair,bm_uniq_pair,pb_uniq_pair)
Melt_sig_pair <- Melt_saa_sig[Melt_saa_sig$pair %in% pairs_order,]
Melt_sig_pair$pvalue[Melt_sig_pair$pvalue > 0.05] <- 0.1
Melt_sig_pair$pvalue[Melt_sig_pair$pvalue <= 0.001] <- 0.001
Melt_sig_pair$tag <- ""
Melt_sig_pair$tag[Melt_sig_pair$pair %in% over_pair] <- "overlap"
Melt_sig_pair$tag[Melt_sig_pair$pair %in% bm_uniq_pair] <- "BM"
Melt_sig_pair$tag[Melt_sig_pair$pair %in% pb_uniq_pair] <- "PB"
Melt_sig_pair$tag <- factor(Melt_sig_pair$tag,levels=c("overlap","BM","PB"))
Melt_sig_pair$source <- gsub(".SAA.",".",fixed=TRUE,Melt_sig_pair$source)
Melt_sig_pair$source <- factor(Melt_sig_pair$source,levels=c("CD4T.Naive.PB","CD4T.Memory.PB","CD8T.Naive.PB","CD8T.Memory.PB","CD8T.Effector.PB",
                         "CD4T.Naive.BM","CD4T.Memory.BM","CD8T.Naive.BM","CD8T.Memory.BM","CD8T.Effector.BM"))
target.cellGroup.order.v2 <- target.cellGroup.order[-c(5:6)]
Melt_sig_pair$expr[Melt_sig_pair$expr >= 2.0] <- 2.0
ggplot(Melt_sig_pair,aes(x=factor(target.cellGroup,levels=target.cellGroup.order.v2),y=pair))+
     geom_point(aes(color=expr,size=-log10(pvalue)),alpha=0.8) +
	 scale_size_continuous(range = c(0.2,1.2)) + 
     scale_colour_gradient(low="#0000CD",high="#FF0000",limits=c(0.0,2.0)) +
     xlab("") + ylab("") +
     theme_bw()+theme(axis.text=element_text(size=5,color="black"),axis.text.x=element_text(angle=45, hjust=1,vjust=1),
     axis.text.y = element_text(size = 5,color="black"),axis.title.y = element_text(size = 5),
     panel.grid.minor = element_blank()) + 
	 facet_grid(tag~source, scales = "free",space="free")
     ggsave(filename =paste0(dir,"/output/dotplot_significant.pair.total.pair.BMvsPB.SAA.pdf"),height=7,width=8.4,onefile=F)

ggplot(Melt_sig_pair,aes(x=factor(target.cellGroup,levels=target.cellGroup.order),y=pair))+
     geom_point(aes(color=expr,size=-log10(pvalue)),alpha=0.8) +
	 scale_size_continuous(range = c(0.2,1.2)) + 
     scale_colour_gradient(low="#0000CD",high="#FF0000",limits=c(0.0,2.0)) +
  	 xlab("") + ylab("") +
     theme_bw()+ 
	 facet_grid(tag~source, scales = "free",space="free") +
	 theme(axis.text=element_text(size=7,color="black"),axis.text.x=element_text(angle=45, hjust=1,vjust=1),
     axis.text.y = element_blank(),axis.title.y = element_blank(),
     panel.grid.minor = element_blank()) + theme(legend.position='none')			 
     ggsave(filename =paste0(dir,"/output/dotplot_significant.pair.total.pair.BMvsPB.SAA.NoLegend.pdf"),height=7.6,width=7,onefile=F)
	 
## VennDiagram
library(VennDiagram)
Melt_saa_sig_pval <- Melt_saa_sig[Melt_saa_sig$pvalue <= 0.05,]
Melt_ctrl_sig_pval <- Melt_ctrl_sig[Melt_ctrl_sig$pvalue <= 0.05,]
for(i in c("CD4T","CD8T")){
     pdf(file=paste0(dir,"/output/SAAvsCtrl_venn diagram.BAvsPB.",i,".pdf"),4,4,onefile=F)
     p <- venn.diagram(list(BMAA=Melt_saa_sig_pval[Melt_saa_sig_pval$source.phenotype==i & Melt_saa_sig_pval$source.Tissue=="BM",]$pair,
	     BMCtrl=Melt_ctrl_sig_pval[Melt_ctrl_sig_pval$source.phenotype==i & Melt_ctrl_sig_pval$source.Tissue=="BM",]$pair,
		 PBAA=Melt_saa_sig_pval[Melt_saa_sig_pval$source.phenotype==i & Melt_saa_sig_pval$source.Tissue=="PB",]$pair,
		 PBCtrl=Melt_ctrl_sig_pval[Melt_ctrl_sig_pval$source.phenotype==i & Melt_ctrl_sig_pval$source.Tissue=="PB",]$pair),
         filename= NULL,
         lwd=1,
         lty=3,
         main = i,
         main.cex =1.2,
         cex = 2,
         cat.cex = 1,
         col=c('steelblue1','darkorange1','#98FB98','#9370DB'),
         fill=c('steelblue1','darkorange1','#98FB98','#9370DB'),
         cat.col="black")
         grid.draw(p)
         dev.off()
	 }