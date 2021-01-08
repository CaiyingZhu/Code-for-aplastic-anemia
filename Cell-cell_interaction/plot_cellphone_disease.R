rm(list=ls())
## set output directory of cellphoneDB as working directory
setwd("/public/home/zhucy/SingleCell/AA/cellphoneDB/")
library(limma)
library(ggplot2)
pwd <- getwd()
dir <- "disease"

## read in expression and pvalue matrix
## CAA represents non-SAA
## HD represents Ctrl
exprMat <- read.table(paste0(dir,"/means.txt"),row.names=2,sep="\t",head=T)
pMat <- read.table(paste0(dir,"/pvalues.txt"),row.names=2,sep="\t",head=T)
colnames(exprMat) <- gsub("HSC.MPP","HSCMPP",colnames(exprMat))
colnames(exprMat) <- gsub(".CAA.",".CAA_",colnames(exprMat))
colnames(exprMat) <- gsub(".SAA.",".SAA_",colnames(exprMat))
colnames(exprMat) <- gsub(".HD.",".HD_",colnames(exprMat))

colnames(pMat) <- gsub("HSC.MPP","HSCMPP",colnames(pMat))
colnames(pMat) <- gsub(".CAA.",".CAA_",colnames(pMat))
colnames(pMat) <- gsub(".SAA.",".SAA_",colnames(pMat))
colnames(pMat) <- gsub(".HD.",".HD_",colnames(pMat))

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
	 colnames(tab) <- paste0(i,c(".phenotype",".cellGroup",".Disease_type"))
     meta <- cbind(meta,tab)
	 }
	 
meta$source.phenotype <- as.character(meta$source.phenotype)
meta$source.cellGroup <- as.character(meta$source.cellGroup)
meta$source.Disease_type <- as.character(meta$source.Disease_type)
meta$target.phenotype <- as.character(meta$target.phenotype)
meta$target.cellGroup <- as.character(meta$target.cellGroup)
meta$target.Disease_type <- as.character(meta$target.Disease_type)
	 
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
		 CountMat[is.na(CountMat)] <- 0
		     }
		 return(CountMat)
	 }
}

#count significant pair number in cluster-cluster
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

dim(pMat)        # 288 2025
dim(exprMat)     # 288 2025
dim(pMat_conver) # 288 2025

## retain T-HSPCs interaction group pairs
## Tcell and HSPCs interaction
Thspc <- meta[meta$source.phenotype %in% c("CD4T","CD8T") & meta$target.phenotype == "hspc",]
hspcT <- meta[meta$source.phenotype == "hspc"  & meta$target.phenotype %in% c("CD4T","CD8T"),]

pMat_filter.Thspc <- pMat_conver[,rownames(Thspc)]
pMat_filter.Thspc <- pMat_filter.Thspc[rowSums(pMat_filter.Thspc) < dim(pMat_filter.Thspc)[2],]

pMat_filter.hspcT <- pMat_conver[,rownames(hspcT)]
pMat_filter.hspcT <- pMat_filter.hspcT[rowSums(pMat_filter.hspcT) < dim(pMat_filter.hspcT)[2],]

dim(pMat_filter.Thspc) # 103 486
dim(pMat_filter.hspcT) # 122 486

pMat_filter.Thspc <- pMat[rownames(pMat_filter.Thspc),colnames(pMat_filter.Thspc)]
pMat_filter.hspcT <- pMat[rownames(pMat_filter.hspcT),colnames(pMat_filter.hspcT)]

expMat_filter.Thspc <- exprMat[rownames(pMat_filter.Thspc),colnames(pMat_filter.Thspc)]
expMat_filter.hspcT <- exprMat[rownames(pMat_filter.hspcT),colnames(pMat_filter.hspcT)]

## replace "_" as ":" to link molecules in a complex
## pvalue matrix
rownames(pMat_filter.Thspc) <- gsub("ACVR_","ACVR:",rownames(pMat_filter.Thspc))
rownames(pMat_filter.Thspc) <- gsub("ACVR1_","ACVR1:",rownames(pMat_filter.Thspc))
rownames(pMat_filter.hspcT) <- gsub("_TGFR",":TGFR",rownames(pMat_filter.hspcT))

## expression matrix
rownames(expMat_filter.Thspc) <- gsub("ACVR_","ACVR:",rownames(pMat_filter.Thspc))
rownames(expMat_filter.Thspc) <- gsub("ACVR1_","ACVR1:",rownames(pMat_filter.Thspc))
rownames(expMat_filter.hspcT) <- gsub("_TGFR",":TGFR",rownames(pMat_filter.hspcT))

## add molecular pair to pair column
pMat_filter.Thspc$pair <- rownames(pMat_filter.Thspc)
expMat_filter.Thspc$pair <- rownames(expMat_filter.Thspc)

## change colnames HSPCs-T cells to T cells-HSPCs
tab <- data.frame(strsplit2(colnames(pMat_filter.hspcT),'_',fixed=TRUE))
colnames(pMat_filter.hspcT) <- paste(tab$X2,tab$X1,sep="_")

tab <- data.frame(strsplit2(colnames(expMat_filter.hspcT),'_',fixed=TRUE))
colnames(expMat_filter.hspcT) <- paste(tab$X2,tab$X1,sep="_")

## change molecular pair name order from HSPCs-T to 
## T-HSPCs interaction, and add pair column
tab <- data.frame(strsplit2(rownames(pMat_filter.hspcT),'_',fixed=TRUE))
pMat_filter.hspcT$pair <- paste(tab$X2,tab$X1,sep="_")
pMat_filter.hspcT <- pMat_filter.hspcT[,colnames(pMat_filter.Thspc)]

tab <- data.frame(strsplit2(rownames(expMat_filter.hspcT),'_',fixed=TRUE))
expMat_filter.hspcT$pair <- paste(tab$X2,tab$X1,sep="_")
expMat_filter.hspcT <- expMat_filter.hspcT[,colnames(expMat_filter.Thspc)]

## meta table order
source <- hspcT[,c("target","target.phenotype","target.cellGroup","target.Disease_type")]
target <- hspcT[,c("source","source.phenotype","source.cellGroup","source.Disease_type")]

hspcT[,c("target","target.phenotype","target.cellGroup","target.Disease_type")] <- target
hspcT[,c("source","source.phenotype","source.cellGroup","source.Disease_type")] <- source

table(rownames(hspcT)==hspcT$group)
tab <- data.frame(strsplit2(rownames(hspcT),'_',fixed=TRUE))
hspcT$group <- paste(tab$X2,tab$X1,sep="_")
rownames(hspcT) <- paste(tab$X2,tab$X1,sep="_")
hspcT <- hspcT[rownames(Thspc),]

## combined T-HSPCs and HSPCs-T interaction pairs
table(colnames(pMat_filter.hspcT)==colnames(pMat_filter.Thspc))
table(colnames(expMat_filter.hspcT)==colnames(expMat_filter.Thspc))

pMat_filter <- rbind(pMat_filter.Thspc, pMat_filter.hspcT)
expMat_filter <- rbind(expMat_filter.Thspc, expMat_filter.hspcT)

Thspc$interaction <- hspcT$count + Thspc$count
Thspc$count <- NULL

cell.orderT <- c("Naive","Memory","Effector")
cell.orderhspc <- c("HSCMPP","LMPP","MEP","MLP","EBM","Neu1","Neu2","MD1","MD2")
Thspc$source.cellGroup <- factor(Thspc$source.cellGroup,levels=cell.orderT,order=T)
Thspc$target.cellGroup <- factor(Thspc$target.cellGroup,levels=cell.orderhspc,order=T)

## remained the interaction within each disease type
type <- c("HD","CAA","SAA")
meta_filter <- c()
for(i in type){
     this.meta <- Thspc[Thspc$source.Disease_type == i & Thspc$target.Disease_type == i,]
     meta_filter <- rbind(meta_filter,this.meta)
	 }
meta_filter$source.Disease_type <- factor(meta_filter$source.Disease_type,levels=type,order=T)
## interaction distribution
ggplot(meta_filter,aes(x=interaction)) + geom_line(stat="density",colour="black")
ggsave(filename =paste0(dir,"/interaction_distribution.pdf"),height=3,width=6,onefile=F)

## remained non-SAA and Ctrl
nonSAACtrl.meta <- meta_filter[meta_filter$source.Disease_type != "SAA",]

## boxplot 
library(ggpubr)
p <- ggboxplot(nonSAACtrl.meta, x = "source.Disease_type", 
     y = "interaction",
	 fill="source.Disease_type",
     color = "black",
     outlier.shape=NA,
	 palette = c("DarkGray","OrangeRed"),
     facet.by = "source.phenotype",
	 title="Number of Interaction Pairs in non-SAA and Ctrl",
     xlab="",
	 order=c("HD","CAA"),
	 ylab="Number of Interaction Pairs",
	 short.panel.labs = TRUE)
     # Use only p.format as label. Remove method name.
     p + stat_compare_means(method="t.test",aes(label = paste0("p = ", ..p.format..)))
     ggsave(filename =paste0(dir,"/output/Interaction number boxplot non-SAA.pdf"),height=3,width=2.4,onefile=F)

pMat_filter_caa <- pMat_filter[,c(rownames(nonSAACtrl.meta[nonSAACtrl.meta$source.Disease_type == "CAA",]),"pair")]
pMat_filter_ctrl <- pMat_filter[,c(rownames(nonSAACtrl.meta[nonSAACtrl.meta$source.Disease_type == "HD",]),"pair")]

## retain significant pairs
## non-SAA
pMat_filter_caa_conver <- pMat_filter_caa
pMat_filter_caa_conver[1:54][pMat_filter_caa_conver[1:54] > 0.05] <- 1
pMat_filter_caa_conver <- pMat_filter_caa_conver[rowSums(pMat_filter_caa_conver[1:54]) < dim(pMat_filter_caa_conver[1:54])[2],]
pMat_filter_caa <- pMat_filter_caa[rownames(pMat_filter_caa_conver),]
expMat_filter_caa <- expMat_filter[rownames(pMat_filter_caa_conver),colnames(pMat_filter_caa)]
dim(pMat_filter_caa)   # 97 55
dim(expMat_filter_caa) # 97 55

## Ctrl
pMat_filter_ctrl_conver <- pMat_filter_ctrl
pMat_filter_ctrl_conver[1:54][pMat_filter_ctrl_conver[1:54] > 0.05] <- 1
pMat_filter_ctrl_conver <- pMat_filter_ctrl_conver[rowSums(pMat_filter_ctrl_conver[1:54]) < dim(pMat_filter_ctrl_conver[1:54])[2],]
pMat_filter_ctrl <- pMat_filter_ctrl[rownames(pMat_filter_ctrl_conver),]
expMat_filter_ctrl <- expMat_filter[rownames(pMat_filter_ctrl_conver),colnames(pMat_filter_ctrl)]
dim(pMat_filter_ctrl)   #   63 55
dim(expMat_filter_ctrl) #   63 55

## melt convert pvalue and expression matrix 
## dot plot
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
		    colnames(tab_split) <- paste0(i,c(".phenotype",".cellGroup",".Disease_type"))
		    Melt_pval_exp <- cbind(Melt_pval_exp,tab_split)
		    }
	     Melt_pval_exp$source.cellGroup <- factor(Melt_pval_exp$source.cellGroup,levels=cell.orderT,order=T)
	     Melt_pval_exp$target.cellGroup <- factor(Melt_pval_exp$target.cellGroup,levels=cell.orderhspc,order=T)
	     Melt_pval_exp$source.phenotype <- factor(Melt_pval_exp$source.phenotype,levels=c("CD4T","CD8T"),order=T)
	     return(Melt_pval_exp)
	   }

Melt_caa_sig <- melt.conver(pMat=pMat_filter_caa, expMat=expMat_filter_caa)
Melt_ctrl_sig <- melt.conver(pMat=pMat_filter_ctrl, expMat=expMat_filter_ctrl)
Melt_rbind_sig <- rbind(Melt_caa_sig,Melt_ctrl_sig)
Melt_rbind_sig_pval <- Melt_rbind_sig[Melt_rbind_sig$pvalue <= 0.05,]

## stat significant interactions
library(plyr)
count.data <- ddply(Melt_rbind_sig_pval, .(source.Disease_type,pair), summarize, counts= length(pair))
count.data <- count.data[order(count.data$counts),]

## Cell type specific pairs
library(plyr)
library(ggplot2)
Tcell.oder <- c("Naive","Memory","Effector")
caa_cell_uniq <- c()
for(i in Tcell.oder){
     for(j in c("CD4T","CD8T")){
	     k <- paste(j,i,sep=".")
         caa_sig <- Melt_caa_sig[Melt_caa_sig$source.cellGroup==i & Melt_caa_sig$source.phenotype==j ,]
    	 ctrl_sig <- Melt_ctrl_sig[Melt_ctrl_sig$source.cellGroup==i & Melt_ctrl_sig$source.phenotype==j,]
         caa_sig_pval <- caa_sig[caa_sig$pvalue <= 0.05,]
         ctrl_sig_pval <- ctrl_sig[ctrl_sig$pvalue <= 0.05,]
	     caa_sig_pval <- caa_sig_pval[!caa_sig_pval$pair %in% ctrl_sig_pval$pair,]
	     caa_cell_uniq <- rbind(caa_cell_uniq,caa_sig_pval)
	     count.data <- ddply(caa_sig_pval, .(pair), summarize, counts= length(pair))
	     count.data <- count.data[order(count.data$counts),]
         pairs <- unique(count.data$pair)
	     caa_sig <- caa_sig[caa_sig$pair %in% pairs,]
	     caa_sig$pair <- factor(caa_sig$pair,levels=pairs,order=T)
	     caa_sig$pvalue[caa_sig$pvalue > 0.05] <- 0.1
	     caa_sig$pvalue[caa_sig$pvalue <= 0.001] <- 0.001
		 caa_sig$expr[caa_sig$expr >= 2.0] <- 2.0
		 if(T){
	     ggplot(caa_sig,aes(x=factor(target.cellGroup,levels=cell.orderhspc),y=pair))+
       	     geom_point(aes(color=expr,size=-log10(pvalue)),alpha=0.8) +
			 scale_size_continuous(range = c(0.2,1.2)) + 
      	     scale_colour_gradient(low="#0000CD",high="#FF0000",limits=c(0.0,2.0)) +
         	 xlab("") + ylab("") +
             theme_bw()+theme(axis.text=element_text(size=7,color="black"),axis.text.x=element_text(angle=45, hjust=1,vjust=1),
             axis.text.y = element_text(size = 7,color="black"),axis.title.y = element_text(size = 7),
             panel.grid.minor = element_blank()) + 
             ggsave(filename =paste0(dir,"/output/dotplot_significant.pair.",j,".",i,".non-SAA.vs.Ctrl.pdf"),height=3.2,width=4,onefile=F)
		     }
		 if(T){
	     ggplot(caa_sig,aes(x=factor(target.cellGroup,levels=cell.orderhspc),y=pair))+
       	     geom_point(aes(color=expr,size=-log10(pvalue)),alpha=0.8) +
			 scale_size_continuous(range = c(0.2,1.2)) + 
      	     scale_colour_gradient(low="#0000CD",high="#FF0000",limits=c(0.0,2.0)) +
         	 xlab("") + ylab("") +
             theme_bw()+ theme(axis.text=element_text(size=7,color="black"),axis.text.x=element_text(angle=45, hjust=1,vjust=1),
     	     axis.text.y = element_blank(),axis.title.y = element_blank(),axis.title.x = element_text(size=7),
             panel.grid.minor = element_blank())+labs(title=paste0(j," ",i)) + theme(legend.position='none')
             ggsave(filename =paste0(dir,"/output/dotplot_significant.pair.",j,".",i,".non-SAA.vs.Ctrl.NoLegend.pdf"),height=3.2,width=0.7,onefile=F)
	         }
		 }
	 }
write.table(caa_cell_uniq,file=paste0(dir,"/output/caa_cell_uniq.xls"),row.names=F,sep="\t",quote=F)

## heatmap significant counts
cellGroup.source <- c("CD4T.Naive.CAA","CD4T.Memory.CAA","CD4T.Effector.CAA","CD8T.Naive.CAA","CD8T.Memory.CAA","CD8T.Effector.CAA")
cellGroup.target <- c("HSCMPP","LMPP","MEP","MLP","EBM","Neu1","Neu2","MD1","MD2")
Mat <- matrix(NA,nrow=length(cellGroup.source),ncol=length(cellGroup.target))
colnames(Mat) <- cellGroup.target
rownames(Mat) <- cellGroup.source
for(i in cellGroup.source){
     for(j in cellGroup.target){
	     tab <- caa_cell_uniq[caa_cell_uniq$source==i & caa_cell_uniq$target.cellGroup==j,]
         Mat[i,j] <- dim(tab)[1]
	     }
	 }
library(corrplot)
dat <- apply(Mat,2,as.numeric)
rownames(dat) <- rownames(Mat)
col <- colorRampPalette(c( 'white','orange','red3'))
pdf(file = paste0(dir,"/output/heatmap.non-SAA pair count.pdf"),7,5,onefile=F)
corrplot(dat, is.corr = FALSE,tl.cex = 1.5,addCoef.col = "grey",
     method="circle",
     addCoefasPercent = FALSE,
     col=col(200),
     cl.cex = 0.8,
     cl.ratio = 0.18,
     tl.col="black")
     dev.off()
	 
## VennDiagram
library(VennDiagram)
Melt_caa_sig_pval <- Melt_caa_sig[Melt_caa_sig$pvalue <= 0.05,]
Melt_ctrl_sig_pval <- Melt_ctrl_sig[Melt_ctrl_sig$pvalue <= 0.05,]

pdf(file=paste0(dir,"/output/non-SAAvsCtrl_venn diagram.pdf"),4,4,onefile=F)
p <- venn.diagram(list(CD4TCAA=Melt_caa_sig_pval[Melt_caa_sig_pval$source.phenotype=="CD4T",]$pair,
     CD4TCtrl=Melt_ctrl_sig_pval[Melt_ctrl_sig_pval$source.phenotype=="CD4T",]$pair,
	 CD8TCAA=Melt_caa_sig_pval[Melt_caa_sig_pval$source.phenotype=="CD8T",]$pair,
	 CD8TCtrl=Melt_ctrl_sig_pval[Melt_ctrl_sig_pval$source.phenotype=="CD8T",]$pair),
         filename= NULL,
         lwd=1,
         lty=3,
         main = "non-SAA vs Ctrl",
         main.cex =1.3,
         cex = 1.2,
         cat.cex = 1,
         col=c('steelblue1','darkorange1','#98FB98','#9370DB'),
         fill=c('steelblue1','darkorange1','#98FB98','#9370DB'),
         cat.col="black")
         grid.draw(p)
         dev.off()

## there were some cell clusters with less than 10 cells in SAA;
## when comparing SAA and non-SAA, we consider to remove those clusters.
## remove celltype with less than 10 cells (CD4T effector, Neu1 and EBM).
meta_filter$source.phenotype.cellGroup <- paste(meta_filter$source.phenotype,meta_filter$source.cellGroup,sep=".")
meta_filter_rm <- meta_filter[meta_filter$source.phenotype.cellGroup != "CD4T.Effector",]
meta_filter_rm <- meta_filter_rm[!meta_filter_rm$target.cellGroup %in% c("Neu1","EBM"),]
dim(meta_filter_rm) # 105  11

## boxplot three condition
comparisons <- list(c("SAA","HD"),c("CAA","HD"),c("SAA","CAA"))
library(ggpubr)
p <- ggboxplot(meta_filter_rm, x = "source.Disease_type", 
     y = "interaction",
	 fill="source.Disease_type",
     color = "black",
     outlier.shape=NA,
	 palette = "npg",
     facet.by = "source.phenotype",
	 order=c("HD","CAA","SAA"),
	 title="Number of Interaction Pairs in non-SAA and Ctrl",
     xlab="",
	 ylab="Number of Interaction Pairs",
	 short.panel.labs = TRUE)
     # Use only p.format as label. Remove method name.
     p + stat_compare_means(comparisons=comparisons,method="t.test",paired=TRUE, aes(label = paste0("p = ", ..p.format..)))
     ggsave(filename =paste0(dir,"/output/Interaction number boxplot in all type.three.condition.pdf"),height=4,width=4,onefile=F)

pMat_filter_saa <- pMat_filter[,c(rownames(meta_filter_rm[meta_filter_rm$source.Disease_type == "SAA",]),"pair")]
pMat_filter_caa <- pMat_filter[,c(rownames(meta_filter_rm[meta_filter_rm$source.Disease_type == "CAA",]),"pair")]
pMat_filter_ctrl <- pMat_filter[,c(rownames(meta_filter_rm[meta_filter_rm$source.Disease_type == "HD",]),"pair")]
pMat_filter_rm <- pMat_filter[,c(rownames(meta_filter_rm),"pair")]
expMat_filter_rm <- expMat_filter[,c(rownames(meta_filter_rm),"pair")]

## retain significant pairs
## SAA
pMat_filter_saa_conver <- pMat_filter_saa
pMat_filter_saa_conver[1:35][pMat_filter_saa_conver[1:35] > 0.05] <- 1
pMat_filter_saa_conver <- pMat_filter_saa_conver[rowSums(pMat_filter_saa_conver[1:35]) < dim(pMat_filter_saa_conver[1:35])[2],]
pMat_filter_saa <- pMat_filter_saa[rownames(pMat_filter_saa_conver),]
expMat_filter_saa <- expMat_filter[rownames(pMat_filter_saa_conver),colnames(pMat_filter_saa)]
dim(pMat_filter_saa)   # 101  36
dim(expMat_filter_saa) # 101  36

## non-SAA
pMat_filter_caa_conver <- pMat_filter_caa
pMat_filter_caa_conver[1:35][pMat_filter_caa_conver[1:35] > 0.05] <- 1
pMat_filter_caa_conver <- pMat_filter_caa_conver[rowSums(pMat_filter_caa_conver[1:35]) < dim(pMat_filter_caa_conver[1:35])[2],]
pMat_filter_caa <- pMat_filter_caa[rownames(pMat_filter_caa_conver),]
expMat_filter_caa <- expMat_filter[rownames(pMat_filter_caa_conver),colnames(pMat_filter_caa)]
dim(pMat_filter_caa)   # 86 36
dim(expMat_filter_caa) # 86 36

## Ctrl
pMat_filter_ctrl_conver <- pMat_filter_ctrl
pMat_filter_ctrl_conver[1:35][pMat_filter_ctrl_conver[1:35] > 0.05] <- 1
pMat_filter_ctrl_conver <- pMat_filter_ctrl_conver[rowSums(pMat_filter_ctrl_conver[1:35]) < dim(pMat_filter_ctrl_conver[1:35])[2],]
pMat_filter_ctrl <- pMat_filter_ctrl[rownames(pMat_filter_ctrl_conver),]
expMat_filter_ctrl <- expMat_filter[rownames(pMat_filter_ctrl_conver),colnames(pMat_filter_ctrl)]
dim(pMat_filter_ctrl)   # 57 36
dim(expMat_filter_ctrl) # 57 36

## melt convert pvalue and expression matrix 
Melt_saa_sig <- melt.conver(pMat=pMat_filter_saa, expMat=expMat_filter_saa)
Melt_caa_sig <- melt.conver(pMat=pMat_filter_caa, expMat=expMat_filter_caa)
Melt_ctrl_sig <- melt.conver(pMat=pMat_filter_ctrl, expMat=expMat_filter_ctrl)
	 
## remain significant pairs in each condition
Melt_saa_sig_pval <- Melt_saa_sig[Melt_saa_sig$pvalue <= 0.05,]
Melt_caa_sig_pval <- Melt_caa_sig[Melt_caa_sig$pvalue <= 0.05,]
Melt_ctrl_sig_pval <- Melt_ctrl_sig[Melt_ctrl_sig$pvalue <= 0.05,]

for(i in c("CD4T","CD8T")){
     pdf(file=paste0(dir,"/output/SAA.vs.nonSAA.vs.Ctrl_venn diagram_",i,".pdf"),4,4,onefile=F)
     p <- venn.diagram(list(SAA=Melt_saa_sig_pval[Melt_saa_sig_pval$source.phenotype==i,]$pair,CAA=Melt_caa_sig_pval[Melt_caa_sig_pval$source.phenotype==i,]$pair,
	      Ctrl=Melt_ctrl_sig_pval[Melt_ctrl_sig_pval$source.phenotype==i,]$pair),
          filename= NULL,
          lwd=1,
          lty=3,
          main = i,
          main.cex =1.3,
          cex = 1.2,
          cat.cex = 1,
          col=c('#DE5419','#4EB1C9','#B5B5B6'),
          fill=c('#DE5419','#4EB1C9','#B5B5B6'),
          cat.col="black")
          grid.draw(p)
          dev.off()
	 }

## remain significant pairs
Melt_sig_total <- melt.conver(pMat=pMat_filter_rm, expMat=expMat_filter_rm)
Melt_sig_total$phenotype.cellGroup <- paste(Melt_sig_total$source.phenotype,Melt_sig_total$source.cellGroup,sep=".")
phenotype.order <- c("CD4T.Naive","CD4T.Memory","CD8T.Naive","CD8T.Memory","CD8T.Effector")
Melt_sig_total$phenotype.cellGroup <- factor(Melt_sig_total$phenotype.cellGroup,levels=phenotype.order,order=T)
Melt_sig_total$phenotype.cellGroup.disease <- paste(Melt_sig_total$source.phenotype,Melt_sig_total$source.cellGroup,Melt_sig_total$source.Disease_type,sep=".")
Melt_sig_total$target.cellGroup <- factor(Melt_sig_total$target.cellGroup,levels=cell.orderhspc[-c(5:6)],order=T)

saa_caa_uniq <- c()
for(i in Tcell.oder){
     for(j in c("CD8T","CD4T")){
	     k <- paste0(j,".",i)
		 if(k != "CD4T.Effector"){
         saa_sig <- Melt_saa_sig_pval[Melt_saa_sig_pval$source.cellGroup==i & Melt_saa_sig_pval$source.phenotype==j ,]
		 caa_sig <- Melt_caa_sig_pval[Melt_caa_sig_pval$source.cellGroup==i & Melt_caa_sig_pval$source.phenotype==j ,]
    	 ctrl_sig <- Melt_ctrl_sig_pval[Melt_ctrl_sig_pval$source.cellGroup==i & Melt_ctrl_sig_pval$source.phenotype==j,]       
		 saa_caa_sig <- rbind(saa_sig,caa_sig)
	     saa_caa_sig_uniq <- saa_caa_sig[!saa_caa_sig$pair %in% ctrl_sig$pair,]
	     saa_caa_uniq <- rbind(saa_caa_uniq,saa_caa_sig_uniq)
	     count.data <- ddply(saa_caa_sig_uniq, .(pair), summarize, counts= length(pair))
	     count.data <- count.data[order(count.data$counts),]
         pairs <- unique(count.data$pair)
	     Melt_sig_pair <- Melt_sig_total[Melt_sig_total$pair %in% pairs & Melt_sig_total$source.Disease_type != "HD" & Melt_sig_total$phenotype.cellGroup == k,]
	     Melt_sig_pair$pair <- factor(Melt_sig_pair$pair,levels=pairs,order=T)
	     Melt_sig_pair$pvalue[Melt_sig_pair$pvalue > 0.05] <- 0.1
	     Melt_sig_pair$pvalue[Melt_sig_pair$pvalue <= 0.001] <- 0.001
		 Melt_sig_pair$source.Disease_type <- factor(Melt_sig_pair$source.Disease_type, levels=c("CAA","SAA"))
		 Melt_sig_pair$expr[Melt_sig_pair$expr >= 2.0] <- 2.0
		 if(T){
	     ggplot(Melt_sig_pair,aes(x=factor(target.cellGroup),y=pair))+
       	         geom_point(aes(color=expr,size=-log10(pvalue)),alpha=0.8) +
				 scale_size_continuous(range = c(0.2,1.2)) + 
      	         scale_colour_gradient(low="#0000CD",high="#FF0000",limits=c(0.0,2.0)) +
         	     xlab("") + ylab("") +
                 theme_bw()+theme(axis.text=element_text(size=7,color="black"),axis.text.x=element_text(angle=45, hjust=1,vjust=1),
     	         axis.text.y = element_text(size = 7,color="black"),axis.title.y = element_text(size = 7),
                 panel.grid.minor = element_blank())+labs(title=paste0(j," ",i)) + 
				 facet_grid(~source.Disease_type, scales = "free",space="free")
                 ggsave(filename =paste0(dir,"/output/dotplot_significant.pair.",j,".",i,".SAA.CAA.pdf"),height=5,width=6,onefile=F)
	             }
		 if(T){
	     ggplot(Melt_sig_pair,aes(x=factor(target.cellGroup),y=pair))+
       	     geom_point(aes(color=expr,size=-log10(pvalue)),alpha=0.8) +
			 scale_size_continuous(range = c(0.2,1.2)) + 
      	     scale_colour_gradient(low="#0000CD",high="#FF0000",limits=c(0.0,2.0)) +
         	 xlab("") + ylab("") +
             theme_bw()+ 
			 facet_grid(~source.Disease_type, scales = "free",space="free") +
			 theme(axis.text=element_text(size=7,color="black"),axis.text.x=element_text(angle=45, hjust=1,vjust=1),
     	     axis.text.y = element_blank(),axis.title.y = element_blank(),
             panel.grid.minor = element_blank())+labs(title=paste0(j," ",i)) + theme(legend.position='none')			 
             ggsave(filename =paste0(dir,"/output/dotplot_significant.pair.",j,".",i,".SAA.vs.CAA.NoLegend.pdf"),height=5,width=1.4,onefile=F)
	         }
			}
	     }
     }

## plot in one figure
count.data <- ddply(saa_caa_uniq, .(source,source.Disease_type,pair), summarize, counts= length(pair))
count.data <- count.data[order(count.data$counts),]
pairs <- unique(count.data$pair)
over_pair <- intersect(as.vector(count.data[count.data$source.Disease_type=="SAA",]$pair),as.vector(count.data[count.data$source.Disease_type=="CAA",]$pair))
over_pair <- pairs[pairs %in% over_pair]
saa_uniq_pair <- setdiff(as.vector(count.data[count.data$source.Disease_type=="SAA",]$pair),as.vector(Melt_caa_sig_pval$pair))
saa_uniq_pair <- pairs[pairs %in% saa_uniq_pair]
caa_uniq_pair <- setdiff(as.vector(count.data[count.data$source.Disease_type=="CAA",]$pair),as.vector(Melt_saa_sig_pval$pair))
caa_uniq_pair <- pairs[pairs %in% caa_uniq_pair]
pairs_order <- c(over_pair,saa_uniq_pair,caa_uniq_pair)
Melt_sig_pair <- Melt_sig_total[Melt_sig_total$pair %in% pairs_order & Melt_sig_total$source.Disease_type != "HD",]
Melt_sig_pair$source.Disease_type <- as.character(Melt_sig_pair$source.Disease_type)
Melt_sig_pair$source.Disease_type <- factor(Melt_sig_pair$source.Disease_type, levels=c("CAA","SAA"))
source.order <- c("CD4T.Naive.CAA","CD4T.Memory.CAA","CD8T.Naive.CAA","CD8T.Memory.CAA","CD8T.Effector.CAA","CD4T.Naive.SAA","CD4T.Memory.SAA","CD8T.Naive.SAA","CD8T.Memory.SAA","CD8T.Effector.SAA")
Melt_sig_pair$source <- factor(Melt_sig_pair$source,levels=source.order,order=T)
Melt_sig_pair$pvalue[Melt_sig_pair$pvalue > 0.05] <- 0.1
Melt_sig_pair$pvalue[Melt_sig_pair$pvalue <= 0.001] <- 0.001
Melt_sig_pair$tag <- ""
Melt_sig_pair$tag[Melt_sig_pair$pair %in% over_pair] <- "overlap"
Melt_sig_pair$tag[Melt_sig_pair$pair %in% saa_uniq_pair] <- "SAA"
Melt_sig_pair$tag[Melt_sig_pair$pair %in% caa_uniq_pair] <- "non-SAA"
Melt_sig_pair$tag <- factor(Melt_sig_pair$tag,levels=c("overlap","SAA","non-SAA"))
Melt_sig_pair$expr[Melt_sig_pair$expr >= 2.0] <- 2.0

ggplot(Melt_sig_pair,aes(x=factor(target.cellGroup,levels=cell.orderhspc[-c(5:6)]),y=pair))+
     geom_point(aes(color=expr,size=-log10(pvalue)),alpha=0.8) +
	 scale_size_continuous(range = c(0.2,1.2)) + 
     scale_colour_gradient(low="#0000CD",high="#FF0000",limits=c(0.0,2.0)) +
     xlab("") + ylab("") +
     theme_bw()+theme(axis.text=element_text(size=7,color="black"),axis.text.x=element_text(angle=45, hjust=1,vjust=1),
     axis.text.y = element_text(size = 7,color="black"),axis.title.y = element_text(size = 7),
     panel.grid.minor = element_blank()) + 
	 facet_grid(tag~source, scales = "free",space="free")
     ggsave(filename =paste0(dir,"/output/dotplot_significant.pair.total.pair.SAA.non-SAA.pdf"),height=8,width=8,onefile=F)