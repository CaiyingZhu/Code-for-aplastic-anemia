## R version 3.5.0; Seurat version 2.3.4
rm(list=ls())
pwd <- getwd()
library(Seurat)
library(ggplot2)
load(paste0(pwd,"/int/02.Seurat.AACtrl.aligned.Rdata"))

## Cell type composition and distribution
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

hspc.Roe <- Roe(hspc.aligned.obj@meta.data, condition="Tissue", cellType="cellGroup", samples="SampleID", ctrl=c("Ctrl1","Ctrl4"))

roe.rescale <- hspc.Roe$samples
roe.rescale$O2E[roe.rescale$O2E > 6] <- 6

ggplot(data=hspc.Roe$group, aes(x=factor(cellType, levels=new.ident.order), y=meanO2E)) +
	geom_bar(stat="identity", alpha=0.8, color="white", fill="grey") +
	geom_point(data=roe.rescale, aes(x=cellType, y=O2E, color=samples, size=-log10(chisq.p)), alpha=0.5) +
	geom_errorbar(aes(ymin=meanO2E-sem, ymax=meanO2E+sem), width=.2,position=position_dodge(.9)) +
	geom_hline(yintercept = 1, color="grey", linetype="dashed") +
	ggtitle("hspc cell type distribution in AA vs. Ctrl") +
	guides(fill=guide_legend(title=""), color=guide_legend(title="")) +
	#guides(fill="none", color="none", size="none") +
	xlab("") + ylab("Ro/e") + theme_classic() +
	theme(axis.text.x=element_text(angle=45),axis.ticks.length = unit(.3, "cm"))
ggsave(paste0(pwd, "/output/hspc.Roe.newOrder.pdf"), width=6, height=7)

## separate SAA and non-SAA
for(i in c("SAA","non-SAA")){
     sub.obj <- SubsetData(object =hspc.aligned.obj, cells.use = hspc.aligned.obj@cell.names[hspc.aligned.obj@meta.data$Disease_types %in% c("Ctrl",i)], subset.raw =T)
     sub.Roe <- Roe(sub.obj@meta.data, condition="Tissue", cellType="cellGroup", samples="SampleID", ctrl=c("Ctrl1","Ctrl4"))
     roe.rescale <- sub.Roe$samples
     roe.rescale$O2E[roe.rescale$O2E > 6] <- 6
     ggplot(data=sub.Roe$group, aes(x=factor(cellType, levels=new.ident.order), y=meanO2E)) +
         geom_bar(stat="identity", alpha=0.8, color="white", fill="grey") +
     	 geom_point(data=roe.rescale, aes(x=cellType, y=O2E, color=samples, size=-log10(chisq.p)), alpha=0.5) +
     	 geom_errorbar(aes(ymin=meanO2E-sem, ymax=meanO2E+sem), width=.2,position=position_dodge(.9)) +
	     geom_hline(yintercept = 1, color="grey", linetype="dashed") +
	     ggtitle(paste0("hspc cell type distribution in ",i," vs. Ctrl")) +
	     guides(fill=guide_legend(title=""), color=guide_legend(title="")) +
	     xlab("") + ylab("Ro/e") + theme_classic() +
	     theme(axis.text.x=element_text(angle=45),axis.ticks.length = unit(.3, "cm"))
         ggsave(paste0(pwd, "/output/hspc.Roe.newOrder.",i,".pdf"), width=4, height=7)
     }
	 
## Fraction of clusters by patients
suppressPackageStartupMessages(library(dplyr))
hspc.counts <- sapply(unique(hspc.aligned.obj@meta.data$SampleID), function(x){
	 table(as.factor(hspc.aligned.obj@meta.data$cellGroup)[ hspc.aligned.obj@meta.data$SampleID == x])
     })
hspc.counts <- as.data.frame(t(prop.table(hspc.counts, 2)), stringsAsFactors=F)
hspc.counts.melt <- reshape2::melt(hspc.counts %>% tibble::rownames_to_column( var = "group" ))
names(hspc.counts.melt) <- c("Sample", "Cluster", "Percentage")
hspc.counts.melt$Cluster <- factor(hspc.counts.melt$Cluster, levels = new.ident.order, ordered=T)
hspc.counts.melt <- hspc.counts.melt[order(hspc.counts.melt$Cluster),]
oneCluster <- "HSC/MPP"
hspc.oneCluster <- hspc.counts.melt[hspc.counts.melt$Cluster == oneCluster & !hspc.counts.melt$Sample %in% c("Ctrl1","Ctrl4"), ]
hspc.sample.order <- c("Ctrl1","Ctrl4", as.character(hspc.oneCluster$Sample)[ order(hspc.oneCluster$Percentage,decreasing=T) ])

ggplot(hspc.counts.melt, aes(x=factor(Sample, levels=hspc.sample.order), y=Percentage)) +
	geom_bar(aes(fill=factor(Cluster, levels=new.ident.order)), stat="identity") +
	geom_hline(yintercept=cumsum(rev(hspc.counts.melt$Percentage[ hspc.counts.melt$Sample == "Ctrl1"])), linetype="dashed") +
	guides(fill=guide_legend("")) +
	xlab("") + theme(axis.text.x=element_text(angle=90))
ggsave(paste0(pwd, "/output/hspc.clusterFrequency_byPatient.pdf"))

## order samples by each cluster
## show the position of SAA and 
## perform wilcox.test to evaluate
## the cellular composition bias in
## SAA vs. non-SAA
sample.SAA <- c("P10", "P11", "P15")
hspc.counts.melt$type <- "non-SAA"
hspc.counts.melt$type[ hspc.counts.melt$Sample %in% sample.SAA] <- "SAA"
hspc.counts.melt$type[ hspc.counts.melt$Sample %in% c("Ctrl1","Ctrl4")] <- "Ctrl"

sample.df <- c()
pvalues <- c()
for( i in new.ident.order){
	 this.idx <- order(hspc.counts.melt$Percentage[ hspc.counts.melt$Cluster == i ])
	 sample.df <- rbind(sample.df, hspc.counts.melt$Sample[ this.idx ])
	 this.saa <- hspc.counts.melt$Percentage[ hspc.counts.melt$type == "SAA" & hspc.counts.melt$Cluster== i]
	 this.nonsaa <- hspc.counts.melt$Percentage[ hspc.counts.melt$type == "non-SAA" & hspc.counts.melt$Cluster== i]
	 this.Pvalue <- wilcox.test(this.saa, this.nonsaa)$p.value
	 pvalues <- c(pvalues, this.Pvalue)
     }
	 
sample.df.digits <- matrix(0, nrow=nrow(sample.df), ncol=ncol(sample.df))
row.names(sample.df.digits) <- paste0(new.ident.order, " (p=", format(pvalues,digits=2), ")")
sample.df.digits[ sample.df %in% sample.SAA] <- 1
sample.df.digits[ sample.df %in% c("Ctrl1","Ctrl4")] <- -1

pdf(paste0(pwd, "/output/sampleOrder_by_cellular_percentage.heatmap.pdf"), height=2.5)
pheatmap::pheatmap(
	 sample.df.digits,
	 color=gplots::colorpanel(n=3, low="blue", mid="grey", high="red"),
	 cluster_rows=F, cluster_cols=F,
	 show_rownames = T, show_colnames = F,
	 fontsize_row = 12, border="white",
	 legend_breaks = c(-1, 0, 1), legend_labels=c("Ctrl", "non-SAA", "SAA")
     )
dev.off()

# percentage by non-SAA and SAA
library(ggplot2)
target.groups <- c("HSC/MPP", "Neu1")
for(i in target.groups){
	 hsc.saa <- hspc.counts.melt$Percentage[ hspc.counts.melt$type == "SAA" & hspc.counts.melt$Cluster==i]
	 hsc.caa <- hspc.counts.melt$Percentage[ hspc.counts.melt$type == "non-SAA" & hspc.counts.melt$Cluster==i]
	 wilcox.test.Pvalue <- format(wilcox.test(hsc.saa, hsc.caa)$p.value, digits=2)
	 ggplot(subset(hspc.counts.melt, Cluster==i & type != "Ctrl"), aes(x=type, y=Percentage)) +
		 geom_boxplot(aes(fill=factor(type, levels=c("non-SAA", "SAA"))), outlier.shape = NA) + theme_classic() + 
		 guides(fill=guide_legend("none")) + 
		 ylab(paste0("Portion of ", i, " in CD34+ cells")) + xlab("") +
		 theme(axis.text=element_text(size=12), axis.title=element_text(size=15)) +
		 ggtitle(paste0(i, "\nwilcox.test, p value=", wilcox.test.Pvalue))
	 ggsave(paste0(pwd,"/output/", sub("/",".",i),".precentage_non-SAA_SAA.pdf"), height=4, width=3)
     }
