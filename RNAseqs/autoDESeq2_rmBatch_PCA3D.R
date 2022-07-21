## 2022-07-21 Kaining Hu
autoDESeq2 <- function(I,C,n){
library(DESeq2)
inputfilename <- I
Allcds <- read.csv(inputfilename,header = T, sep="\t",row.names = 1,stringsAsFactors = F)
congroupfile <- C
Congroup <- read.csv(congroupfile,sep="\t",header = T,stringsAsFactors = T)
uniqCongroup <- unique(as.character(Congroup$Group))
Congroups <- combn(uniqCongroup,2)
dds <- DESeqDataSetFromMatrix(Allcds,colData = Congroup,design =~ Group + Batch)
dds <- DESeq(dds)
write.table(counts(dds,normalized=T),"normalizedallCDS.txt", sep="\t", col.names = NA)
write.table(sizeFactors(dds),"sizefactors.txt", sep="\t", col.names = NA)
rld <- rlog(dds,blind=F)
write.table(assay(rld),"rlog.txt",sep="\t", col.names = NA)
write.table(cor(Allcds),"corAllcds.txt",sep="\t", col.names = NA)
write.table(cor(counts(dds)),"corcounts.txt",sep="\t", col.names = NA)
library("pheatmap")
library("RColorBrewer")
sampleDists <- dist( t( assay(rld) ) )
sampleDistMatrix <- as.matrix( sampleDists )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
compicname<-"CorHeatmap"
pdfname<-paste(compicname,".pdf")
pdf(pdfname)
pheatmap(sampleDistMatrix,
clustering_distance_rows=sampleDists,
clustering_distance_cols=sampleDists,
col=colors, fontsize = 4)
dev.off()
tiffname<-paste(compicname,".tiff")
tiff(tiffname, width = 2000,height = 2000)
pheatmap(sampleDistMatrix,
clustering_distance_rows=sampleDists,
clustering_distance_cols=sampleDists,
col=colors, fontsize = 14)
dev.off()
pngname<-paste(compicname,".png")
png(pngname, width = 2000,height = 2000)
pheatmap(sampleDistMatrix,
clustering_distance_rows=sampleDists,
clustering_distance_cols=sampleDists,
col=colors, fontsize = 14)
dev.off()


autoPCA <- function(i,n,d,G){
compicname<-paste("PCA",i)
txtname <- paste(compicname,".txt")
write.table(plotPCA(d,intgroup=G,ntop=n,returnData=T),file=txtname,sep="\t", col.names = NA)
pcap <- plotPCA(object = d,intgroup=G,ntop=n)
pdfname<-paste(compicname,".pdf")
pdf(pdfname)
print(pcap)
dev.off()
tiffname<-paste(compicname,".tiff")
tiff(tiffname)
print(pcap)
dev.off()
pngname<-paste(compicname,".png")
png(pngname)
print(pcap)
dev.off()
}

#library("btools")

plotPCA3D <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE){
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(PC1 = pca$x[, 1],
                  PC2 = pca$x[, 2],
                  PC3 = pca$x[, 3],
                  group = group,
                  intgroup.df,
                  name = colnames(object))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:3]
    return(d)
  }
  message("Generating plotly plot")
  p <- plotly::plot_ly(data = d,
                       x = ~PC1,
                       y = ~PC2,
                       z = ~PC3,
                       color = group,
                       mode = "markers",
                       type = "scatter3d")
  return(p)
}


autoPCA3D <- function(i,n,d,G){
  compicname<-paste("PCA3D",i)
  txtname <- paste(compicname,".txt")
  write.table(plotPCA3D(d,intgroup=G,ntop=n,returnData=T),file=txtname,sep="\t", col.names = NA)
  pcap <- plotPCA3D(object = d,intgroup=G,ntop=n)
  pdfname<-paste(compicname,".pdf")
  pdf(pdfname)
  print(pcap)
  dev.off()
  tiffname<-paste(compicname,".tiff")
  tiff(tiffname)
  print(pcap)
  dev.off()
  pngname<-paste(compicname,".png")
  png(pngname)
  print(pcap)
  dev.off()
}


autoPCA("rawGroup",n,rld,"Group")
autoPCA("rawGroupBatch",n,rld,c("Group","Batch"))
autoPCA3D("rawGroup",n,rld,"Group")
autoPCA3D("rawGroupBatch",n,rld,c("Group","Batch"))


vsd <- vst(dds, blind=FALSE)
write.table(assay(vsd),"vsd.txt",sep="\t", col.names = NA)
mat <- assay(vsd)
mm <- model.matrix(~Group, colData(vsd))
mat <- limma::removeBatchEffect(mat, batch=vsd$Batch, design=mm)
assay(vsd) <- mat
write.table(assay(vsd),"vsd_rm_batch.txt",sep="\t",col.names = NA)

autoPCA("rmBatch_Group",n,vsd,"Group")
autoPCA("rmBatch_GroupBatch",n,vsd,c("Group","Batch"))
autoPCA3D("rmBatch_Group",n,vsd,"Group")
autoPCA3D("rmBatch_GroupBatch",n,vsd,c("Group","Batch"))



sampleDists <- dist( t( assay(vsd) ) )
sampleDistMatrix <- as.matrix( sampleDists )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
compicname<-"CorHeatmap_rm_Batch"
pdfname<-paste(compicname,".pdf")
pdf(pdfname)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         
         col=colors, fontsize = 4)
dev.off()
tiffname<-paste(compicname,".tiff")
tiff(tiffname, width = 2000,height = 2000)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
        
         col=colors, fontsize = 24)
dev.off()
pngname<-paste(compicname,".png")
png(pngname, width = 2000,height = 2000)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
        
         col=colors, fontsize = 24)
dev.off()







autores <- function(G,A,B){
comtmp <- paste("res",A,"vs",B,".txt",sep="_")
compicname <- paste("res",A,"vs",B,"MA",sep="_")
comname <- comtmp
res <-results(dds,contrast=c(G,A,B),independentFiltering =F)
write.table(res,comname,sep="\t",col.names=NA)
p <- plotMA(res,alpha=0.05,ylim=c(-10,10))
pdfname<-paste(compicname,".pdf")
pdf(pdfname)
plotMA(res,alpha=0.05,ylim=c(-10,10))
dev.off()
tiffname<-paste(compicname,".tiff")
tiff(tiffname)
plotMA(res,alpha=0.05,ylim=c(-10,10))
dev.off()
pngname<-paste(compicname,".png")
png(pngname)
plotMA(res,alpha=0.05,ylim=c(-10,10))
dev.off()
return(comname)
}

for (i in 1:(length(Congroups)/2)) {
    if (file.exists( paste("res",Congroups[1,i],"vs",Congroups[2,i],".txt",sep="_"))){
    next
}
autores("Group",Congroups[1,i],Congroups[2,i])
}
}