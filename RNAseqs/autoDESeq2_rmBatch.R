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
    write.table(counts(dds,normalized=T),"normalizedallCDS.txt",sep="\t")
    write.table(sizeFactors(dds),"sizefactors.txt")
    rld <- rlog(dds,blind=F)
    write.table(assay(rld),"rlog.txt",sep="\t")
    write.table(cor(Allcds),"corAllcds.txt",sep="\t")
    write.table(cor(counts(dds)),"corcounts.txt",sep="\t")
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
    autoPCA <- function(i,n,m){
        compicname<-paste(i,"PCA")
        txtname <- paste(compicname,".txt")
        write.table(plotPCA(m,intgroup="Group",ntop=n,returnData=T),file=txtname,sep="\t")
        pcap <- plotPCA(object = m,intgroup="Group",ntop=n)
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
    
    
    
    autoPCA("raw",n,rld)
    
    vsd <- vst(dds, blind=FALSE)
    mat <- assay(vsd)
    mm <- model.matrix(~Group, colData(vsd))
    mat <- limma::removeBatchEffect(mat, batch=vsd$Batch, design=mm)
    assay(vsd) <- mat
    autoPCA("rmBatch",n,vsd)
    write.table(assay(vsd),"vsd.txt",sep="\t",col.names = NA)
    
    
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