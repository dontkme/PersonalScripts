### AutoIR.R. Kaining Hu edited on 2021-01-27

library(DESeq2)
# library(dplyr)
library(tidyverse)
source("../DESeq2Constructor.R")  # Load IRFinder-related function



results = read.table("filePaths.txt")
paths = as.vector(results$V1)                                            # File names must be saved in a vector
experiment = read.table("experiment.txt",header=T)                       
experiment$Condition=factor(experiment$Condition,levels=c("WT","KO"))    # Set WT as the baseline in the analysis
rownames(experiment)=NULL                                                # Force removing rownames

# WARNING: make sure the rownames of `experiment` is set to NULL. 
# WARNING: users MUST check if the order of files in the `path` matches the order of samples in `experiment` before continue  

metaList=DESeqDataSetFromIRFinder(filePaths=paths, designMatrix=experiment, designFormula=~1)
# The above line generate a meta list contains four slots
# First slot is a DESeq2 Object that can be directly pass to DESeq2 analysis.  
# Second slot is a matrix for trimmed mean of intron depth
# Third slot  is a matrix for correct splicing depth flanking introns
# Fourth slot is a matrix for maximum splicing reads at either ends of introns
# We build a “null” regression model on the interception only. 
# A “real” model can be assigned either here directly, or in the downstream step. See below

dds = metaList$DESeq2Object                       # Extract DESeq2 Object with normalization factors ready
#colData(dds)                                      # Check design of matrix


# Please note that sample size has been doubled and one additional column "IRFinder" has been added.
# This is because IRFinder considers each sample has two sets of counts: one for reads inside intronic region and one for reads at splice site, indicating by "IR" and "Splice" respectively.
# "IRFinder" is considered as an additional variable in the GLM model.
# Please also be aware that size factors have been set to 1 for all samples. Re-estimation of size factors is NOT recommended and is going to bias the result.
# More details at the end of the instruction.

design(dds) = ~Condition + Condition:IRFinder     # Build a formula of GLM. Read below for more details. 
dds = DESeq(dds)                                  # Estimate parameters and fit to model
#resultsNames(dds)                                 # Check the actual variable name assigned by DESeq2
#write.table(counts(dds,normalized=T),"AllCountsNormalized.txt",sep="\t")
write.table(counts(dds,normalized=F),"AllCounts.txt",sep="\t",col.names = NA)


res.WT = results(dds, name = "ConditionWT.IRFinderIR")
#write.table(res.WT,"res.WT.txt",sep="\t",col.names = NA)
# This tests if the number of IR reads are significantly different from normal spliced reads, in the WT samples.
# We might only be interested in the "log2FoldChange" column, instead of the significance.
# This is because "log2FoldChange" represents log2(number of intronic reads/number of normal spliced reads).
# So we the value of (intronic reads/normal spliced reads) by

WT.IR_vs_Splice=2^res.WT$log2FoldChange

# As IR ratio is calculated as (intronic reads/(intronic reads+normal spliced reads))
# We can easily convert the above value to IR ratio by

IRratio.WT = WT.IR_vs_Splice/(1+WT.IR_vs_Splice)

# Similarly, we can get IR ratio in the KO samples
res.KO = results(dds, name = "ConditionKO.IRFinderIR")
write.table(res.KO,"res.KO.txt",sep="\t",col.names = NA)

KO.IR_vs_Splice=2^res.KO$log2FoldChange
IRratio.KO = KO.IR_vs_Splice/(1+KO.IR_vs_Splice)

# Finally we can test the difference of (intronic reads/normal spliced reads) ratio between WT and KO
res.diff = results(dds, contrast=list("ConditionKO.IRFinderIR","ConditionWT.IRFinderIR")) 
write.table(res.diff,"res.diff.txt",sep="\t",col.names = NA)


# We can plot the changes of IR ratio with p values
# In this example we defined significant IR changes as
# 1) IR changes no less than 10% (both direction) and 
# 2) with adjusted p values less than 0.05

IR.change = IRratio.KO - IRratio.WT
png("IR.change.png", width = 3000,height = 2000)
plot(IR.change,col=ifelse(res.diff$padj < 0.05 & abs(IR.change)>=0.1, "red", "black"),cex.axis = 2)
dev.off()

IRratioall <- cbind(Intron=row.names(res.diff),IRratio.WT,IRratio.KO,IR.change,baseMean=res.diff$baseMean,log2FoldChange=res.diff$log2FoldChange,lfcSE=res.diff$lfcSE,stat=res.diff$stat,pvalue=res.diff$pvalue,padj=res.diff$padj)
IRratioall <- as.data.frame(IRratioall)
IRratioall$padj <- as.numeric(as.character(IRratioall$padj))
IRratioall$log2FoldChange <- as.numeric(as.character(IRratioall$log2FoldChange))

IRratioall <- IRratioall %>%mutate(PadjSig=ifelse(padj<0.05,"TRUE","FALSE"),UPSig=ifelse(padj<0.05 & log2FoldChange >0,"TRUE","FALSE"),DownSig=ifelse(padj<0.05 & log2FoldChange <0,"TRUE","FALSE"))
write.table(IRratioall[,11:13]  %>% map(table),"IRSigTable.txt",sep="\t",col.names = NA)
write.table(IRratioall,"IRratioAll.txt", col.names = NA,sep="\t")
# IRPadjLT0.5 <- IRratioall %>% filter(padj<0.05)
# IRPadjLT0.5UP <- IRPadjLT0.5 %>% filter(log2FoldChange>0)
# IRPadjLT0.5Down <- IRPadjLT0.5 %>% filter(log2FoldChange<0)
# write.table(IRPadjLT0.5,"IRratioAll_padjLT0.05.txt",sep="\t",col.names = NA)
# write.table(IRPadjLT0.5UP,"IRratioAll_padjLT0.05_UP.txt",sep="\t",col.names = NA)
# write.table(IRPadjLT0.5Down,"IRratioAll_padjLT0.05_Down.txt",sep="\t",col.names = NA)
# PLEASE NOTE    
# You might see dots at the same horizontal level (same IR change) either marked as significant (red) or not (black)
# This is because the Wald test applied above is testing on fold changes instead of absolute changes
# For example, both changes from 0.01 to 0.11 and from 0.8 to 0.9 have absolute changes of 0.1.
# But the fold changes for them are 11 folds and 1.1 folds respectively, and lead to different p values