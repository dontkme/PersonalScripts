ENS2GS <- function(I,L,O){
  library(dplyr)
  Inputfile <- read.table(I, sep="\t", header = T)
  Ens2GSTable <- read.table(L, sep="\t", header = T)
  outputfile <- left_join(Inputfile, Ens2GSTable, by = c("X"="Geneid"))
  outfilename<- paste(O,"addGS.txt",sep = ".")
  write.table(outputfile, file = outfilename, sep="\t",row.names = F)
  
}