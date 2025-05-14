message("Converting FreeBayes.txt file into InformativeSNPs.txt and MpileupInput. Filtering for presence in 1000 Genomes, SNPs, and Heterozygosity")
Data <- read.table("FB.txt", sep = "\t", header = F, stringsAsFactors = F, quote = "", skip = 1)
Data$V17[is.na(Data$V17)] <- 0
Columns <- ncol(Data)
Samples <- Columns - 56
cat(Samples, "Sample(s) Found")

if (Samples > 1) {
  for (i in 57:Columns){
    q <- i-56
    #Data[,i] <- sapply(strsplit(Data[,i],"\\\""), "[", 2)
    GTSample <- paste("GT", q, sep = "")
    assign(GTSample, unlist(strsplit(Data[1,i],":")))
    GTlength <- as.numeric(length(get(GTSample)))
    GTMutpos <- GTlength-2
    GTTotal <- paste("GT_Total", q, sep = "_")
    GTMut <- paste("GT_Mut", q, sep = "_")
    assign(GTTotal, as.numeric(sapply(strsplit(Data[,i],":"), "[", GTlength)))
    assign(GTMut, as.numeric(sapply(strsplit(Data[,i],":"), "[", GTMutpos)))
    Dot1 <- paste("Dot_1", q, sep = "_")
    Dot2 <- paste("Dot_2", q, sep = "_")
    assign(Dot2, sapply(strsplit(sapply(strsplit(Data[,i],":"), "[", 1),"/"), "[", 2))
    assign(Dot1, sapply(strsplit(sapply(strsplit(Data[,i],":"), "[", 1),"/"), "[", 1))
  }
  DotPres <- ifelse((Dot_1_1 == "." | (Dot_2_1 == "." | (is.na(Dot_2_1) == "TRUE"))), "Yes", "No")
  for (j in 2:Samples) {
    GTMutTert <- paste("GT_Mut", j, sep = "_")
    GTTotalTert <- paste("GT_Total", j, sep = "_")
    GT_Mut_1 <- GT_Mut_1+get(GTMutTert)
    GT_Total_1 <- GT_Total_1+get(GTTotalTert)
    Dot1Tert <- paste("Dot_1", j, sep = "_")
    Dot2Tert <- paste("Dot_2", j, sep = "_")
    DotPres <- ifelse(DotPres == "Yes" | get(Dot1Tert) == "." | (get(Dot2Tert) == "." | (is.na(get(Dot2Tert)) == "TRUE")), "Yes", "No")
}
  MAF <- GT_Mut_1/GT_Total_1
  Data1000 <- Data[c(which(Data$V17 >= .01 & (Data$V4 == "C" | Data$V4 == "T" | Data$V4 == "G" | Data$V4 == "A") & (Data$V5 == "C" | Data$V5 == "T" | Data$V5 == "G" | Data$V5 == "A") & MAF >= .30 & MAF <= .70 & DotPres == "No")), ]
} else {
  GT <- unlist(strsplit(Data$V57[1],":"))
  GTlength <- length(GT)
  GTpercent <- GTlength-1
  Dot1 <- sapply(strsplit(sapply(strsplit(Data[,57],":"), "[", 1),"/"), "[", 1)
  Dot2 <- sapply(strsplit(sapply(strsplit(Data[,57],":"), "[", 1),"/"), "[", 2)
  DotPres <- ifelse((Dot1 == "." | Dot2 == "."), "Yes", "No")
  MAF <- as.numeric(sapply(strsplit(Data[,57],":"), "[", GTpercent))
  Data1000 <- Data[c(which(Data$V17 >= .01 & (Data$V4 == "C" | Data$V4 == "T" | Data$V4 == "G" | Data$V4 == "A") & (Data$V5 == "C" | Data$V5 == "T" | Data$V5 == "G" | Data$V5 == "A") & MAF >= 30 & MAF <= 70 & DotPres == "No")), ]
}

write.table(Data1000, "InformativeSNPs.txt", sep = "\t", row.names = F, col.names = T, quote = F, na = "")
Data1000$V1 <- sapply(strsplit(Data1000$V1, "r"), "[", 2)
Data1000$V58 <- rep("dummy", nrow(Data1000))
Mpileup <- rbind(1:8, 1:8, Data1000[, c(58, 1:3, 58, 58, 4:6)])
write.table(Mpileup, "MpileupInput.txt", sep = "\t", row.names = F, col.names = F, quote = F, na = "")
