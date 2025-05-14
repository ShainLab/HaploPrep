message("Reading MpileupOutput.txt and Filtering InformativeSNPs.txt and converting to VCF + MpileupInformativeSNP")
Data <- read.csv("InformativeSNPs.txt", sep = "\t", header = T, stringsAsFactors = F, quote = "")
Mpileup <- read.csv("MpileupOutput.txt", sep = "\t", header = F, stringsAsFactors = F, quote = "")
Ref <- Mpileup[-1,9]
Alt <- Mpileup[-1,10]
MAF <- Alt/(Ref + Alt)
Total <- Ref+Alt
Data1 <- Data[c(which(MAF >= .4 & MAF <= .6 & Total >= 10)), ]
library("vcfR")
library("memisc")
library("R.utils")
FB_VCF <- read.vcfR("FB.vcf")
fixnames <- colnames(FB_VCF@fix)
gtnames <- colnames(FB_VCF@gt)
charnumb <- as.character(Data1[, 49])
charqual <- as.character(Data1[, 53])
colnum <- as.numeric(ncol(Data1))
message("Creating InformativeSNPs.vcf")
FB_VCF@fix <- as.matrix(cbind(Data1[,48],charnumb, Data1[,50:52], charqual, Data1[,54:55]))
FB_VCF@gt <- as.matrix(Data1[,c(56:colnum)])
colnames(FB_VCF@fix) <- fixnames
colnames(FB_VCF@gt) <- gtnames
#FB_VCF <- as.matrix(FB_VCF)
write.vcf(FB_VCF, file = "InformativeSNPs.vcf.gz")
gunzip("InformativeSNPs.vcf.gz", remove = F)
file.remove("InformativeSNPs.vcf.gz")
Data1$V1 <- sapply(strsplit(Data1$V1, "r"), "[", 2)
Data1$dummy <- rep("dummy", nrow(Data1))
coldum <- colnum+1
Mpileup <- rbind(1:8, 1:8, Data1[, c(coldum, 1, 2, 3, coldum, coldum, 4, 5, coldum)])
write.table(Mpileup, "MpileupInformativeSNPInput.txt", sep = "\t", row.names = F, col.names = F, quote = F, na = "")
