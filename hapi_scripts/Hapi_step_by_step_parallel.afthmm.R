options(stringsAsFactors=F)

args=commandArgs(T)

input1 <- args[1]
input2 <- args[2]
output <- args[3]

gmt <- read.table(input1,header=T)

library(HMM)
library(Hapi)
library(parallel)

### covert A/T/C/G to 0/1
hetDa <- gmt[,1:4]
ref <- hetDa$ref
alt <- hetDa$alt

gmtDa <- gmt[,-(1:4)]
gmtDa <- base2num(gmt = gmtDa, ref = ref, alt = alt)
head(gmtDa)

### define HMM probabilities
#hmm = initHMM(States=c("S","D"), Symbols=c("s","d"), 
#    transProbs=matrix(c(0.99999,0.00001,0.00001,0.99999),2),
#    emissionProbs=matrix(c(0.99,0.01,0.01,0.99),2), 
#    startProbs = c(0.5,0.5))
#hmm


### filter out genotyping errors
#gmtDa <- hapiFilterError(gmt = gmtDa, hmm = hmm, cl = cl)

gmtDa <- read.table(input2,header=T)

#filter_out <- paste(output,"HmmFilter",sep=".")
#write.table(gmtDa,filter_out,row.names=T,col.names=T,sep="\t",quote=F)

#stopCluster(cl)
## load Hapi again ##
#library(Hapi)

### select a subset of high-quality markers
gmtFrame <- hapiFrameSelection(gmt = gmtDa, n = 10) ###

### imputation
imputedFrame <- hapiImupte(gmt = gmtFrame, nSPT = 2, allowNA = 20)

head(imputedFrame)

### majority voting
draftHap <- hapiPhase(gmt = imputedFrame) ###
head(draftHap)

### check positions with cv-links
draftHap[draftHap$cvlink>=1,]

### identification of clusters of cv-links
cvCluster <- hapiCVCluster(draftHap = draftHap, cvlink = 2)
cvCluster


### determine hetSNPs in small regions involving multiple cv-links
filter <- c()
for (i in 1:nrow(cvCluster)) {
    filter <- c(filter, which (rownames(draftHap) >= cvCluster$left[i] & 
    rownames(draftHap) <= cvCluster$right[i]))
}

length(filter)

### filter out hetSNPs in complex regions and infer new draft haplotypes
if (length(filter) > 0) {
    imputedFrame <- imputedFrame[-filter, ]
    draftHap <- hapiPhase(imputedFrame)
} 

finalDraft <- hapiBlockMPR(draftHap = draftHap, gmtFrame = gmtFrame, cvlink = 1)

head(finalDraft)

consensusHap <- hapiAssemble(draftHap = finalDraft, gmt = gmtDa)

head(consensusHap)

outRdata <- paste(output,"RData",sep=".")
save.image(file = outRdata)

consensusHap <- hapiAssembleEnd(gmt = gmtDa, draftHap = finalDraft, 
                                consensusHap = consensusHap, k = 300)


### Haplotype 1
hap1 <- sum(consensusHap$hap1==0)
### Haplotype 2
hap2 <- sum(consensusHap$hap1==1)
### Number of unphased hetSNPs
hap7 <- sum(consensusHap$hap1==7)

### Accuracy
max(hap1, hap2)/sum(hap1, hap2)



### find hetSNP overlaps
snp <- which(rownames(hetDa) %in% rownames(consensusHap))

ref <- hetDa$ref[snp]
alt <- hetDa$alt[snp]

### convert back to A/T/C/G
consensusHap <- num2base(hap = consensusHap, ref = ref, alt = alt)
head(consensusHap)


### output all the information
hapOutput <- data.frame(gmt[snp,], consensusHap)
head(hapOutput)

write.table(hapOutput,output,row.names=T,col.names=T,sep="\t",quote=F)





