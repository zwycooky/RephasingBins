options(stringsAsFactors=F)

args=commandArgs(T)

input <- args[1]
output <- args[2]
pwddir <- args[3]
#input <- "ptg000025l.vcf"
#output <- "res"

gmt <- read.table(input,header=T, colClasses = c ("character"))

library(HMM)
library(Hapi)
library(parallel)
cl <- makeCluster(5)

hapiFilterError <- function(gmt, hmm=NULL, cl) {
    gmt <- data.frame(gmt)
    total <- nrow(gmt)
    if (is.null(hmm)) {
        hmm <- initHMM(States=c('F','M'), Symbols=c('f','m'), 
            transProbs=matrix(c(0.99999,0.00001,0.00001,0.99999),2),
            emissionProbs=matrix(c(0.99,0.01,0.01,0.99),2), 
            startProbs = c(0.5,0.5))
    }
    
    ####
    filterErrorFun <- function(gmt1, gmt2, hmm=NULL) {
        if (is.null(hmm)) {
            hmm <- initHMM(States=c('F','M'), Symbols=c('f','m'), 
                transProbs=matrix(c(0.99999,0.00001,0.00001,0.99999),2),
                emissionProbs=matrix(c(0.99,0.01,0.01,0.99),2), 
                startProbs = c(0.5,0.5))
        }
        idComp <- gmt1 == gmt2
        idKnownPos <- which(!is.na(idComp))
    
        if (length(idKnownPos)<=1) {
            genoError <- fastCorrectIdentityFun(c(1,1), position=NULL, hmm=hmm)
            return (as.numeric(idKnownPos[genoError]))
        }
    
        idComp <- as.numeric(idComp[!is.na(idComp)])
        genoError <- fastCorrectIdentityFun(idComp, position=NULL, hmm=hmm)
    
        #if (length(genoError)==0) {
        #  return (NULL)
        #}
        return (as.numeric(idKnownPos[genoError]))
    }


    ###############


    #####
    fastCorrectIdentityFun <- function(genoIdentity, position=NULL, hmm=NULL) {
        if (is.null(hmm)) {
            hmm <- initHMM(States=c('F','M'), Symbols=c('f','m'), 
                transProbs=matrix(c(0.99999,0.00001,0.00001,0.99999),2),
                emissionProbs=matrix(c(0.99,0.01,0.01,0.99),2), 
                startProbs = c(0.5,0.5))
        }
    
        genoSymbol <- ifelse(genoIdentity==0,hmm$Symbols[1],hmm$Symbols[2])
    
        library(HMM)

        correctGeno <- viterbi(hmm, genoSymbol)
        correctGeno <- ifelse(correctGeno==hmm$States[1], 0, 1)
    
        genoError <- correctGeno-genoIdentity
        genoError <- which(genoError != 0)
    
        return (genoError)
    }
    
    nSNP <- 10000000
	round_filter <- 1
    while (nrow(gmt) < nSNP | round_filter <= 2) {
        nSNP <- nrow(gmt)
		if (round_filter == 1) {
			for (i in 1:ncol(gmt)) {
				genoError <- parLapply(cl,gmt[,-i], function(x) 
					filterErrorFun(gmt[,i],x,hmm=hmm))
				#filter <- sort(unique(unlist(genoError)))
			
				filter2NA <- as.numeric(names(table(unlist(genoError)))[table(unlist(genoError)) > 2])
			
				#message (paste(length(filter), 'hetSNPs 
				#with potential genotyping errors are filtered out !',
				#sep=' '))
            
				if (length(filter2NA) > 0) {
					gmt <- gmt
					gmt[filter2NA,i] <- NA
				}
			}
			round_filter <- 2
		}else{
			round_filter <- 3
			for (i in 1:ncol(gmt)) {
				genoError <- parLapply(cl,gmt[,-i], function(x) 
					filterErrorFun(gmt[,i],x,hmm=hmm))
				filter <- sort(unique(unlist(genoError)))
			
				#message (paste(length(filter), 'hetSNPs 
				#with potential genotyping errors are filtered out !',
				#sep=' '))
            
				if (length(filter) > 0) {
					gmt <- gmt[-filter,]
				}
			}
		}
    }
    
    final <- nrow(gmt)
    
    message (paste(total-final, 
        'hetSNPs with potential genotyping errors are filtered out !',
        sep=' '))
    
    return (gmt)
}


### covert A/T/C/G to 0/1
hetDa <- gmt[,1:4]
ref <- hetDa$ref
alt <- hetDa$alt

gmtDa <- gmt[,-(1:4)]
gmtDa <- base2num(gmt = gmtDa, ref = ref, alt = alt)
head(gmtDa)

### define HMM probabilities
hmm = initHMM(States=c("S","D"), Symbols=c("s","d"), 
    transProbs=matrix(c(0.99999,0.00001,0.00001,0.99999),2),
    emissionProbs=matrix(c(0.99,0.01,0.01,0.99),2), 
    startProbs = c(0.5,0.5))
hmm


### filter out genotyping errors
gmtDa <- hapiFilterError(gmt = gmtDa, hmm = hmm, cl = cl)

final_filter <- NULL
for(i in 1:nrow(gmtDa)) {
	if (sum(gmtDa[i,] == 1,na.rm=T) > 4 & sum(gmtDa[i,] == 0,na.rm=T) > 4) {
		final_filter <- c(final_filter,i)
	}
}
gmtDa <- gmtDa[final_filter,]

filter_out <- paste(output,"HmmFilter",sep=".")
write.table(gmtDa,filter_out,row.names=T,col.names=T,sep="\t",quote=F)
stopCluster(cl)

## some errors were turn to NA ##
tmp_snp <- which(rownames(gmt) %in% rownames(gmtDa))
tmp_ref <- hetDa$ref[tmp_snp]
tmp_alt <- hetDa$alt[tmp_snp]

# imputation #
## load Hapi again ##
library(Hapi)

### select a subset of high-quality markers
gmtFrame <- hapiFrameSelection(gmt = gmtDa, n = 10) ###

### imputation
imputedFrame <- hapiImupte(gmt = gmtFrame, nSPT = 2, allowNA = 97)
tmp_gmtDa <- imputedFrame

head(imputedFrame)

for(i in 1:ncol(tmp_gmtDa)) {
	tmp_geno <- as.character(tmp_gmtDa[,i])
	if (sum(tmp_geno == "0",na.rm=T) > 0) {
		tmp_geno[which(tmp_geno == "0")] <- tmp_ref[which(tmp_geno == "0")]
	}
	if (sum(tmp_geno == "1",na.rm=T) > 0) {
		tmp_geno[which(tmp_geno == "1")] <- tmp_alt[which(tmp_geno == "1")]
	}
	tmp_gmtDa[,i] <- tmp_geno
}

gmt[tmp_snp,-(1:4)] <- tmp_gmtDa


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

hapiBlockMPR2_script <- paste(pwddir,"/hapi_scripts/hapiBlockMPR2.0.R",sep="")
source(hapiBlockMPR2_script)

### select a subset of high-quality markers
gmtFrame <- hapiFrameSelection(gmt = gmtDa, n = 10) ###

finalDraft <- hapiBlockMPR2(draftHap = draftHap, gmtFrame = gmtFrame, cvlink = 2, smallBlock = 5)

head(finalDraft)

## filter gmtDa for NA ##
remove_NA_col <- NULL
for (i in 1:ncol(gmtDa)) {
	if (sum(is.na(gmtDa[,i])) == nrow(gmtDa)) {
		remove_NA_col <- c(remove_NA_col, i)
	}
}
if (length(remove_NA_col > 0)) {
	gmt_filter <- gmtDa[,-remove_NA_col]
}else{
	gmt_filter <- gmtDa
}

hapiAssemble2_script <- paste(pwddir,"/hapi_scripts/hapiAssemble2.0.R",sep="")
source(hapiAssemble2_script)
consensusHap <- hapiAssemble(draftHap = finalDraft, gmt = gmt_filter)

outRdata <- paste(output,"RData",sep=".")
save.image(file = outRdata)

head(consensusHap)

#consensusHap <- hapiAssembleEnd(gmt = gmt_filter, draftHap = finalDraft, 
#                                consensusHap = consensusHap, k = 300)

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





