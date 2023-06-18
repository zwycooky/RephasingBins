options(stringsAsFactors=F)

binmap <- read.table("Fudingdabai_bins_re_anchored_acood_with_direction.txt",header=F)
co_file <- read.table("co_final.filter.filter.txt",header=F)
chr_len <- read.table("chromosome_length_re_anchored.txt",header=T)
phasing <- read.table("merged_binMarkers.re-phasing.co8.missing30.txt",header=T)
F1_gmap <- read.table("Fuding_gmapF1_pos.filter.txt",header=F)
fuding_genes <- read.table("fuding_gene_acoord.txt",header=F)

rownames(binmap) <- binmap[,2]
co_reg <- cbind(binmap[co_file[,3],3:4],binmap[co_file[,4],3:4])

co_file$co_res <- co_reg[,3] - co_reg[,2]
co_file$co_pos <- (co_reg[,2] + co_reg[,3])/2

F1_gmap$pos <- (F1_gmap[,3] + F1_gmap[,4])/2

sliding_window_cMMb <- function (chr_len, co_pos, mywin = 3000000, mystep = 1000000) {
	lft <- seq(1,chr_len,mystep)
	rht <- lft + mywin
	swmat <- cbind(lft,rht)
	res <- matrix(0,nrow=nrow(swmat),ncol=2)
	for (i in 1:nrow(swmat)) {
		conum <- sum(co_pos > swmat[i,1] & co_pos <= swmat[i,2])
		cocM <- round( ((conum / 109) * 100)/3 ,2)
		if (conum != 0) {
			pos <- sum(swmat[i,])/2
			res[i,] <- c(pos,cocM)
		}else{
			pos <- sum(swmat[i,])/2
			res[i,] <- c(pos,0)
		}
	}
	#res <- res[!is.na(res[,2]),] 
	return(res)
}

mat <- F1_gmap[F1_gmap[,1] == "chr1",6:7]
cM_Mb_GMAP <- function (mat) {
	mat1 <- mat[-1,]
	mat2 <- mat[-nrow(mat),]
	cM <- abs(mat1[,1] - mat2[,1])
	reg <- abs(mat1[,2] - mat2[,2]) / 1000000
	pos <- (mat1[,2] + mat2[,2]) / 2
	cM_reg <- round(cM / reg,2)
	res <- cbind(pos,cM_reg)
	return(res)
}

co_count <- NULL
gmapF1_cM <- NULL
gene_count <- NULL
chr_vec <- paste("chr",1:15,sep="")
for (i in 1:length(chr_vec)) {
	tmp <- co_file[co_file[,1] == chr_vec[i],]
	tmp <- tmp[order(tmp$co_pos),]
	tmp_len <- chr_len[chr_vec[i],2]
	res <- sliding_window_cMMb(tmp_len, tmp$co_pos)
	res <- cbind(rep(i,nrow(res)),res)
	co_count <- rbind(co_count,res)
	
	tmp_gene <- fuding_genes[fuding_genes[,2] == chr_vec[i],]
	gene_pos <- (tmp_gene[,3] + tmp_gene[,4]) / 2
	gene_res <- sliding_window_cMMb(tmp_len, gene_pos)
	gene_res <- cbind(rep(i,nrow(res)),gene_res)
	gene_count <- rbind(gene_count,gene_res)
	
	tmp_mat <- F1_gmap[F1_gmap[,1] == chr_vec[i],6:7]
	res_F1 <- cM_Mb_GMAP(tmp_mat)
	res_F1 <- cbind(rep(i,nrow(res_F1)), res_F1)
	gmapF1_cM <- rbind(gmapF1_cM,res_F1)
}

## start plot ##
## height: 100, weight: 185
## bottom margin 5, chr: 5-95
## left margin: 5, co_count: 10, interval between 2 chr: 2 (30)
## format physical chr length ##
max_chrlen <- max(chr_len[,2])
max_co_count1 <- max( co_count[,3] )
max_co_count2 <- max( gmapF1_cM[,3] )
max_gene_count <- max( gene_count[,3] )
rel_len <- 90

# format co_count #
conum_vec <- NULL
format_co_count <- matrix(0,nrow=nrow(co_count),ncol=3)
for (i in 1:nrow(co_count)) {
	chr <- co_count[i,1]
	pos <- (co_count[i,2] / max_chrlen) * rel_len
	conum <- (co_count[i,3] / max_co_count1) * 7.5
	conum_vec <- c(conum_vec,conum)
	x <- 5 + (chr -1) * 12 + conum
	y <- 100 - pos
	format_co_count[i,] <- c(chr,x,y)
}

format_gene_count <- matrix(0,nrow=nrow(gene_count),ncol=3)
for (i in 1:nrow(gene_count)) {
	chr <- gene_count[i,1]
	pos <- (gene_count[i,2] / max_chrlen) * rel_len
	conum <- (gene_count[i,3] / max_gene_count) * 7.5
	x <- 5 + (chr -1) * 12 + conum
	y <- 100 - pos
	format_gene_count[i,] <- c(chr,x,y)
}


format_F1_count <- matrix(0,nrow=nrow(gmapF1_cM),ncol=3)
for (i in 1:nrow(gmapF1_cM)) {
	chr <- gmapF1_cM[i,1]
	pos <- (gmapF1_cM[i,2] / max_chrlen) * rel_len
	conum <- (gmapF1_cM[i,3] / max_co_count2) * 7.5
	x <- 5 + (chr -1) * 12 + conum
	y <- 100 - pos
	format_F1_count[i,] <- c(chr,x,y)
}

for (gm in 1:ncol(phasing) ) {
target_gmt <- colnames(phasing)[gm]
output <- paste("results/",target_gmt,".co",".pdf",sep="")

pdf(output, width=8.53, height=3.8)
par(mar=c(0,0,0,0))
plot (0, axes=F, type="n", xlim=c(-10,180), ylim=c(0,118), xlab="", ylab="")

for (i in 1:15) {
	tmp <- format_co_count[format_co_count[,1] == i,]
	#tmpF1 <- format_F1_count[format_F1_count[,1] == i,]
	#tmpgene <- format_gene_count[format_gene_count[,1] == i,]
	lines(x=tmp[,2],y=tmp[,3], col="#9EB7BB", lwd=1)
	#lines(x=tmpF1[,2],y=tmpF1[,3], col="#DEB774", lwd=1)
	#lines(x=tmpgene[,2],y=tmpgene[,3],col="#DEB774", lwd=1)
	
	# plot co bars #
	cobar_xlft <- 5 + (i-1) * 12
	cobar_xright <- cobar_xlft + 7.5
	lines(x=c(cobar_xlft,cobar_xright),y=c(101.5,101.5))
	lines(x=c(cobar_xlft,cobar_xlft),y=c(101.5,102))
	lines(x=c(cobar_xright,cobar_xright),y=c(101.5,102))
	lines(x=c((cobar_xlft+cobar_xright)/2,(cobar_xlft+cobar_xright)/2),y=c(101.5,102))
	text(x=cobar_xlft,y=104.5,labels=0,cex=0.6)
	text(x=cobar_xright,y=104.5,labels=round(max_co_count1,0),cex=0.6)
	text(x=(cobar_xlft+cobar_xright)/2,y=104.5,labels=round(max_co_count1/2,1),cex=0.6)
	
	text(x=cobar_xlft,y=108.5,labels=0.4,cex=0.6)
	text(x=cobar_xright,y=108.5,labels=1,cex=0.6)
	text(x=(cobar_xlft+cobar_xright)/2,y=108.5,labels=0.7,cex=0.6)
	
	chr_pos <- (cobar_xlft + cobar_xright) / 2
	text (chr_pos,116,labels=paste("chr",i,sep=""), cex = 0.7)
	
	if (i == 1) {
		text(x=cobar_xlft-2,y=104.5,labels="single sperm:",cex=0.7,adj=1, col="#9EB7BB")
		text(x=cobar_xlft-2,y=108.5,labels="frequency of bin:",cex=0.7,adj=1, col="#DEB774")
	}
	
	# plot hot co spot #
#	tmp_co_count <- co_count[format_co_count[,1] == i,]
#	hot1 <- tmp[tmp_co_count[,3] > 1.5,]
#	if (length(nrow(hot1)) != 0) {
#		points(x=rep(cobar_xright+0.5,nrow(hot1)), y=hot1[,3], pch=17, cex=0.6, col="#9EB7BB")
#	}else{
#		points(x=cobar_xright+0.5, y=hot1[3], pch=17, cex=0.6, col="#9EB7BB")
#	}
	
	#tmp_co_count2 <- gmapF1_cM[format_F1_count[,1] == i,]
	#hot2 <- tmpF1[tmp_co_count2[,3] > 1.5,]
	#if (length(nrow(hot2)) != 0) {
	#	points(x=rep(cobar_xright+1.5,nrow(hot2)), y=hot2[,3], pch=17, cex=0.6, col="#DEB774")
	#}else{
	#	points(x=cobar_xright+1.5, y=hot2[3], pch=17, cex=0.6, col="#DEB774")
	#}	
	
	# plot phasing #
	linkage_file = paste("linkage_group/","chr",i,".map.ordered.txt",sep="")
	gmap <- read.table(linkage_file,header=F)
	col_vec <- phasing[gmap[,1],target_gmt]
	col_vec <- sub("A","#2E5A71",col_vec)
	col_vec <- sub("B","#9EB7BB",col_vec)
	col_vec <- sub("U","#FFFFFF",col_vec)
	
	# plot freq #
	tmp_ph_gmap <- phasing[gmap[,1],]
	freq_vec <- rep(0,nrow(tmp_ph_gmap))
	for (fq in 1:nrow(tmp_ph_gmap)){
		freqA <- sum(phasing[fq,] == 'A') / (sum(phasing[fq,] == 'A') + sum(phasing[fq,] == 'B')) - 0.4
		freq_vec[fq] <- (freqA / 0.6) * 7.5 + 5 + (i-1) * 12
	}
	
	phasing_plot_mat <- binmap[gmap[,1],3:4]
	phasing_plot_mat <- 100 - (phasing_plot_mat/ max_chrlen) * rel_len
	phasing_plot_mat_pos <- rep(0,nrow(phasing_plot_mat))
	for (pm in 1:nrow(phasing_plot_mat)) {
		chr_lft <- 5 + (i-1) * 12
		rect(chr_lft-1.5,phasing_plot_mat[pm,2],chr_lft-0.1,phasing_plot_mat[pm,1], border=NA, col=col_vec[pm])
		phasing_plot_mat_pos[pm] <- (phasing_plot_mat[pm,2] + phasing_plot_mat[pm,1]) / 2
	}
	
	lines(x=freq_vec,y=phasing_plot_mat_pos, col="#DEB774")
	
	# plot chr #
	#chr_lft <- 5 + (i-1) * 12
	#chr_bot <- (chr_len[i,2] / max_chrlen) * rel_len
	#rect(chr_lft-1.5,100-chr_bot,chr_lft-0.1,100, border=NA, col= "#2E5A71")
}

## plot left axis ##
lft_axis <- seq(0,max_chrlen,length=10)
lines(x=c(-1,-1),y=c(100,10))
for (i in 1:length(lft_axis)) {
	
	rel_pos <- 100 - (lft_axis[i] / max_chrlen) * rel_len
	
	lines(x=c(-1,-1.5),y=c(rel_pos,rel_pos))
	mb <- round(lft_axis[i] / 1000000,1)
	text(x=-3,y=rel_pos,labels=mb,cex=0.7,adj=1)
}

## plot bars ##
rect(5 + 3*12, 20, 5 + 3*12 + 4, 22, border="#2E5A71", col="#2E5A71")
text(5 + 3*12 + 4+1, y= 21, labels="hapA", adj=0, cex = 0.7)
rect(5 + 3*12 + 15, 20, 5 + 3*12 + 4 + 15, 22, border="#9EB7BB", col="#9EB7BB")
text(5 + 3*12 + 4 + 15+1, y= 21, labels="hapB", adj=0, cex = 0.7)

lines(x=c(5 + 3*12 + 15 + 15, 5 + 3*12 + 15 + 15 + 4), y=c(21,21), lwd=2, col="#9EB7BB")
text(5 + 3*12 + 15 + 15 + 4 + 1, y=21, labels="single sperm (cM/Mb)", adj=0, cex=0.7)

#lines(x=c(5 + 3*12 + 15 + 15 + 35, 5 + 3*12 + 15 + 15 + 35 + 4), y=c(21,21), lwd=2, col="#DEB774")
#text(5 + 3*12 + 15 + 15 + 35 + 4 + 1, y=21, labels="F1 (cM/Mb)", adj=0, cex=0.7)

#points(x=5 + 3*12 + 15 + 15 + 35, y=21, pch=17, col="#9EB7BB")
#text(5 + 3*12 + 15 + 15 + 35 + 2, y=21, labels="single sperm co hot region", adj=0, cex=0.7)

#points(x=5 + 3*12 + 15 + 15 + 35 + 30, y=21, pch=17, col="#DEB774")
#text(5 + 3*12 + 15 + 15 + 35 + 30 + 2, y=21, labels="F1 co hot region", adj=0, cex=0.7)

lines(x=c(5 + 3*12 + 15 + 15 + 35, 5 + 3*12 + 15 + 15 + 35 + 4), y=c(21,21), lwd = 2, col="#DEB774")
text(5 + 3*12 + 15 + 15 + 35 + 4 + 1, y=21, labels="frequency of bins", adj=0, cex=0.7)

text(x=-12,y=59, srt=90 ,labels="Physical position (Mb)", cex = 0.7)

dev.off()

}
