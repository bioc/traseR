
################################################
### plot SNP functional class distribution

plotContext<-function(snpdb,region=NULL,keyword=NULL,pvalue=1e-3){
	
	snpdb=snpdb[snpdb$p.value<=pvalue]
	
	if(is.null(keyword)){
		snp=unique(as.data.frame(snpdb[,c("SNP_ID","Context")]))
		snp.meta=snp[,setdiff(colnames(snp),c("seqnames","start","end","strand","width"))]
		snp=makeGRangesFromDataFrame(snp[,c("seqnames","start","end","strand")])
		mcols(snp)=snp.meta
	}else{
		ind=grep(tolower(keyword),tolower(snpdb$Trait))
		if(length(ind)==0){
			stop("No match key word found in Trait!")
		}else{
			snp=unique(as.data.frame(snpdb[ind,c("SNP_ID","Trait","Context")]))
			snp.meta=snp[,setdiff(colnames(snp),c("seqnames","start","end","strand","width"))]
			snp=makeGRangesFromDataFrame(snp[,c("seqnames","start","end","strand")])
			mcols(snp)=snp.meta
		}
	}
	if(is.null(region)){
		a=snp$Context
	}else{
		if(class(region)!="data.frame" & class(region)!="GRanges")
			stop("Query genomic interval should be either data.frame or GRange!")
		if(class(region)=="data.frame"){
			if(!identical(colnames(region),c("chr","start","end")))
				stop("Interval colnames should be chr,start,end!")
			region=makeGRangesFromDataFrame(region)
		}
		seqlevel=c(paste("chr",1:22,sep=""),"chrX")
		region=region[!is.na(match(seqnames(region),seqlevel))]
		end(snp)=end(snp)+1
		o=findOverlaps(region,snp)
		if(length(o@subjectHits)==0)
			stop("No SNP overlapped with the query region!")
		a=snp$Context[unique(o@subjectHits)]
	}
	a=table(a)
	pa=round(a/sum(a)*100)
	lbls=paste(names(a), pa)
	lbls=paste(lbls,"%",sep="")
	pie(pa,labels=lbls,col=rainbow(length(lbls)),main="Pie Chart of Context")
}

#####################################################
### plot GWAS pvalues density distribution for snp

plotPvalue<-function(snpdb,region=NULL,keyword=NULL,plot.type=c("densityplot","boxplot"),pvalue=1e-3,xymax=50){

	snpdb=snpdb[snpdb$p.value<=pvalue & snpdb$p.value!=0]
	allpval=-log(snpdb$p.value,base=10)

	if(!is.null(keyword)){
		ind=grep(tolower(keyword),tolower(snpdb$Trait))
		if(length(ind)==0){
			stop("No match key word found in Trait!")
		}else{
			snp=snpdb[ind]
		}
		pval1=unique(allpval[ind])
		if(is.null(region)){
			if(plot.type=="boxplot"){
				plotdata=list(allpval,pval1)
				names(plotdata)=c("SNP(all)","SNP(keyword)")
				boxplot(plotdata,xlab="-log10Pvalue",main="-log10Pvalue Distribution",col=c("yellow","red"),ylim=c(0,xymax))
			}else if(plot.type=="densityplot"){
				d=density(allpval)
				d1=density(pval1)
				ymax=max(d$y,d1$y)
				plot(density(allpval),xlab="",ylab="","-log10Pvalue Distribution",col="yellow",lwd=2,xlim=c(0,xymax))
				lines(d1,col="red",lwd=2)
				legend("topright",legend=c("SNP(all)","SNP(keyword)"),lty=c(1,1),lwd=2,col=c("yellow","red"))
			}else{
				stop("Undefined plot type!")
			}
		}else{
			if(class(region)!="data.frame" & class(region)!="GRanges")
				stop("Query genomic interval should be either data.frame or GRange!")
			if(class(region)=="data.frame"){
				if(!identical(colnames(region),c("chr","start","end")))
					stop("Interval colnames should be chr,start,end!")
				region=makeGRangesFromDataFrame(region)
			}
			seqlevel=c(paste("chr",1:22,sep=""),"chrX")
			region=region[!is.na(match(seqnames(region),seqlevel))]
			end(snp)=end(snp)+1
			o=findOverlaps(region,snp)
			if(length(o@subjectHits)==0)
				stop("No SNP overlapped with genomic intervals!")
			pval2=unique(allpval[o@subjectHits])		
			if(plot.type=="boxplot"){
				plotdata=list(allpval,pval1,pval2)
				names(plotdata)=c("SNP(all)","SNP(keyword)","SNP(overlapped)")
				boxplot(plotdata,xlab="-log10Pvalue",main="-log10Pvalue Distribution",col=c("yellow","red","blue"),ylim=c(0,xymax))
			}else if(plot.type=="densityplot"){
				d=density(allpval)
				d1=density(pval1)
				d2=density(pval2)
				ymax=max(d$y,d1$y,d2$y)
				plot(d,xlab="",ylab="","-log10Pvalue Distribution",col="yellow",lwd=2,xlim=c(0,xymax))
				lines(d1,col="red",lwd=2)
				lines(d2,col="blue",lwd=2)
				legend("topright",legend=c("SNP(all)","SNP(keyword)","SNP(overlapped)"),lty=c(1,1),lwd=2,col=c("yellow","red","blue"))

			}else{
				stop("Undefined plot type!")
			}
		}	
	}else{
		if(is.null(region)){
			if(plot.type=="boxplot"){
				boxplot(allpval,xlab="-log10Pvalue",main="-log10Pvalue Distribution",col="yellow",ylim=c(0,xymax))
			}else if(plot.type=="densityplot"){
				plot(density(allpval),xlab="",ylab="","-log10Pvalue Distribution",col="yellow",lwd=2,xlim=c(0,xymax))
			}else{
				stop("Undefined plot type!")
			}
		}else{
			if(class(region)!="data.frame" & class(region)!="GRanges")
				stop("Query genomic interval should be either data.frame or GRange!")
			if(class(region)=="data.frame"){
				if(!identical(colnames(region),c("chr","start","end")))
					stop("Interval colnames should be chr,start,end!")
				region=makeGRangesFromDataFrame(region)
			}
			seqlevel=c(paste("chr",1:22,sep=""),"chrX")
			region=region[!is.na(match(seqnames(region),seqlevel))]
			end(snp)=end(snp)+1
			o=findOverlaps(region,snp)
			if(length(o@subjectHits)==0)
				stop("No SNP overlapped with the genomic intervals!")
			pval2=unique(allpval[o@subjectHits])	
			if(plot.type=="boxplot"){
				plotdata=list(allpval,pval2)
				names(plotdata)=c("SNP(all)","SNP(Overlapped)")
				boxplot(plotdata,xlab="-log10Pvalue",main="-log10Pvalue Distribution",col=c("yellow","blue"),ylim=c(0,xymax))
			}else if(plot.type=="densityplot"){
				d=density(allpval)
				d2=density(pval2)
				ymax=max(d$y,d2$y)
				plot(d,xlab="",ylab="","-log10Pvalue Distribution",col="yellow",lwd=2,xlim=c(0,xymax))
				lines(d2,col="blue",lwd=2)
				legend("topright",legend=c("SNP(all)","SNP(overlapped)"),lty=c(1,1),lwd=2,col=c("yellow","blue"))
			}else{
				stop("Undefined plot type!")
			}
		}
	}
}




#################################################################
### plot gene along with possible functional snp on chromosome

plotGene<-function(snpdb,gene,ext=10000){
	if(is.null(gene)){
		stop("Must specify one gene!")
	}
	if(length(gene)!=1){
		stop("Only one gene is allowed for display!")
	}
	
	ind=grep(tolower(gene),tolower(snpdb$GENE_NAME))
	if(length(ind)==0){
		stop("No matched gene is found!")
	}else{
		snp=unique(as.data.frame(snpdb[ind,c("SNP_ID","p.value","Trait")]))
		snp.meta=snp[,setdiff(colnames(snp),c("seqnames","start","end","strand","width"))]
		snp=makeGRangesFromDataFrame(snp[,c("seqnames","start","end","strand")])
		mcols(snp)=snp.meta
		gene.loci=snpdb[ind,c("GENE_NAME","GENE_START","GENE_END","GENE_STRAND")]	
		gene.loci=unique(as.data.frame(gene.loci))
		gene.loci=gene.loci[,setdiff(colnames(gene.loci),c("start","end","strand","width"))]
	}
	
	min.pos=gene.loci$GENE_START-ext
	max.pos=gene.loci$GENE_END+ext
	center.pos=round((min.pos+max.pos)/2)
	
	ind0=snp$p.value==0
	ind1=snp$p.value!=0
	if(sum(ind1)>0 & sum(ind0)==0){
		maxp=max(-log10(snp$p.value))
	}else if(sum(ind1)>0 & sum(ind0)>0){
		maxp=max(-log10(snp$p.value[ind1]))+2
	}else if(sum(ind1)==0 & sum(ind0)>0){
		maxp=10
	}
	
	plot(0,xlim=c(min.pos, max.pos),ylim=c(0,maxp),xlab="", ylab="", axes=FALSE)
	mtext(paste("Chr",gene.loci$seqnames, " position (bp)", sep=""), side=1, line=2.5)
	axis(1, at=c(min.pos,center.pos,max.pos), labels=c(min.pos,center.pos,max.pos), las=1) 
	axis(2, at=seq(0,maxp,2), labels=seq(0,maxp,2), las=1) 
	mtext("-log10Pvalue", side=2, line=2)
	box()

	if(sum(ind1)>0 & sum(ind0)==0){
		points(start(snp), -(log10(snp$p.value)), pch=23, cex=2, bg="red")
		text(start(snp), -(log10(snp$p.value))-0.5, labels=c(paste(as.matrix(snp$SNP_ID),"(",snp$Trait,")",sep="")), pos=3, offset=1,cex=1)
	}else if(sum(ind1)>0 & sum(ind0)>0){
		points(start(snp)[ind1], -(log10(snp$p.value[ind1])), pch=23, cex=2, bg="red")
		text(start(snp)[ind1], -(log10(snp$p.value[ind1]))-0.5, labels=c(paste(as.matrix(snp$SNP_ID[ind1]),"(",snp$Trait[ind1],")",sep="")), pos=3, offset=1,cex=0.5)
		points(start(snp)[ind0], maxp+1, pch=23, cex=2, bg="red")
		text(start(snp)[ind0], maxp+0.5, labels=c(paste(as.matrix(snp$SNP_ID[ind0]),"(",snp$Trait[ind0],")",sep="")), pos=3, offset=1,cex=1)
	}else if(sum(ind1)==0 & sum(ind0)>0){
		points(start(snp),1, 9, pch=23, cex=2, bg="red")
		text(start(snp), 9-0.5, labels=c(paste(as.matrix(snp$SNP_ID),"(",snp$Trait,")",sep="")), pos=3, offset=1,cex=1)
	}
	
	if(gene.loci$GENE_STRAND == "+"){
		arrows(max(gene.loci$GENE_START, min.pos),0, min(gene.loci$GENE_END, max.pos),0, length=0.05, lwd=2, code=2, lty="solid", col="darkgreen")
	}else if(gene.loci$GENE_STRAND == "-"){		
		arrows(max(gene.loci$GENE_START, min.pos),0, min(gene.loci$GENE_END, max.pos),0, length=0.05, lwd=2, code=1, lty="solid", col="darkgreen")
	}else{
		segments(max(gene.loci$GENE_START, min.pos),0, min(gene.loci$GENE_END, max.pos),0,lwd=2,lty="solid", col="darkgreen")
	}
	if(!is.na(gene.loci$GENE_NAME)) {
		text((gene.loci$GENE_START+gene.loci$GENE_END)/2, 0.2, labels=gene.loci$GENE_NAME, cex=1)
	}
}





#############################################################
### plot SNP along with possible nearby genes on chromosome


plotSNP<-function(snpdb,snpid,ext=10000){
	ind=grep(snpid,snpdb$SNP_ID)
	if(length(ind)==0){
		stop(paste(snpid,"cannot be found"))	
	}else{
		snp=unique(as.data.frame(snpdb[ind,c("SNP_ID","p.value","Trait")]))
		snp.meta=snp[,setdiff(colnames(snp),c("seqnames","start","end","strand","width"))]
		snp=makeGRangesFromDataFrame(snp[,c("seqnames","start","end","strand")])
		mcols(snp)=snp.meta
		gene.loci=snpdb[ind,c("GENE_NAME","GENE_START","GENE_END","GENE_STRAND")]	
		gene.loci=unique(as.data.frame(gene.loci))
		gene.loci=gene.loci[,setdiff(colnames(gene.loci),c("start","end","strand","width"))]
	}
		
	min.pos=min(start(snp))-ext
	max.pos=max(start(snp))+ext
	center.pos=round((min.pos+max.pos)/2)
	
	ind0=snp$p.value==0
	ind1=snp$p.value!=0
	if(sum(ind1)>0 & sum(ind0)==0){
		maxp=max(-log10(snp$p.value))
	}else if(sum(ind1)>0 & sum(ind0)>0){
		maxp=max(-log10(snp$p.value[ind1]))+2
	}else if(sum(ind1)==0 & sum(ind0)>0){
		maxp=10
	}
	
	plot(0,xlim=c(min.pos, max.pos),ylim=c(0,maxp),xlab="", ylab="", axes=FALSE)
	mtext(paste("Chr",gene.loci$seqnames, " position (bp)", sep=""), side=1, line=2.5)
	axis(1, at=c(min.pos,center.pos,max.pos), labels=c(min.pos,center.pos,max.pos), las=1) 
	axis(2, at=seq(0,maxp,2), labels=seq(0,maxp,2), las=1) 
	mtext("-log10Pvalue", side=2, line=2)
	box()

	if(sum(ind1)>0 & sum(ind0)==0){
		points(start(snp), -(log10(snp$p.value)), pch=23, cex=2, bg="red")
		text(start(snp), -(log10(snp$p.value))-0.5, labels=c(paste(as.matrix(snp$SNP_ID),"(",snp$Trait,")",sep="")), pos=3, offset=1,cex=1)
	}else if(sum(ind1)>0 & sum(ind0)>0){
		points(start(snp)[ind1], -(log10(snp$p.value[ind1])), pch=23, cex=2, bg="red")
		text(start(snp)[ind1], -(log10(snp$p.value[ind1]))-0.5, labels=c(paste(as.matrix(snp$SNP_ID[ind1]),"(",snp$Trait[ind1],")",sep="")), pos=3, offset=1,cex=0.5)
		points(start(snp)[ind0], maxp+1, pch=23, cex=2, bg="red")
		text(start(snp)[ind0], maxp+0.5, labels=c(paste(as.matrix(snp$SNP_ID[ind0]),"(",snp$Trait[ind0],")",sep="")), pos=3, offset=1,cex=1)
	}else if(sum(ind1)==0 & sum(ind0)>0){
		points(start(snp),1, 9, pch=23, cex=2, bg="red")
		text(start(snp), 9-0.5, labels=c(paste(as.matrix(snp$SNP_ID),"(",snp$Trait,")",sep="")), pos=3, offset=1,cex=1)
	}
	
	if(gene.loci$GENE_STRAND == "+"){
		arrows(max(gene.loci$GENE_START, min.pos),0, min(gene.loci$GENE_END, max.pos),0, length=0.05, lwd=2, code=2, lty="solid", col="darkgreen")
	}else if(gene.loci$GENE_STRAND == "-"){		
		arrows(max(gene.loci$GENE_START, min.pos),0, min(gene.loci$GENE_END, max.pos),0, length=0.05, lwd=2, code=1, lty="solid", col="darkgreen")
	}else{
		segments(max(gene.loci$GENE_START, min.pos),0, min(gene.loci$GENE_END, max.pos),0,lwd=2,lty="solid", col="darkgreen")
	}
	if(!is.na(gene.loci$GENE_NAME)) {
		text((gene.loci$GENE_START+gene.loci$GENE_END)/2, 0.2, labels=gene.loci$GENE_NAME, cex=1)
	}
	
}



######################################################################
### plot possible functional snp & genes for a given genomic interval

plotInterval<-function(snpdb,interval,ext=10000){
	
	interval=makeGRangesFromDataFrame(interval)
	end(snpdb)=end(snpdb)+1
	o=findOverlaps(interval,snpdb)
	if(length(o@subjectHits)==0)
		stop("The query region is not overlapped with any function SNP")
	
	snp=unique(as.data.frame(snpdb[unique(o@subjectHits),c("SNP_ID","p.value","Trait")]))
	snp.meta=snp[,setdiff(colnames(snp),c("seqnames","start","end","strand","width"))]
	snp=makeGRangesFromDataFrame(snp[,c("seqnames","start","end","strand")])
	mcols(snp)=snp.meta
	genes=snpdb[unique(o@subjectHits),c("GENE_NAME","GENE_START","GENE_END","GENE_STRAND")]	
	genes=unique(as.data.frame(genes))
	genes=genes[,setdiff(colnames(genes),c("start","end","strand","width"))]

	min.pos=start(interval)-ext
	max.pos=end(interval)+ext
	center.pos=round((min.pos+max.pos)/2)
	ind0=snp$p.value==0
	ind1=snp$p.value!=0
	if(sum(ind1)>0 & sum(ind0)==0){
		maxp=max(-log10(snp$p.value))
	}else if(sum(ind1)>0 & sum(ind0)>0){
		maxp=max(-log10(snp$p.value[ind1]))+2
	}else if(sum(ind1)==0 & sum(ind0)>0){
		maxp=10
	}
	
	plot(0,xlim=c(min.pos, max.pos),ylim=c(0,maxp),xlab="", ylab="", axes=FALSE)
	mtext(paste("Chr",as.character(seqnames(interval)), " position (bp)", sep=""), side=1, line=2.5)
	axis(1, at=c(min.pos,center.pos,max.pos), labels=c(min.pos,center.pos,max.pos), las=1) 
	axis(2, at=seq(0,maxp,2), labels=seq(0,maxp,2), las=1) 
	mtext("-log10Pvalue", side=2, line=2)
	box()

	if(sum(ind1)>0 & sum(ind0)==0){
		points(start(snp), -(log10(snp$p.value)), pch=23, cex=2, bg="red")
		text(start(snp), -(log10(snp$p.value))-0.5, labels=c(paste(as.matrix(snp$SNP_ID),"(",snp$Trait,")",sep="")), pos=3, offset=1,cex=1)
	}else if(sum(ind1)>0 & sum(ind0)>0){
		points(start(snp)[ind1], -(log10(snp$p.value[ind1])), pch=23, cex=2, bg="red")
		text(start(snp)[ind1], -(log10(snp$p.value[ind1]))-0.5, labels=c(paste(as.matrix(snp$SNP_ID[ind1]),"(",snp$Trait[ind1],")",sep="")), pos=3, offset=1,cex=0.5)
		points(start(snp)[ind0], maxp+1, pch=23, cex=2, bg="red")
		text(start(snp)[ind0], maxp+0.5, labels=c(paste(as.matrix(snp$SNP_ID[ind0]),"(",snp$Trait[ind0],")",sep="")), pos=3, offset=1,cex=1)
	}else if(sum(ind1)==0 & sum(ind0)>0){
		points(start(snp),1, 9, pch=23, cex=2, bg="red")
		text(start(snp), 9-0.5, labels=c(paste(as.matrix(snp$SNP_ID),"(",snp$Trait,")",sep="")), pos=3, offset=1,cex=1)
	}
	
	for(i in 1:nrow(genes)){ 
		if(genes[i,]$GENE_STRAND == "+"){
			arrows(max(genes[i,]$GENE_START, min.pos),0, min(genes[i,]$GENE_END, max.pos),0, length=0.05, lwd=2, code=2, lty="solid", col="darkgreen")
		}else if(genes[i,]$GENE_STRAND == "-"){		
			arrows(max(genes[i,]$GENE_START, min.pos),0, min(genes[i,]$GENE_END, max.pos),0, length=0.05, lwd=2, code=1, lty="solid", col="darkgreen")
		}else{
			segments(max(genes[i,]$GENE_START, min.pos),0, min(genes[i,]$GENE_END, max.pos),0,lwd=2,lty="solid", col="darkgreen")
		}
		if(!is.na(genes[i,]$GENE_NAME)) {
			text((genes[i,]$GENE_START+genes[i,]$GENE_END)/2, 0.2, labels=genes[i,]$GENE_NAME, cex=1)
		}
	}

}





