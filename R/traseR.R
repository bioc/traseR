################################################################
#### Core function to perform hypothesis testing.
###  It is a nested function of enrichTest, no export for users.

enrichTest<-function(params,flag){

		if(flag){
			tb=matrix(c(params$hits,params$nsnp-params$hits,params$region.size-params$hits,
				params$genome.size-params$region.size-params$nsnp+params$hits),2,2)
		}else{
			tb=matrix(c(params$hits1,params$hits2,params$nsnp1-params$hits1,params$nsnp2-params$hits2),2,2)
		}
		rownames(tb)=c("region.in","region.out")
		colnames(tb)=c("taSNPs","non-taSNPs")
		
		odds.ratio=((1+tb[1,1])/(1+tb[1,ncol(tb)]))/((1+tb[nrow(tb),1])/(1+tb[nrow(tb),ncol(tb)]))	
		
		
		if(params$test.method=="nonparametric"){	
			if(params$alternative=="greater"){
				pvalue=sum(params$hits<=params$crit)/params$nbatch/params$ntimes
			}else if(params$alternative=="less"){
				pvalue=sum(params$hits>=params$crit)/params$nbatch/params$ntimes
			}else{
				stop("nonparametric test only support one.sided test!")
			}
		}else if(params$test.method=="chisq"){
			a=chisq.test(tb)
			pvalue=a$p.value
		}else if(params$test.method=="binomial"){
			if(flag){
				a=binom.test(params$hits,params$region.size,
					p=params$nsnp/params$genome.size,alternative=params$alternative)
			}else{
				a=binom.test(params$hits1+1, params$hits1+params$hits2+1,
					p=params$nsnp1/(params$nsnp1+params$nsnp2),alternative=params$alternative)
			}
			pvalue=a$p.value
		}else if(params$test.method=="fisher"){
			if(params$alternative=="two.sided") stop("Fisher' exact test only support one.sided test!")
			if(flag){
				pvalue= phyper(params$hits,params$nsnp,params$genome.size-params$nsnp,params$region.size,
						lower.tail=params$alternative=="less")
			}else{
				pvalue= phyper(params$hits1,params$nsnp1,params$nsnp2,params$hits2+params$hits1,
						lower.tail=params$alternative=="less")
			}
		}
		list(tb=tb,pvalue=pvalue,odds.ratio=odds.ratio)
}


##########################################################################
#### Core function of traseR, perform SNP enrichment Test for each trait


traseR<-function(snpdb,region,snpdb.bg=NULL,keyword=NULL,
	rankby=c("pvalue","odds.ratio"),	
	test.method=c("binomial","fisher","chisq","nonparametric"),
	alternative = c("greater","less","two.sided"),
	ntimes=100, nbatch=1,
	trait.threshold=0,traitclass.threshold=0,pvalue=1e-3){
	
	test.method=match.arg(test.method)
	alternative=match.arg(alternative)
	rankby=match.arg(rankby)
	
	if(missing(snpdb))
		stop("Functional SNP database must be specified!")	
	if(missing(region))
		stop("Need to specify the genomic intervals!")
	if(!is.null(snpdb$p.value)){
		snpdb=snpdb[snpdb$p.value<=pvalue]
	}
	flag=is.null(snpdb.bg)	
		
	### only autosomal chromosome and X chromosome regions
	if(class(region)!="data.frame" & class(region)!="GRanges")
		stop("Query genomic interval should be either data.frame or GRanges!")
	if(class(region)=="data.frame"){
		if(!identical(colnames(region),c("chr","start","end")))
			stop("Interval colnames should be chr,start,end!")
		region=makeGRangesFromDataFrame(region)
	}
	
	seqlevel=c(paste("chr",1:22,sep=""),"chrX")
	region=region[!is.na(match(seqnames(region),seqlevel))]
	nregion=length(region)
	region.size=sum(as.numeric(end(region))-as.numeric(start(region)))

	### region size and hg19 genome size
	seqlength=seqlengths(BSgenome.Hsapiens.UCSC.hg19)
	seqlength=seqlength[match(seqlevel,names(seqlength))]
	genome.size=sum(as.numeric(seqlength))
	
	
	### load SNP database in NHGRI format,calculate overlap with query region
	#trait-specific SNP dataframe
	snp=unique(as.data.frame(snpdb[,c("SNP_ID","Trait")]))
	snp.meta=snp[,setdiff(colnames(snp),c("seqnames","start","end","strand","width"))]
	snp=makeGRangesFromDataFrame(snp[,c("seqnames","start","end","strand")])
	mcols(snp)=snp.meta
	
	#trait-class-specific SNP dataframe
	snpclass=unique(as.data.frame(snpdb[,c("SNP_ID","Trait_Class")]))
	snpclass.meta=snpclass[,setdiff(colnames(snpclass),c("seqnames","start","end","strand","width"))]
	snpclass=makeGRangesFromDataFrame(snpclass[,c("seqnames","start","end","strand")])
	mcols(snpclass)=snpclass.meta
	
	
	if(!is.null(keyword)){
		ind=grep(tolower(keyword),tolower(snp$Trait))
		if(length(ind)==0)
			stop("No match key word found in Trait!")
		snp=snp[ind]
		ind=grep(tolower(keyword),tolower(snp$Trait_Class))
		if(length(ind)==0)
			stop("No match key word found in Trait Class!")
		snpclass=snpclass[ind]
	}
	
	### trait-specific snp GRanges
	snp=snp[order(snp$Trait)]
	traits=unique(snp$Trait)
	traits=unique(names(table(snp$Trait))[table(snp$Trait)>trait.threshold])
	snp=snp[!is.na(match(snp$Trait,traits))]
	
	o=findOverlaps(region,snp)	
	allhits=length(unique(o@subjectHits))
	if(allhits==0){
		stop("Overall SNPs are not overlapped with any genomic intervals!")
	}
	nsnp=length(unique(snp$SNP_ID))
	tt=snp[unique(o@subjectHits)]
	ntraits=length(traits)
	
	
	### trait-class-specific snp GRanges
	snpclass=snpclass[order(snpclass$Trait_Class)]
	traitclass=unique(snpclass$Trait_Class)
	traitclass=unique(names(table(snpclass$Trait_Class))[table(snpclass$Trait_Class)>traitclass.threshold])
	snpclass=snpclass[!is.na(match(snpclass$Trait_Class,traitclass))]
	
	o=findOverlaps(region,snpclass)	
	allhits=length(unique(o@subjectHits))
	if(allhits==0){
		stop("Overall SNPs are not overlapped with any genomic intervals!")
	}
	nsnpclass=length(unique(snpclass$SNP_ID))
	ttclass=snpclass[unique(o@subjectHits)]
	ntraitclass=length(traitclass)
	
	
	### print out the size of region and number of traits associated with the region
	message(paste("There are ",region.size,"bp in the query region, "),paste("accounting for ",region.size/genome.size," of the genome."))
	message(paste("There are ",ntraits,"traits in the analysis."))
	message(paste("There are ",ntraitclass,"trait class in the analysis."))

	
	### Calculate background SNP overlap
	if(!flag){
 		ind=match(snpdb.bg$SNP_ID,snp$SNP_ID)
		snp.bg=unique(snpdb.bg[is.na(ind)])
		o=findOverlaps(region,snp.bg)	
		allhits.bg=length(unique(o@subjectHits))
		nsnp.bg=length(unique(snp.bg$SNP_ID))
	}


	### parameters pass to enrichTest
	params=list()	
	params$region.size=region.size
	params$genome.size=genome.size
	params$test.method=test.method
	params$alternative=alternative
	
	
	### performance criteria
	pvalues=odds.ratios=hits=snpnums=rep(-1,ntraits)
	pvalues.class=odds.ratios.class=hits.class=snpnums.class=rep(-1,ntraitclass)
	booleans=olist=booleans.class=olist.class=list()
	crit=crit.class=numeric(ntimes*nbatch)


	### if it is nonparametric test
	if(test.method=="nonparametric"){
		len.order=seqlength[(match(seqnames(region),seqlevel))@values]
		len.order=rep(len.order,ntimes)
		rchr=rep(seqnames(region),ntimes)
		region.width=end(region)-start(region)
		region.width=rep(region.width,ntimes)
		for(i in seq_len(nbatch)){
			rstart=floor((len.order-region.width)*runif(length(len.order)))
			rend=rstart+region.width		
			region.permute=GRanges(seqnames = Rle(rchr), ranges = IRanges(start=rstart,end=rend))
			o=findOverlaps(region,region.permute)
			region.permute=region.permute[-o@subjectHits]
			o=findOverlaps(region.permute,snp)
			oclass=findOverlaps(region.permute,snpclass)
			for(j in seq_len(ntimes)){
				booleans[[(i-1)*ntimes+j]]=o@queryHits<=j*nregion & o@queryHits>=(j-1)*nregion
				crit[(i-1)*ntimes+j]=length(unique(o@subjectHits[booleans[[(i-1)*ntimes+j]]]))
				booleans.class[[(i-1)*ntimes+j]]=oclass@queryHits<=j*nregion & oclass@queryHits>=(j-1)*nregion
				crit.class[(i-1)*ntimes+j]=length(unique(oclass@subjectHits[booleans.class[[(i-1)*ntimes+j]]]))
			}
			olist[[i]]=o
			olist.class[[i]]=oclass
		}
		params$ntimes=ntimes
		params$nbatch=nbatch
		params$booleans=booleans
		params$crit=crit
		params$booleans.class=booleans.class
		params$crit.class=crit.class
	}
	
	if(flag){
		### test overall trait-associated SNP enrichment
		params$nsnp=nsnp
		params$hits=allhits
		b=enrichTest(params,flag)
		tb.all=data.frame(Trait="All",p.value=b$pvalue,odds.ratio=b$odds.ratio,hits=allhits,taSNP.num=nsnp)
		colnames(tb.all)=c("Trait","p.value","odds.ratio","taSNP.hits","taSNP.num")
		rownames(tb.all)=NULL
		
		#trait-specific testing
		for(i in seq_len(ntraits)){
			if(i%%100==0) message(i," traits have been tested!")
			hits[i]=sum(tt$"Trait"==traits[i])
			snpnums[i]=sum(snp$"Trait"==traits[i])
			params$nsnp=snpnums[i]
			params$hits=hits[i]
			if(test.method=="nonparametric"){
				for(j in seq_len(nbatch)){
					for(k in seq_len(ntimes)){
						crit[(j-1)*ntimes+k]=length(unique(  olist[[j]]@subjectHits[ booleans[[(j-1)*ntimes+k]] & snp$Trait[olist[[j]]@subjectHits]==traits[i] ]   ))
					}
				}
				params$crit=crit
			}
			b=enrichTest(params,flag)
			pvalues[i]=b$pvalue
			odds.ratios[i]=b$odds.ratio
		}
		
		#trait-class-specific testing
		for(i in seq_len(ntraitclass)){	
			if(i%%10==0) message(i," trait class have been tested!")
			hits.class[i]=sum(ttclass$"Trait_Class"==traitclass[i])
			snpnums.class[i]=sum(snpclass$"Trait_Class"==traitclass[i])
			params$nsnp=snpnums.class[i]
			params$hits=hits.class[i]
			if(test.method=="nonparametric"){
				for(j in seq_len(nbatch)){
					for(k in seq_len(ntimes)){
						crit.class[(j-1)*ntimes+k]=length(unique(  olist.class[[j]]@subjectHits[ booleans.class[[(j-1)*ntimes+k]] & snpclass$Trait_Class[olist.class[[j]]@subjectHits]==traitclass[i] ]   ))
					}
				}
				params$crit=crit.class
			}
			b=enrichTest(params,flag)
			pvalues.class[i]=b$pvalue
			odds.ratios.class[i]=b$odds.ratio
		}
		
		qvalues=p.adjust(pvalues,"fdr")
		tb1=data.frame(traits,pvalues,qvalues,odds.ratios,hits,snpnums)
		colnames(tb1)=c("Trait","p.value","q.value","odds.ratio","taSNP.hits","taSNP.num")
		rownames(tb1)=NULL
		
		qvalues.class=p.adjust(pvalues.class,"fdr")
		tb2=data.frame(traitclass,pvalues.class,qvalues.class,odds.ratios.class,hits.class,snpnums.class)
		colnames(tb2)=c("Trait_Class","p.value","q.value","odds.ratio","taSNP.hits","taSNP.num")
		rownames(tb2)=NULL
		
	
	}else{
	
		### test overall trait-associated SNP enrichment
		params$nsnp1=nsnp
		params$hits1=allhits
		params$nsnp2=nsnp.bg
		params$hits2=allhits.bg
		b=enrichTest(params,flag)
		tb.all=data.frame(Trait="All",p.value=b$pvalue,odds.ratio=b$odds.ratio,hits=allhits,taSNP.num=nsnp)
		colnames(tb.all)=c("Trait","p.value","odds.ratio","taSNP.hits","taSNP.num")
		rownames(tb.all)=NULL
		
		### trait-specific SNP enrichment test
		for(i in seq_len(ntraits)){
			if(i%%100==0) message(i," traits have been tested!")
			hits[i]=sum(tt$"Trait"==traits[i])
			snpnums[i]=sum(snp$"Trait"==traits[i])
			params$nsnp1=snpnums[i]
			params$hits1=hits[i]
			params$nsnp2=nsnp.bg
			params$hits2=allhits.bg
			if(test.method=="nonparametric"){
				for(j in seq_len(nbatch)){
					for(k in seq_len(ntimes)){
						crit[(j-1)*ntimes+k]=length(unique(  olist[[j]]@subjectHits[ booleans[[(j-1)*ntimes+k]] & snp$Trait[olist[[j]]@subjectHits]==traits[i] ]   ))
					}
				}
				params$crit=crit
			}
			b=enrichTest(params,flag)
			pvalues[i]=b$pvalue
			odds.ratios[i]=b$odds.ratio
		}
		
		
		#trait-class-specific testing
		for(i in seq_len(ntraitclass)){	
			if(i%%10==0) message(i," trait class have been tested!")
			hits.class[i]=sum(ttclass$"Trait_Class"==traitclass[i])
			snpnums.class[i]=sum(snpclass$"Trait_Class"==traitclass[i])
			params$nsnp1=snpnums.class[i]
			params$hits1=hits.class[i]
			params$nsnp2=nsnp.bg
			params$hits2=allhits.bg
			if(test.method=="nonparametric"){
				for(j in seq_len(nbatch)){
					for(k in seq_len(ntimes)){
						crit.class[(j-1)*ntimes+k]=length(unique(  olist.class[[j]]@subjectHits[ booleans.class[[(j-1)*ntimes+k]] & snpclass$Trait_Class[olist.class[[j]]@subjectHits]==traitclass[i] ]   ))
					}
				}
				params$crit=crit.class
			}
			b=enrichTest(params,flag)
			pvalues.class[i]=b$pvalue
			odds.ratios.class[i]=b$odds.ratio
		}
		
		qvalues=p.adjust(pvalues,"fdr")
		tb1=data.frame(traits,pvalues,qvalues,odds.ratios,hits,snpnums)
		colnames(tb1)=c("Trait","p.value","q.value","odds.ratio","taSNP.hits","taSNP.num")
		rownames(tb1)=NULL
		
		qvalues.class=p.adjust(pvalues.class,"fdr")
		tb2=data.frame(traitclass,pvalues.class,qvalues.class,odds.ratios.class,hits.class,snpnums.class)
		colnames(tb2)=c("Trait_Class","p.value","q.value","odds.ratio","taSNP.hits","taSNP.num")
		rownames(tb2)=NULL
				
	}
	
	### sort trait by either pvalue or odds ratio
	if(rankby=="pvalue"){
		tb1=tb1[order(tb1[,"p.value"],-tb1[,"odds.ratio"],decreasing=FALSE),]
		tb2=tb2[order(tb2[,"p.value"],-tb2[,"odds.ratio"],decreasing=FALSE),]
	}else if(rankby=="odds.ratio"){
		tb1=tb1[order(tb1[,"odds.ratio"],-tb1[,"p.value"],decreasing=TRUE),]
		tb2=tb2[order(tb2[,"p.value"],-tb2[,"odds.ratio"],decreasing=FALSE),]
	}
	x=list(tb.all=tb.all,tb1=tb1,tb2=tb2,ntraits=ntraits,ntraitclass=ntraitclass)
	structure(x,class="traseR")
}




#################################################################################
### print.traseR
### print the overall SNP enrichment 
### print the trait-specific SNP enrichment above bonferroni correction threshold


print.traseR<-function(x,isTopK1=FALSE,topK1=10,isTopK2=FALSE,topK2=10,
	trait.threshold=10,traitclass.threshold=10,...){
	message("There are ",x$ntraits," traits in the test.")
	message("The overall functional SNP enrichment test results are:")
	print(x$tb.all)
	message("The trait-associated SNP enrichment test results are:")
	if(isTopK1){
		print(x$tb1[seq_len(topK1),])
	}else{
		print(x$tb1[x$tb1$p.value<0.05/x$ntraits & x$tb1$taSNP.num>trait.threshold,])
	}
	message("The trait-class-associated SNP enrichment test results are:")
	if(isTopK2){
		print(x$tb2[seq_len(topK2),])
	}else{
		print(x$tb2[x$tb2$p.value<0.05/x$ntraitclass & x$tb2$taSNP.num>traitclass.threshold,])
	}
}


