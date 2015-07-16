################################################################
#### Core function to perform hypothesis testing.
###  It is a nested function of enrichTest, no export for users.

enrichTest<-function(params,flag){

		if(flag){
			tb=matrix(c(params$hits,params$nsnp-params$hits,params$region.size-params$hits,
				params$genome.size-params$region.size-params$nsnp+params$hits),2,2)
			rownames(tb)=c("region.in","region.out")
			colnames(tb)=c("snp.yes","snp.no")
		}else{
			tb=matrix(c(params$hits1,params$hits2,params$nsnp1-params$hits1,params$nsnp2-params$hits2),2,2)
			rownames(tb)=c("region.in","region.out")
			colnames(tb)=c("snp.trait","snp.background")
		}
		
		odds.ratio=((1+tb[1,1])/(1+tb[1,ncol(tb)]))/((1+tb[nrow(tb),1])/(1+tb[nrow(tb),ncol(tb)]))	

		if(params$test.method=="chisq"){
			a=chisq.test(tb)
			pvalue=a$p.value
		}else if(params$test.method=="binomial"){
			if(flag){
				a=binom.test(params$hits,params$region.size,
					p=params$nsnp/params$genome.size,alternative=params$alternative)
			}else{
				a=binom.test(params$hits1,params$nsnp1,
					p=(params$hits1+params$hits2)/(params$nsnp1+params$nsnp2),alternative=params$alternative)
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
	test.method=c("chisq","binomial","fisher"),
	alternative = c("greater","less","two.sided"),
	trait.threshold=10,pvalue=1e-3){
	
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
	region.size=sum(end(region)-start(region))

	### region size and hg19 genome size
	seqlength=seqlengths(BSgenome.Hsapiens.UCSC.hg19)
	seqlength=seqlength[match(seqlevel,names(seqlength))]
	genome.size=sum(as.numeric(seqlength))
	
	
	### load SNP database in NHGRI format,calculate overlap with query region
	snp=unique(as.data.frame(snpdb[,c("SNP","Trait")]))
	snp.meta=snp[,setdiff(colnames(snp),c("seqnames","start","end","strand","width"))]
	snp=makeGRangesFromDataFrame(snp[,c("seqnames","start","end","strand")])
	mcols(snp)=snp.meta
	
	if(!is.null(keyword)){
		ind=grep(tolower(keyword),tolower(snp$Trait))
		if(length(ind)==0)
			stop("No match key word found in Trait!")
		snp=snp[ind]
	}
	
	snp=snp[order(snp$Trait)]
	traits=unique(snp$Trait)
	traits=unique(names(table(snp$Trait))[table(snp$Trait)>trait.threshold])
	snp=snp[!is.na(match(snp$Trait,traits))]
	end(snp)=end(snp)+1
	
	o=findOverlaps(region,snp)	
	allhits=length(unique(o@subjectHits))
	if(allhits==0){
		stop("Overall SNPs are not overlapped with any genomic intervals!")
	}
	nsnp=length(unique(snp$SNP))
	tt=snp[unique(o@subjectHits)]
	ntraits=length(traits)
	
	### print out the size of region and number of traits associated with the region
	message(paste("There are ",region.size,"bp in the query region, "),paste("accounting for ",region.size/genome.size," of the genome."))
	message(paste("There are ",ntraits,"traits in the analysis."))

	
	### Calculate background SNP overlap
	if(!flag){
 		ind=match(snpdb.bg$SNP,snp$SNP)
		snp.bg=unique(snpdb.bg[is.na(ind)])
		end(snp.bg)=end(snp.bg)+1
		o=findOverlaps(region,snp.bg)	
		allhits.bg=length(unique(o@subjectHits))
		nsnp.bg=length(unique(snp.bg$SNP))
	}


	### parameters pass to enrichTest
	params=list()	
	params$region.size=region.size
	params$genome.size=genome.size
	params$test.method=test.method
	params$alternative=alternative
	### performance criteria
	pvalues=odds.ratios=hits=snpnums=rep(-1,ntraits)
	

	if(flag){
		### test overall trait-associated SNP enrichment
		params$nsnp=nsnp
		params$hits=allhits
		b=enrichTest(params,flag)
		tb.all=data.frame(Trait="All",p.value=b$pvalue,odds.ratio=b$odds.ratio,hits=allhits,SNP_Num=nsnp)
		rownames(tb.all)=NULL
		for(i in seq_len(ntraits)){
			hits[i]=sum(tt$"Trait"==traits[i])
			snpnums[i]=sum(snp$"Trait"==traits[i])
			params$nsnp=snpnums[i]
			params$hits=hits[i]
			b=enrichTest(params,flag)
			pvalues[i]=b$pvalue
			odds.ratios[i]=b$odds.ratio
		}
		qvalues=p.adjust(pvalues,"fdr")
		tb=data.frame(traits,pvalues,qvalues,odds.ratios,hits,snpnums)
		colnames(tb)=c("Trait","p.value","q.value","odds.ratio","hits","SNP_Num")
		rownames(tb)=NULL
	
	}else{
	
		### test overall trait-associated SNP enrichment
		params$nsnp1=nsnp
		params$hits1=allhits
		params$nsnp2=nsnp.bg
		params$hits2=allhits.bg
		b=enrichTest(params,flag)
		tb.all=data.frame(Trait="All",p.value=b$pvalue,odds.ratio=b$odds.ratio,hits=allhits,SNP_Num=nsnp)
		rownames(tb.all)=NULL
		
		### trait-specific SNP enrichment test
		for(i in seq_len(ntraits)){
			hits[i]=sum(tt$"Trait"==traits[i])
			snpnums[i]=sum(snp$"Trait"==traits[i])
			params$nsnp1=snpnums[i]
			params$hits1=hits[i]
			params$nsnp2=nsnp.bg
			params$hits2=allhits.bg
			b=enrichTest(params,flag)
			pvalues[i]=b$pvalue
			odds.ratios[i]=b$odds.ratio
		}
		qvalues=p.adjust(pvalues,"fdr")
		tb=data.frame(traits,pvalues,qvalues,odds.ratios,hits,snpnums)
		colnames(tb)=c("Trait","p.value","q.value","odds.ratio","hits","SNP_Num")
		rownames(tb)=NULL
	}
	
	### sort trait by either pvalue or odds ratio
	if(rankby=="pvalue"){
		tb=tb[order(tb[,"p.value"],-tb[,"odds.ratio"],decreasing=FALSE),]
	}else if(rankby=="odds.ratio"){
		tb=tb[order(tb[,"odds.ratio"],-tb[,"p.value"],decreasing=TRUE),]
	}
	x=list(tb.all=tb.all,tb=tb,ntraits=ntraits)
	structure(x,class="traseR")
}


#################################################################################
### print.traseR
### print the overall SNP enrichment 
### print the trait-specific SNP enrichment above bonferroni correction threshold

print.traseR<-function(x,topK=5,...){
	message("There are ",x$ntraits," traits in the test.")
	message("The overall functional SNP enrichment test results are:")
	print(x$tb.all)
	if(sum(x$tb$p.value<0.05)>topK){
		message(paste("TopK",topK,"trait-associated SNP enrichment test results are:"))
		print(x$tb[x$tb$p.value<0.05/x$ntraits,][seq_len(topK),])
	}else{
		message("The trait-associated SNP enrichment test results are:")
		print(x$tb[x$tb$p.value<0.05/x$ntraits,])
	}
}





