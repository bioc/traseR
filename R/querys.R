
##################################################
### Return SNPs given query key words

queryKeyword<-function(snpdb,region=NULL,keyword=NULL,returnby=c("SNP_ID","trait"),pvalue=1e-3){

	returnby=match.arg(returnby)
	snpdb=snpdb[snpdb$p.value<=pvalue]
	snp=unique(as.data.frame(snpdb[,c("SNP_ID","p.value","Trait")]))
	snp.meta=snp[,setdiff(colnames(snp),c("seqnames","start","end","strand","width","Trait_Class"))]
	snp=makeGRangesFromDataFrame(snp[,c("seqnames","start","end","strand")])
	mcols(snp)=snp.meta
		
	if(!is.null(region)){
		if(class(region)!="data.frame" & class(region)!="GRanges")
			stop("Query genomic interval should be either data.frame or GRanges!")
		if(class(region)=="data.frame"){
			if(!identical(colnames(region),c("chr","start","end")))
				stop("Interval colnames should be chr,start,end!")
			region=makeGRangesFromDataFrame(region)
		}
		seqlevel=c(paste("chr",1:22,sep=""),"chrX")
		region=region[!is.na(match(seqnames(region),seqlevel))]
		o=findOverlaps(region,snp)
		if(length(subjectHits(o))==0){
			stop("No SNP overlapped with the query region!")
		}else{
			snp=unique(as.data.frame(snp[unique(subjectHits(o)),c("SNP_ID","p.value","Trait")]))
		}
	}
	
	if(returnby=="trait"){
		if(!is.null(keyword)){
			ind=grep(tolower(keyword),tolower(snp$Trait))
			if(length(ind)==0){
				stop("No matched key word found in Trait!")
			}else{								
				snp=unique(snp[ind,c("SNP_ID","seqnames","start","Trait")])		
			}
		}
		a=by(snp[,c("SNP_ID","Trait")], snp$Trait, function(x) x[1])
		trait.name=names(a)
		snp.name=sapply(a,function(x) paste(x,collapse=";"))
		snp.num=sapply(a,length)
		ind=match(trait.name,snp$Trait)
		names(snp.num)=names(snp.name)=NULL
		tb=data.frame(Trait=snp[ind,"Trait"],taSNP.num=snp.num,taSNP.name=snp.name)
		colnames(tb)=c("Trait","taSNP.num","taSNP.name")
		rownames(tb)=NULL	
		
	}else if(returnby=="SNP_ID"){
		if(!is.null(keyword)){
			ind=grep(tolower(keyword),tolower(snp$Trait))
			if(length(ind)==0){
				stop("No matched key word found in Trait!")
			}else{
				snp=unique(snp[ind,c("SNP_ID","seqnames","start","Trait")])
			}
		}
		a=by(snp[,c("SNP_ID","Trait")], snp$SNP_ID, function(x) x[2])
		snp.name=names(a)
		trait.name=sapply(a,function(x) paste(x,collapse=";"))
		trait.num=sapply(a,length)
		ind=match(snp.name,snp$SNP_ID)
		tb=cbind(snp[ind,c("SNP_ID","seqnames","start")],trait.num,trait.name)
		colnames(tb)=c("SNP_ID","Chr","Position","Trait.num","Trait.name")
	}
	tb
}




##################################################
### Return nearby SNPs given query gene

queryGene<-function(snpdb,genes=NULL){
	snpdb=as.data.frame(snpdb)
	if(is.null(genes)){
		snp=unique(snpdb[,c("seqnames","GENE_NAME","GENE_START","GENE_END","GENE_STRAND","SNP_ID","Trait")])
	}else{
		ind=match(tolower(genes),tolower(snpdb$GENE_NAME))
		if(length(ind)==0){
			stop("No matched gene found!")
		}else{
			snp=unique(snpdb[ind,c("seqnames","GENE_NAME","GENE_START","GENE_END","GENE_STRAND","SNP_ID","Trait")])
		}
	}
	snp1=unique(snp[,c("GENE_NAME","SNP_ID")])
	snp2=unique(snp[,c("GENE_NAME","Trait")])
	a=by(snp1,snp1$GENE_NAME, function(x) x[2])
	gene.name=names(a)
	snp.name=sapply(a,function(x) paste(x,collapse=";"))
	snp.num=sapply(a,length)
	a=by(snp2,snp2$GENE_NAME, function(x) x[2])
	trait.name=sapply(a,function(x) paste(x,collapse=";"))
	trait.num=sapply(a,length)
	gene.loci=snp[match(gene.name,snp$GENE_NAME),c("seqnames","GENE_START","GENE_END","GENE_STRAND")]
	tb=cbind(gene.name,gene.loci,trait.num,trait.name,snp.num,snp.name)
	rownames(tb)=NULL
	colnames(tb)=c("GENE_NAME","Chr","GENE_START","GENE_END","GENE_STRAND","Trait.num","Trait.name","taSNP.num","taSNP.name")
	tb=makeGRangesFromDataFrame(tb,seqnames.field="Chr",start.field="GENE_START",end.field="GENE_END",
	strand.field="GENE_STRAND",keep.extra.columns = TRUE)
	tb
}



##################################################
### Return nearby gene given query SNP

querySNP<-function(snpdb,snpid,region=NULL){
	if(!is.null(snpid)){
		ind=match(tolower(snpid),tolower(snpdb$SNP_ID))
		if(length(ind)==0){
			stop("No matched SNP found!")
		}else{
			snp=snpdb[ind]
		}
	}
	if(!is.null(region)){
		if(class(region)!="data.frame" & class(region)!="GRanges")
			stop("Query genomic interval should be either data.frame or GRanges!")
		if(class(region)=="data.frame"){
			if(!identical(colnames(region),c("chr","start","end")))
				stop("Interval colnames should be chr,start,end!")
			region=makeGRangesFromDataFrame(region)
		}
		seqlevel=c(paste("chr",1:22,sep=""),"chrX")
		region=region[!is.na(match(seqnames(region),seqlevel))]
		o=findOverlaps(region,snp)
		if(length(subjectHits(o))==0){
			stop("No SNP overlapped with the query region!")
		}else{
			snp=snp[unique(subjectHits(o))]
		}
	}
	snp
}




















