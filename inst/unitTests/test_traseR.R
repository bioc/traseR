data(Tcell)
data(taSNP)

test.traseR<-function(){
    x=traseR(snpdb=taSNP,region=Tcell,rankby="pvalue",test.method="binomial",alternative ="greater")
    checkIdentical("traseR",as.character(class(x)))
    checkIdentical("data.frame",as.character(class(x$tb.all)))
    checkIdentical("data.frame",as.character(class(x$tb1)))
    checkIdentical("data.frame",as.character(class(x$tb2)))
    checkIdentical("integer",as.character(class(x$ntraits)))
    checkIdentical("integer",as.character(class(x$ntraitclass)))
    checkEquals(length(x),5)
    checkEquals(colnames(x$tb.all),c("Trait" ,"p.value","odds.ratio" ,"taSNP.hits" , "taSNP.num") )
    checkEquals(colnames(x$tb1),c("Trait" ,"p.value" ,"q.value","odds.ratio" ,"taSNP.hits" ,"taSNP.num") )
    checkEquals(colnames(x$tb2),c("Trait_Class" ,"p.value" ,"q.value","odds.ratio" ,"taSNP.hits" ,"taSNP.num") )
    checkEquals(ncol(x$tb.all),5)
	checkEquals(ncol(x$tb1),6)
	checkEquals(ncol(x$tb2),6)
	checkException(traseR(snpdb=taSNP,region=data.frame(chr="chr1",start=1,end=100),rankby="pvalue",test.method="binomial",alternative ="greater"))
 }
 
 

 test.queryKeyword<-function(){
 	 x=queryKeyword(snpdb=taSNP,region=Tcell,keyword="autoimmune",returnby="SNP_ID")
     checkIdentical("data.frame",as.character(class(x)))
     checkEquals(ncol(x),5)
     checkEquals(names(x), c("SNP_ID","Chr","Position","Trait.num","Trait.name"))
     checkException(queryKeyword(snpdb=taSNP,region=data.frame(chr="chr1",start=1,end=100),keyword="autoimmune",returnby="SNP_ID"))
     x=queryKeyword(snpdb=taSNP,region=Tcell,keyword="autoimmune",returnby="trait")
     checkIdentical("data.frame",as.character(class(x)))
     checkEquals(ncol(x),3)
     checkEquals(names(x), c("Trait","taSNP.num","taSNP.name"))
     checkException(queryKeyword(snpdb=taSNP,region=data.frame(chr="chr1",start=1,end=100),keyword="autoimmune",returnby="trait"))
 }
 
 
 
 test.querySNP<-function(){
 	 x=querySNP(snpdb=taSNP,snpid=c("rs3766178","rs880051"))
     checkIdentical("GRanges",as.character(class(x)))
     checkEquals(names(mcols(x)),c("Trait","SNP_ID","p.value","Context","GENE_NAME","GENE_START","GENE_END","GENE_STRAND","Trait_Class"))
     checkException(querySNP(snpdb=taSNP,snpid=c("rs1","rs2")))
 }
 
 
 
 test.queryGene<-function(){
 	  x=queryGene(snpdb=taSNP,genes=c("AGRN","UBE2J2","SSU72"))
      checkIdentical("GRanges",as.character(class(x)))
      checkEquals(names(mcols(x)),c("GENE_NAME","Trait.num","Trait.name","taSNP.num","taSNP.name"))
      checkException(queryGene(snpdb=taSNP,genes=c("gene1","gene2")))
 }
 
 

 
