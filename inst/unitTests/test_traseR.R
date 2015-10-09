data(Tcell)
data(taSNPDB)

test.traseR<-function(){
    x=traseR(snpdb=taSNPDB,region=Tcell,rankby="pvalue",test.method="binomial",alternative ="greater")
    checkIdentical("traseR",as.character(class(x)))
    checkIdentical("data.frame",as.character(class(x$tb.all)))
    checkIdentical("data.frame",as.character(class(x$tb)))
    checkIdentical("integer",as.character(class(x$ntraits)))
    checkEquals(length(x),3)
    checkEquals(colnames(x$tb.all),c("Trait" ,"p.value","odds.ratio" ,"taSNP.hits" , "taSNP.num") )
    	checkEquals(colnames(x$tb),c("Trait" ,"p.value" ,"q.value","odds.ratio" ,"taSNP.hits" ,"taSNP.num") )
    	checkEquals(ncol(x$tb.all),5)
	checkEquals(ncol(x$tb),6)
	checkException(traseR(snpdb=taSNPDB,region=data.frame(chr="chr1",start=1,end=100),rankby="pvalue",test.method="binomial",alternative ="greater"))
 }
 


 test.queryKeyword<-function(){
 	 x=queryKeyword(snpdb=taSNPDB,region=Tcell,keyword="autoimmune",returnby="SNP_ID")
     checkIdentical("data.frame",as.character(class(x)))
     checkEquals(ncol(x),5)
     checkEquals(names(x), c("SNP_ID","Chr","Position","Trait.num","Trait.name"))
     checkException(queryKeyword(snpdb=taSNPDB,region=data.frame(chr="chr1",start=1,end=100),keyword="autoimmune",returnby="SNP_ID"))
     x=queryKeyword(snpdb=taSNPDB,region=Tcell,keyword="autoimmune",returnby="trait")
     checkIdentical("data.frame",as.character(class(x)))
     checkEquals(ncol(x),3)
     checkEquals(names(x), c("Trait","taSNP.num","taSNP.name"))
     checkException(queryKeyword(snpdb=taSNPDB,region=data.frame(chr="chr1",start=1,end=100),keyword="autoimmune",returnby="trait"))
 }
 
 
 
 test.querySNP<-function(){
 	 x=querySNP(snpdb=taSNPDB,snpid=c("rs3766178","rs880051"))
     checkIdentical("GRanges",as.character(class(x)))
     checkEquals(names(mcols(x)),c("Trait","SNP_ID","p.value","Context","GENE_NAME","GENE_START","GENE_END","GENE_STRAND"))
     checkException(querySNP(snpdb=taSNPDB,snpid=c("rs1","rs2")))
 }
 
 
 
 test.queryGene<-function(){
 	  x=queryGene(snpdb=taSNPDB,genes=c("AGRN","UBE2J2","SSU72"))
      checkIdentical("GRanges",as.character(class(x)))
      checkEquals(names(mcols(x)),c("GENE_NAME","Trait.num","Trait.name","taSNP.num","taSNP.name"))
      checkException(queryGene(snpdb=taSNPDB,genes=c("gene1","gene2")))
 }
 
 

 
