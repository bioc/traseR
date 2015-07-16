data(Tcell)
data(ASB)

test.traseR<-function(){
    x=traseR(snpdb=ASB,region=Tcell,rankby="pvalue",test.method="binomial",alternative ="greater")
    checkIdentical("traseR",as.character(class(x)))
    checkIdentical("data.frame",as.character(class(x$tb.all)))
    checkIdentical("data.frame",as.character(class(x$tb)))
    checkIdentical("integer",as.character(class(x$ntraits)))
    checkEquals(length(x),3)
    checkEquals(colnames(x$tb.all),c("Trait" ,"p.value","odds.ratio" ,"hits" , "SNP_Num") )
    	checkEquals(colnames(x$tb),c("Trait" ,"p.value" ,"q.value","odds.ratio" ,"hits" ,"SNP_Num") )
    	checkEquals(ncol(x$tb.all),5)
	checkEquals(ncol(x$tb),6)
	checkException(traseR(snpdb=ASB,region=data.frame(chr="chr1",start=1,end=100),rankby="pvalue",test.method="binomial",alternative ="greater"))
 }
 


 test.queryKeyword<-function(){
 	 x=queryKeyword(snpdb=ASB,region=Tcell,keyword="autoimmune",returnby="SNP")
     checkIdentical("data.frame",as.character(class(x)))
     checkEquals(ncol(x),5)
     checkEquals(names(x), c("SNP","Chr","Position","Trait_Num","Trait_Name"))
     checkException(queryKeyword(snpdb=ASB,region=data.frame(chr="chr1",start=1,end=100),keyword="autoimmune",returnby="SNP"))
     x=queryKeyword(snpdb=ASB,region=Tcell,keyword="autoimmune",returnby="trait")
     checkIdentical("data.frame",as.character(class(x)))
     checkEquals(ncol(x),3)
     checkEquals(names(x), c("Trait","SNP_Num","SNP_Name"))
     checkException(queryKeyword(snpdb=ASB,region=data.frame(chr="chr1",start=1,end=100),keyword="autoimmune",returnby="trait"))
 }
 
 
 
 test.querySNP<-function(){
 	 x=querySNP(snpdb=ASB,snpid=c("rs3766178","rs880051"))
     checkIdentical("GRanges",as.character(class(x)))
     checkEquals(names(mcols(x)),c("Trait","SNP","p.value","Context","GENE_NAME","GENE_START","GENE_END","GENE_STRAND"))
     checkException(querySNP(snpdb=ASB,snpid=c("rs1","rs2")))
 }
 
 
 
 test.queryGene<-function(){
 	  x=queryGene(snpdb=ASB,genes=c("AGRN","UBE2J2","SSU72"))
      checkIdentical("GRanges",as.character(class(x)))
      checkEquals(names(mcols(x)),c("GENE_NAME","Trait_Num","Trait_Name","SNP_Num","SNP_Name"))
      checkException(queryGene(snpdb=ASB,genes=c("gene1","gene2")))
 }
 
 
