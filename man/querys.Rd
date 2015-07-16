\name{querys}
\alias{queryKeyword}
\alias{queryGene}
\alias{querySNP}


\title{
  Retrieve trait-associated SNPs based
}
\description{
  These are a group of functions to retrieve the
  trait-associated SNPs based on input
}

\usage{
	
queryKeyword(snpdb, region=NULL, keyword = NULL, returnby = c("SNP", "trait"), pvalue = 1e-3)

queryGene(snpdb, genes = NULL)

querySNP(snpdb, snpid, region = NULL)

}


\arguments{
  \item{snpdb}{ A GRange object or data frame, which is GWAS trait-associated SNPs downloaded from up-to-date dbGaP and NHGRI public database. It is maintained to be updated to the latest version. The data frame contains the following columns,\code{Trait},\code{SNP},\code{p.value},\code{Chr},\code{Position},\code{Context},\code{GENE_NAME},\code{GENE_START},\code{GENE_END},\code{GENE_STRAND}. The data frame is in \code{data} subdirectory. Users are free to add more SNP records to the data frame for practical use.
  	}
 \item{region}{
A data frame, which is genomic intervals with three columns, chromosome, genomic start position, genomic end position.
}
  \item{keyword}{
The keyword is used when specific trait is of interest. If \code{keyword} is specified, only the SNPs associated to the trait are used for analyses. Otherwise, all traits will be analyzed.
}
  \item{snpid}{SNP name}
  \item{genes}{Gene names}
  \item{pvalue}{
SNPs with p-value less than this threshold are used for analyses.
}
   \item{returnby}{
  	Either SNP or trait. If \code{returnby} is specified as 'SNP', a data frame based on 'SNP' is returned. If \code{returnby} is specified as 'trait', a data frame based on 'trait' is returned.
}

}

\author{
  Li Chen <li.chen@emory.edu>, Zhaohui Qin<zhaohui.qin@emory.edu>
}

\value{
	queryKeyword: Return a data frame of traits containing the keyword
	queryGene: Return a data frame of traits associated with the gene
	querySNP: Return a data frame of traits associated with the SNP
}

\examples{
	
	data(ASB)
	data(Tcell)
	x=queryKeyword(snpdb=ASB,region=Tcell,keyword="Autoimmune",returnby="SNP")
	x=queryGene(snpdb=ASB,genes=c("AGRN","UBE2J2","SSU72"))
	x=querySNP(snpdb=ASB,snpid=c("rs3766178","rs880051"))
	
}



