\name{plots}
\alias{plotContext}
\alias{plotGene}
\alias{plotSNP}
\alias{plotPvalue}
\alias{plotInterval}

\title{
  Visualize of trait-associated SNPs
}
\description{
  These are a group of functions to generate plot to visualize the
  trait-associated SNPs. 
}

\usage{
	
plotContext(snpdb, region=NULL, keyword = NULL, pvalue = 1e-3)

plotPvalue(snpdb, region=NULL, keyword = NULL, plot.type = c("densityplot", "boxplot"), pvalue = 1e-3, xymax = 50)

plotSNP(snpdb, snpid, ext = 10000)

plotGene(snpdb, gene, ext = 10000)

plotInterval(snpdb,interval,ext = 10000)

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
  \item{gene}{Gene name}
  \item{pvalue}{
SNPs with p-value less than this threshold are used for analyses.
}
  \item{plot.type}{Either "densityplot" or "boxplot"}
  \item{ext}{Bp extended upstream and downstream}
  \item{xymax}{The maximum range on x-axis and y-axis}
  \item{interval}{A data frame, genomic interval:chromosome, genomic start position, genomic end position}

}

\value{
	\item{plotContext}{A pie plot with the distribution of SNP function class}
	\item{plotPvalue}{A density plot of -logPvalue of trait-associated SNPs}
	\item{plotSNP}{A plot of trait-associated SNP on chromosome}
	\item{plotGene}{A plot with the gene and possible nearby trait-associated SNPs }
	\item{plotInterval}{A plot with chromosome interval with possible nearby genes and trait-associated SNPs }
}

\author{
  Li Chen <li.chen@emory.edu>, Zhaohui Qin<zhaohui.qin@emory.edu>
}


\examples{
	
	data(ASB)
	plotContext(snpdb=ASB,keyword="Autoimmune")
	plotGene(snpdb=ASB,gene="ZFP92",ext=50000)
	plotSNP(snpdb=ASB,snpid="rs766420",ext=50000)
	plotInterval(snpdb=ASB,data.frame(chr="chrX",start=152633780,end=152737085))
}






