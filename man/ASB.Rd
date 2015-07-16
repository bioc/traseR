\name{ASB}
\alias{ASB}
\title{
trait-associated SNPs in dbGaP and NHGRI downloaded from Association Results Browser
}
\description{
A GRange object \code{ASB} contains trait-associated SNPs from dbGaP and NHGRI downloaded from Association Results Browser.
}

\usage{
data(ASB)
}

\value{
	The data frame \code{ASB} contains the following columns
	\item{Trait}{Trait}
	\item{SNP}{SNP}
	\item{p.value}{GWAS SNP p-value}
	\item{seqnames}{Chromosome}
	\item{ranges}{Chromosome position}
	\item{Context}{SNP functional class}
	\item{GENE_NAME}{Gene name}
	\item{GENE_START}{Gene start genomic position}
	\item{GENE_END}{Gene end genomic position}
	\item{GENE_STRAND}{Gene strand}
}


\author{
Li Chen <li.chen@emory.edu>, Zhaohui S.Qin<zhaohui.qin@emory.edu>
}



