\name{taSNPDB}
\alias{taSNPDB}
\title{
trait-associated SNPs in dbGaP and NHGRI downloaded from Association Results Browser
}
\description{
A GRange object \code{taSNPDB} contains trait-associated SNPs from dbGaP and NHGRI downloaded from Association Results Browser.
}

\usage{
data(taSNPDB)
}

\value{
	The data frame \code{taSNPDB} contains the following columns
	\item{Trait}{Trait}
	\item{SNP_ID}{SNP rs number}
	\item{p.value}{GWAS SNP p-value}
	\item{seqnames}{Chromosome}
	\item{ranges}{Chromosome position}
	\item{Context}{SNP functional class}
	\item{GENE_NAME}{Nearest gene name}
	\item{GENE_START}{Gene start genomic position}
	\item{GENE_END}{Gene end genomic position}
	\item{GENE_STRAND}{Gene strand}
}


\author{
Li Chen <li.chen@emory.edu>, Zhaohui S.Qin<zhaohui.qin@emory.edu>
}



