\name{taSNP}
\alias{taSNP}
\title{
trait-associated SNPs in dbGaP and NHGRI downloaded from Association Results Browser
}
\description{
A GRange object \code{taSNP} contains trait-associated SNPs from dbGaP and NHGRI downloaded from Association Results Browser.
}

\usage{
data(taSNP)
}

\value{
	The data frame \code{taSNP} contains the following columns
	\item{Trait}{Trait}
	\item{Trait Class}{Trait class which is formed based on the phenotype tree. Close traits are grouped together to form one class}
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



