\name{taSNPLD}
\alias{taSNPLD}
\title{
linkage disequilibrium (>0.8) within 100kb SNPs of all trait-associated SNPs from dbGaP and NHGRI
}
\description{
A GRange object \code{taSNPLD} contains linkage disequilibrium (>0.8) SNPs of all trait-associated SNPs from dbGaP and NHGRI.}

\usage{
data(taSNPLD)
}

\value{	
	The data frame \code{taSNPLD} contains four columns,
	\item{SNP_ID}{SNP rs number}
	\item{seqnames}{Chromosome number associated with rs number}
	\item{ranges}{Chromosomal position, in base pairs, associated with rs number}
	\item{Trait}{Trait the SNP is associated with}
	\item{Trait Class}{Trait class which is formed based on the phenotype tree. Close traits are grouped together to form one class}
}

\author{
Li Chen <li.chen@emory.edu>, Zhaohui S.Qin<zhaohui.qin@emory.edu>
}


