\name{print.traseR}
\alias{print.traseR}
\title{
	Print the outcome of taSNP enrichment analyses
}
\description{
	Print the outcome of taSNP enrichment analyses. Print the overall taSNP enrichment, trait-specific taSNP enrichment,trait-class-specific taSNP enrichment.
}

\usage{
	\method{print}{traseR}(x,isTopK1=FALSE,topK1=10,isTopK2=FALSE,topK2=10,trait.threshold=10,traitclass.threshold=10,...)
}

\arguments{
	\item{x}{Object returned from \code{traseR}}
	\item{isTopK1}{If isTopK1 is TRUE, topK1 traits are printed; otherwise, traits with p-value below Bonferroni correction threshold are printed. Default is FALSE.}
	\item{topK1}{Top K1 traits are printed. Default is 10.}
	\item{isTopK2}{If isTopK2 is TRUE, topK2 trait class are printed; otherwise, trait class with p-value below Bonferroni correction threshold are printed. Default is FALSE.}
	\item{topK2}{Top K2 trait class are printed. Default is 10.}
	\item{trait.threshold}{Traits above this threshold are reported. Default is 10.}
	\item{traitclass.threshold}{Trait class above this threshold are reported. Default is 10.}
  	\item{\dots}{Other parameters to \code{print}}
}

\examples{
	data(taSNP)
	data(Tcell)
	x=traseR(snpdb=taSNP,region=Tcell)
	print(x)
}

\value{
	Print a data frame of traits ranked by p-value
}

\author{
Li Chen <li.chen@emory.edu>, Zhaohui S.Qin<zhaohui.qin@emory.edu>
}




