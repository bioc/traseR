\name{print.traseR}
\alias{print.traseR}
\title{
	Print the outcome of taSNP enrichment analyses
}
\description{
	Print the outcome of taSNP enrichment analyses. Print the overall taSNP enrichment and trait-specific taSNP enrichment.
}

\usage{
	\method{print}{traseR}(x,isTopK=FALSE,topK=10,trait.threshold=10,...)
}

\arguments{
	\item{x}{Object returned from \code{traseR}}
	\item{isTopK}{If isTopK is TRUE, topK traits are printed; otherwise, traits with p-value below Bonferroni correction threshold are printed. Default is FALSE.}
	\item{topK}{Top K traits are printed. Default is 10.}
	\item{trait.threshold}{Traits above this threshold are reported. Default is 10.}
  	\item{\dots}{Other parameters to \code{print}}
}

\examples{
	data(taSNPDB)
	data(Tcell)
	x=traseR(snpdb=taSNPDB,region=Tcell)
	print(x)
}

\value{
	Print a data frame of traits ranked by p-value
}

\author{
Li Chen <li.chen@emory.edu>, Zhaohui S.Qin<zhaohui.qin@emory.edu>
}

