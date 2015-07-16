\name{print.traseR}
\alias{print.traseR}
\title{
	Print the outcome of SNP enrichment analyses
}
\description{
	Print the outcome of SNP enrichment analyses. Print the overall SNP enrichment and trait-specific SNP enrichment.
	Only traits above Bonferroni correction threshold are printed.
}

\usage{
	\method{print}{traseR}(x,topK=5,...)
}

\arguments{
	\item{x}{Object returned from \code{traseR}}
	\item{topK}{Top number of traits are printed. Default is 5.}
  	\item{\dots}{Other parameters to \code{print}}
}

\examples{
	data(ASB)
	data(Tcell)
	x=traseR(snpdb=ASB,region=Tcell,rankby="pvalue",test.method="binomial",alternative ="greater")
	print(x)
}

\value{
	Print a data frame of traits ranked by pvalue
}

\author{
Li Chen <li.chen@emory.edu>, Zhaohui S.Qin<zhaohui.qin@emory.edu>
}

