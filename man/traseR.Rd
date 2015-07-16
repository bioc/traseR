\name{traseR}
\alias{traseR}
\title{
TRait-Associated SNP EnRichment analyses
}
\description{
Perform GWAS trait-associated SNP enrichment analyses in genomic intervals using different approaches
}
\usage{
traseR(snpdb, region, snpdb.bg=NULL, keyword = NULL, rankby = c("pvalue", "odds.ratio"), 
test.method = c("chisq", "binomial", "fisher"), alternative = c("greater", "less", "two.sided"), trait.threshold = 10, pvalue = 1e-3)
}

\arguments{
  \item{snpdb}{ A GRange object.
 It could be GWAS trait-associated SNPs downloaded from up-to-date dbGaP and NHGRI public database. It is maintained to be updated to the latest version. The data frame contains the following columns,\code{Source},\code{Trait},\code{SNP},\code{p.value},\code{Chr},\code{Position},\code{Context},\code{GENE_NAME},\code{GENE_START},\code{GENE_END},\code{GENE_STRAND}. The data frame is in \code{data} subdirectory. Users are free to add more SNP records to the data frame for practical use. It could also be a data frame with columns as, \code{SNP},\code{Chr},\code{Position}.
}
  \item{region}{
A GRange object or data frame, which is genomic intervals with three columns, chromosome, genomic start position, genomic end position.
}
  \item{snpdb.bg}{
  	 A GRange object contains non-trait-associated SNPs. They are treated as background for statistical testing instead of whole genome as background if specified.
}
  \item{keyword}{
The keyword is used when specific trait is of interest. If \code{keyword} is specified, only the SNPs associated to the trait are used for analyses. Otherwise, all traits will be analyzed.
}
  \item{rankby}{
Traits could be ranked by either p-value or adds.ratio based on the enrichment level of trait-associated SNPs in genomic intervals.
}
  \item{test.method}{
Several hypothesis testing options are provided: \code{chisq}(Chi-squared test),\code{binomial}(binomial test),\code{fisher}(Fisher's exact test)
}
  \item{alternative}{
Indicate the alternative hypothesis. If \code{greater}, test if the genomic intervals are enriched in trait-associated SNPs than background. If \code{less}, test if the genomic intervals are depleted in trait-associated SNPs than background. If \code{two.sided}, test if there is difference between the enrichment of trait-associated SNPs in genomic intervals and in background.
}
  \item{trait.threshold}{
Test traits with number of SNPs more than the threshold.
}
  \item{pvalue}{
SNPs with p-value less than this threshold are used for analyses.
}
}
\details{
Return a list that contains two data frames and one integer. One data frame \code{tb.all} contains the results of 
enrichment analyses for all trait-associated SNPs in genomic intervals.
Another data frame \code{tb} contains the results of enrichment analyses for each trait-associated SNPs
in gneomic intervals separately. One integer indicates how many traits are analyzed.
}
\value{
The data frame \code{tb} has columns,
  \item{Trait }{Name of trait}
  \item{p.value }{P-value calculated from hypothesis testing}
  \item{q.value}{Adjusted p-value from multiple testing using FDR correction}
  \item{odds.ratio}{Odds ratio calculated based on number of trait-associated SNPs in genomic intervals,
  	number of trait-associated SNPs across whole genome, genomic intervals size (bps) and genome size (bps)}
  \item{hits}{Number of trait-associated SNPs in genomic intervals}
  \item{SNP_Num}{Number of SNPs for specific trait}
}

\author{
Li Chen <li.chen@emory.edu>, Zhaohui S.Qin<zhaohui.qin@emory.edu>
}

\seealso{
	print.traseR
}
\examples{
	data(ASB)
	data(Tcell)
	x=traseR(snpdb=ASB,region=Tcell,rankby="pvalue",test.method="binomial",alternative ="greater")
	print(x)
	data(CEU)
	x=traseR(snpdb=ASB,region=Tcell,snpdb.bg=CEU,rankby="pvalue",test.method="binomial",alternative ="greater")
	print(x)
}




