\name{extractTopSwitches}
\alias{extractTopSwitches}

\title{
Extract Top Isoform Switches.
}
\description{
This function allows the user extract the (top) switching genes/isoforms (with functional consequences).
}
\usage{
extractTopSwitches(
    switchAnalyzeRlist,
    filterForConsequences=FALSE,
    extractGenes=TRUE,
    alpha=0.05,
    dIFcutoff = 0.1,
    n=10,
    inEachComparison=FALSE,
    sortByQvals=TRUE
)
}

\arguments{
\item{switchAnalyzeRlist}{
A \code{switchAnalyzeRlist} object.
}

\item{extractGenes}{
A logic indicating whether to extract the (top) switching isoforms (if FALSE) or top swiching genes (if TRUE). Default is TRUE (extract genes).
}
\item{filterForConsequences}{
A logical indicating whether to filter for genes with functional consequences. Requtires that analyzeSwitchConsequences() have been run on the switchAnalyzeRlist. Default is FALSE.
}

\item{alpha}{
The cutoff which the (callibrated) fdr correct p-values must be smaller than for calling significant switches. Defualit is 0.05.
}
\item{dIFcutoff}{
The cutoff which the changes in (absolute) isoform usage must be larger than before an isoform is considered switching. This cutoff can remove cases where isoforms with (very) low dIF values are deemed significant and thereby included in the downstream analysis. This cutoff is analogous to having a cutoff on log2 fold change in a normal differential expression analysis of genes to ensure the genes have a certain effect size. Default is 0.1 (10\%).
}
\item{n}{
The number of switching features (genes/isoforms) to return. Use Inf to return all significant results (NA will internally be conveted to Inf for backward compatability). Default is 10.
}
\item{inEachComparison}{
A logic indicating whether to extract top n in each comparison (if TRUE) or from the all analysis (if FALSE). Default is FALSE.
}
\item{sortByQvals}{
A logic indicating whether to the top \code{n} features are defined be smallest q-values (if \code{sortByQvals=TURE}) or the largest changes in isoform usage (absolute dIF) which are still significant (if \code{sortByQvals=FALSE}). The dIF values for genes are considered as the total change within the gene calculated as \code{sum(abs(dIF))} for each gene. If set to NA no sorting is perfomed. Default is TRUE (sort by p-values).
}
}

\value{
A \code{data.frame} containing the top \code{n} switching genes or isoforms as controlled  by the \code{extractGenes} argument, sorted by q-values or dIF values as controlled  by the \code{sortByQvals} argument.
}

\references{
Vitting-Seerup et al. The Landscape of Isoform Switches in Human Cancers. Mol. Cancer Res. (2017).
}

\author{
Kristoffer Vitting-Seerup
}
\seealso{
\code{\link{preFilter}}\cr
\code{\link{isoformSwitchTestDEXSeq}}\cr
\code{\link{isoformSwitchTestDRIMSeq}}\cr
\code{\link{analyzeSwitchConsequences}}
}
\examples{
# Load example data and prefilter
data("exampleSwitchList")
exampleSwitchList <- preFilter(exampleSwitchList)

# Perfom test
exampleSwitchListAnalyzed <- isoformSwitchTestDEXSeq(exampleSwitchList)

# extract summary of number of switching features
extractSwitchSummary(exampleSwitchListAnalyzed)


### Filter for functional consequences (identified via analyzeSwitchConsequences() )
data("exampleSwitchListAnalyzed")
switchingIso <- extractTopSwitches(
    exampleSwitchListAnalyzed,
    filterForConsequences = TRUE,
)
dim(switchingIso)
head(switchingIso,2)

}
