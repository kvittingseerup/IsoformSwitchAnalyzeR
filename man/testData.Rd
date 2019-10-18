\name{exampleData}
\alias{exampleSwitchList}
\alias{exampleSwitchListIntermediary}
\alias{exampleSwitchListAnalyzed}
\docType{data}

\title{
Example data for IsoformSwitchAnalyzeR
}

\description{
Three switchAnalyzeRlist corresponding to a switchAnalyzeRlist in different stages of an isoform switch analyzer workflow.
}

\usage{
data("exampleSwitchList")

data("exampleSwitchListIntermediary")

data("exampleSwitchListAnalyzed")
}

\format{
see \code{?createSwitchAnalyzeRlist} for detailed format of an \code{switchAnalyzeRlist}
}

\details{
The three example switchAnalyzeRlist are:
\itemize{
    \item { \code{exampleSwitchList} : Which corresponds to a newly created \code{switchAnalyzeRlist} such as one would get by using either of the import* function (such as \code{importCufflinksData}) or by using \code{createSwitchAnalyzeRlist} on your own data. Not this is a small subset to allow for fast example generation.}
    \item { \code{exampleSwitchListIntermediary} : Which corresponds to the exampleSwitchList data (see above) which have been analyzed with the \code{isoformSwitchAnalysisPart1} function meaning that it have been filtered, tested for isoform switches, ORF have been predicted and both nucleotide and ORF amino acid sequences have been added to the \code{switchAnalyzeRlist}. Not this is a small subset to allow for fast example generation.}
    \item { \code{exampleSwitchListAnalyzed} : Which corresponds to a subset of two of the TCGA Cancer types analyzed in Vitting-Seerup et al 2017 which have been analyzed with the full switch analysis workflow (including external sequence analysis of protein domains (via Pfam), coding potential (via CPAT) and signal peptides (via SignalP)). Note that the nucleotide and amino acid sequences normally added to the switchAnalyzeRlist have been removed from the switchAnalyzeRlist (but also that they can easily be added again with the \link{extractSequence} function).}
}
}

\source{
exampleSwitchList and exampleSwitchListIntermediary is a modified subset of a dataset comparing human Embryonic Stemm Cells (hESC) vs induced Pluripotent Cells (iPS) and mature cells (Fibroblast) originally released with the cummeRbund R package. This data is only included to provide examples for usage of function. As it is modified to illustrate the package it should not be considered real and no biological conclusions should be made from it.

The exampleSwitchListAnalyzed is a subset of two of the TCGA Cancer types analyzed in Vitting-Seerup et al 2017 and are unmodified meaning results are real!
}

\references{
Vitting-Seerup et al. The Landscape of Isoform Switches in Human Cancers. Mol. Cancer Res. (2017).
}

\examples{
### Summarize newly created switchAnalyzeRList
data("exampleSwitchList")
summary(exampleSwitchList)
}
\keyword{datasets}
