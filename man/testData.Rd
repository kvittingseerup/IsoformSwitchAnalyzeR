\name{exampleData}
\alias{exampleSwitchList}
\alias{exampleSwitchListIntermediary}
\alias{exampleSwitchListAnalyzed}
\alias{isoformRepExpression}
\alias{isoformExonAnnoation}
\alias{designMatrix}
\docType{data}

\title{
Example data for IsoformSwitchAnalyzeR
}

\description{
The test data included in the IsoformSwitchAnalyzeR package is:
\itemize{
    \item{\code{1}: Three switchAnalyzeRlist corresponding to a switchAnalyzeRlist in different stages of an isoform switch analyzer workflow (see format for details).}
    \item{\code{2}: A set of standard R object which is needed to used the \code{importRdata} function.}
}
}

\usage{
data("exampleSwitchList")

data("exampleSwitchListIntermediary")

data("exampleSwitchListAnalyzed")

data('exampleRdata')
}

\format{
see \code{?createSwitchAnalyzeRlist} for detailed format of an \code{switchAnalyzeRlist}
}

\details{
The 'exampleRdata' loads the following objects, which are all obtained (and massaged) from the ballgown example data (see vignette('ballgown')).
\itemize{
\item { \code{isoformRepExpression} : Data.frame with replicate isoform expression.}
\item { \code{isoformExonAnnoation} : GRange object with the isoform annotation.}
\item { \code{designMatrix} : A design matrix for the data indicating two conditions 'WT' and 'KD'.}
}
More detials about the specification of each of these can be found in \code{?importRdata}.\cr

The three example switchAnalyzeRlist are:
\itemize{
    \item { \code{exampleSwitchList} : Which corresponds to a newly created \code{switchAnalyzeRlist} such as one would get by using either of the import* function (such as \code{importCufflinksData}) or by using \code{createSwitchAnalyzeRlist} on your own data. Not this is a small subset to allow for fast example generation.}
    \item { \code{exampleSwitchListIntermediary} : Which corresponds to the exampleSwitchList data (see above) which have been analyzed with the \code{isoformSwitchAnalysisPart1} function meaning that it have been filtered, tested for isoform switches, ORF have been predicted and both nucleotide and ORF amino acid sequences have been added to the \code{switchAnalyzeRlist}. Not this is a small subset to allow for fast example generation.}
    \item { \code{exampleSwitchListAnalyzed} : Which corresponds to the a full isoform switch analysis workflow (including external sequence analysis of protein domains (via Pfam), coding potential (via CPAT) and signal peptides (via SignalP)). This is the result of with the \code{isoformSwitchAnalysisPart1} function, performing the external sequence analysis and using with the \code{isoformSwitchAnalysisPart2} to run the rest of the pipeline. Note this is a larger subset to illustrate the combined analysis and that the nucleotide and amino acid sequences normally added to the switchAnalyzeRlist have been removed from the switchAnalyzeRlist to save space (they can easily be added again with the \link{extractSequence} function).}
}
}

\source{
This data is a modified subset of a dataset where two splice factors were knocked down. This data is only included to provide examples for usage of function. As it is modified it should not be considered real.
}

\examples{
### Summarize newly created switchAnalyzeRList
data("exampleSwitchList")
summary(exampleSwitchList)

### Visualize isoform switch in HOXC13 gene
data("exampleSwitchListAnalyzed")
switchPlot(exampleSwitchListAnalyzed, gene='HOXC13', condition1='Ctrl', condition2='KD1')

### Analyze general functional consequences pattern from fully analyzed switchAnalyzeRList
# visually
extractConsequenceSummary(exampleSwitchListAnalyzed, asFractionTotal = TRUE)

# Summarize as table instead
extractConsequenceSummary(exampleSwitchListAnalyzed, plot = FALSE, returnResult = TRUE)
}
\keyword{datasets}
