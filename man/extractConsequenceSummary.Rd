\name{extractConsequenceSummary}
\alias{extractConsequenceSummary}
\title{
Analyze Switch Consequences
}
\description{
This functions function summarizes the individual types of consequences for each gene or the pairwise switches and plots and/or returns a data.frame with the information
}
\usage{
extractConsequenceSummary(
    switchAnalyzeRlist,
    consequencesToAnalyze='all',
    includeCombined=FALSE,
    asFractionTotal=FALSE,
    alpha=0.05,
    dIFcutoff=0.1,
    plot=TRUE,
    plotGenes=FALSE,
    simplifyLocation = TRUE,
    removeEmptyConsequences = FALSE,
    localTheme=theme_bw(),
    returnResult=FALSE
)
}

\arguments{
  \item{switchAnalyzeRlist}{
A \code{switchAnalyzeRlist} object where \code{analyzeSwitchConsequences()} have been run to identify consequences of isoform switches
}
\item{consequencesToAnalyze}{
A string indicating which consequences should be considered. See detail section of \code{\link{analyzeSwitchConsequences}} for description . Default is all consequences analyzed with \code{analyzeSwitchConsequences}.
}
\item{includeCombined}{
A logic indicating whether an analysis of how many (how large a fraction) of genes have any type of functional consequence.
}
  \item{asFractionTotal}{
A logic indicating whether the consequences should be summarized calculated as numbers (if FALSE) or as a fraction of the total number of switches/genes (as indicated by \code{plotGenes}). Default is FALSE.
}
\item{alpha}{
The cutoff which the FDR correct p-values must be smaller than for calling significant switches. Default is 0.05.
}
\item{dIFcutoff}{
The cutoff which the changes in (absolute) isoform usage must be larger than before an isoform is considered switching. This cutoff can remove cases where isoforms with (very) low dIF values are deemed significant and thereby included in the downstream analysis. This cutoff is analogous to having a cutoff on log2 fold change in a normal differential expression analysis of genes to ensure the genes have a certain effect size. Default is 0.1 (10\%).
}
  \item{plot}{
A logic indicting whether the analysis should be plotted. If TRUE and \code{returnResult = FALSE} the ggplot2 object will be returned instead. Default is TRUE.
}
  \item{plotGenes}{
A logic indicating whether to plot the number/fraction of genes with (if TRUE) or isoforms (if FALSE) involved with isoform switches with functional consequences (both filtered via \code{alpha} and \code{dIFcutoff}.
}

\item{simplifyLocation}{
A logic indicating whether to simplify the switches involved in changes in sub-cellular localizations (due the the hundreds of possible combinations). Done by only considering where the isoform used more has a location switch to. Default is TRUE.
}

\item{removeEmptyConsequences}{
A logic indicating whether to remove consequenses analyzed but where no differences was found (those showing zero in the plot). Default is FALSE.
}

  \item{localTheme}{
General ggplo2 theme with which the plot is made, see \code{?ggplot2::theme} for more info. Default is theme_bw().
}
  \item{returnResult}{
A logic indicating whether the summarized results should be returned as a data.frame. If FALSE (and \code{plot=TRUE}) the ggplot2 object will be returned instead. Default is FALSE.
}
}

\details{
A less detailed version just summarizing the number of switches with functional consequences can be obtained by setting \code{filterForConsequences=TRUE} in the \code{extractSwitchSummary} function.

For details on the arguments passed to \code{consequencesToAnalyze} please see details section of \link{analyzeSwitchConsequences}.
}

\value{
If \code{returnResult=TRUE} a data.frame with the number (and fraction) of switches with specific consequences in each condition is returned. If \code{plot=TRUE} a plot summarizing the number (or fraction) of switches with specific consequences is created.
}

\references{
Vitting-Seerup et al. The Landscape of Isoform Switches in Human Cancers. Mol. Cancer Res. (2017).
}
\author{
Kristoffer Vitting-Seerup
}
\seealso{
\code{\link{analyzeSwitchConsequences}}\cr
\code{\link{extractConsequenceEnrichment}}\cr
\code{\link{extractConsequenceEnrichmentComparison}} \cr
\code{\link{extractConsequenceGenomeWide}}
}
\examples{
### Prepare example data
data("exampleSwitchListAnalyzed")

### Summarize switch consequences
consequenceSummary <- extractConsequenceSummary(
    exampleSwitchListAnalyzed,
    returnResult = TRUE,        # return data.frame with summary
    plotGenes = TRUE            # plot summary
)

dim(consequenceSummary)

subset(consequenceSummary, featureCompared=='Domains identified')

}
