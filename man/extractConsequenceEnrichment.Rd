\name{extractConsequenceEnrichment}
\alias{extractConsequenceEnrichment}
\title{
Analyze data for enrichment of specific consequences
}
\description{
This functions analyzes the genome wide enrichment of specific consequences by for each set of opposing consequences (e.g.. domain gain vs loss) analyzing the fraction of events belonging to one of them.
}
\usage{
extractConsequenceEnrichment(
    switchAnalyzeRlist,
    consequencesToAnalyze = 'all',
    alpha=0.05,
    dIFcutoff = 0.1,
    countGenes = TRUE,
    analysisOppositeConsequence=FALSE,
    plot=TRUE,
    localTheme = theme_bw(base_size = 12),
    minEventsForPlotting = 10,
    returnResult=TRUE,
    returnSummary=TRUE
)
}

\arguments{
  \item{switchAnalyzeRlist}{
A \code{switchAnalyzeRlist} object where \code{analyzeSwitchConsequences()} have been run to identify consequences of isoform switches
}
\item{consequencesToAnalyze}{
A string indicating which consequences should be considered. See detail section of \code{\link{analyzeSwitchConsequences}} for description . Default is all consequences analyzed with \code{analyzeSwitchConsequences}.
}
\item{alpha}{
The cutoff which the FDR correct p-values must be smaller than for calling significant switches. Default is 0.05.
}
\item{dIFcutoff}{
The cutoff which the changes in (absolute) isoform usage must be larger than before an isoform is considered switching. This cutoff can remove cases where isoforms with (very) low dIF values are deemed significant and thereby included in the downstream analysis. This cutoff is analogous to having a cutoff on log2 fold change in a normal differential expression analysis of genes to ensure the genes have a certain effect size. Default is 0.1 (10\%).
}

\item{countGenes}{
A logic indicating whether it is the number of genes (if TRUE) or isoform switches (if FALSE) which primary result in gain/loss that are counted. Default is TRUE.
}

\item{analysisOppositeConsequence}{
A logic indicating whether reverse the analysis meaning if "Domain gains"" are analyze using default parameters setting \code{analysisOppositeConsequence=TRUE} will case the analysis to be performed on "Domain loss". The main effect is for the visual appearance of plot which will be mirrored (around the 0.5 fraction). Default is FALSE.
}
  \item{plot}{
A logic indicting whether the analysis should be plotted. If TRUE and \code{returnResult = FALSE} the ggplot2 object will be returned instead. Default is TRUE.
}
  \item{localTheme}{
General ggplo2 theme with which the plot is made, see \code{?ggplot2::theme} for more info. Default is theme_bw(base_size = 14).
}

 \item{minEventsForPlotting}{
The minimum number of events (total gain/loss) must be present before the result is visualized. Default is 10.
}

  \item{returnResult}{
A logic indicating whether the analysis should be returned as a data.frame. If FALSE (and \code{plot=TRUE}) the ggplot2 object will be returned instead. Default is TRUE.
}
 \item{returnSummary}{
A logic indicating whether to return the statistical summary (if TRUE) or the underlying data (if FALSE). Depends on \code{returnResult=TRUE} Default is TRUE.
}
}

\details{
The significance test is performed with R's build in \code{binom.test()} with default parameters and resulting p-values are corrected via p.adjust() using FDR (Benjamini-Hochberg).
}

\value{
If \code{plot=TRUE} a plot summarizing the proportions is also created of switches with specific consequences is created.

If \code{returnResult=TRUE} a data.frame with the statistical summary for each opposing consequences in each comparison. This data.frame will have the following collumns:
\itemize{
    \item{condition_1: Condition 1.}
    \item{condition_2: Condition 2.}
    \item{conseqPair: The set of oposite consequences consudedered.}
    \item{feature: Description of which consequence the calculations are done from the perspective of (with the opposite mention in parentheses)}
    \item{propOfRelevantEvents: Proportion of total number of genes (of genes affected by the consequence described in the conseqPair column) being of the type described in the feature column.}
    \item{propCiLo: The lower boundary of the confidence interval of the propOfRelevantEvents.}
    \item{propCiHi: The high boundary of the confidence interval of the propOfRelevantEvents.}
    \item{propPval: The p-value associated with the null hypothesis that propOfRelevantEvents is 0.5.}
    \item{nUp: The number of genes with the consequence described in the feature column.}
    \item{nDown: The number of genes with the opposite consequence of what is described in the feature column (as described in the parentheses of the feature column.}
    \item{propQval: The q-values resulting when p-values are corrected via p.adjust() using FDR (Benjamini-Hochberg).}
}
}

\references{
\itemize{
    \item{Vitting-Seerup et al. The Landscape of Isoform Switches in Human Cancers. Mol. Cancer Res. (2017).}
    \item{Vitting-Seerup et al. IsoformSwitchAnalyzeR: Analysis of changes in genome-wide patterns of alternative splicing and its functional consequences. Bioinformatics (2019).}
}
}
\author{
Kristoffer Vitting-Seerup
}
\seealso{
\code{\link{analyzeSwitchConsequences}}\cr
\code{\link{extractSwitchSummary}} \cr
\code{\link{extractConsequenceEnrichmentComparison}} \cr
\code{\link{extractConsequenceGenomeWide}}
}
\examples{
### Load exampled data
data("exampleSwitchListAnalyzed")

extractConsequenceEnrichment( exampleSwitchListAnalyzed)
}
