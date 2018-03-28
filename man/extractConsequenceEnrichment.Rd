\name{extractConsequenceEnrichment}
\alias{extractConsequenceEnrichment}
\title{
Analyze data for enrichment of specific consequences
}
\description{
This functions analyzes for enrichment of specific consequnces by for each set of opposing consequnces (fx. domain gain vs loss), by analyzing the fraction of events belonging to one of them.
}
\usage{
extractConsequenceEnrichment(
    switchAnalyzeRlist,
    consequencesToAnalyze = 'all',
    alpha=0.05,
    dIFcutoff = 0.1,
    analysisOppositeConsequence=FALSE,
    plot=TRUE,
    localTheme = theme_bw(base_size = 12),
    returnResult=FALSE,
    returnSummary=TRUE
)
}

\arguments{
  \item{switchAnalyzeRlist}{
A \code{switchAnalyzeRlist} object where \code{analyzeSwitchConsequences()} have been run to identify consequences of isoform switches
}
\item{consequencesToAnalyze}{
A string indicating which consequnces should be considered. See detail section of \code{\link{analyzeSwitchConsequences}} for describtion . Default is all consequences analyzed with \code{analyzeSwitchConsequences}.
}
\item{alpha}{
The cutoff which the (callibrated) fdr correct p-values must be smaller than for calling significant switches. Default is 0.05.
}
\item{dIFcutoff}{
The cutoff which the changes in (absolute) isoform usage must be larger than before an isoform is considered switching. This cutoff can remove cases where isoforms with (very) low dIF values are deemed significant and thereby included in the downstream analysis. This cutoff is analogous to having a cutoff on log2 fold change in a normal differential expression analysis of genes to ensure the genes have a certain effect size. Default is 0.1 (10\%).
}
\item{analysisOppositeConsequence}{
A logic indicating whether reverse the analysis meaing if "Domain gains"" are analyze using default parameters setting \code{analysisOppositeConsequence=TRUE} will case the analysis to be performed on "Domain loss". The main effect is for the visual appearence of plot which will be mirrored (around the 0.5 fraction). Default is FALSE.
}
  \item{plot}{
A logic indicting whether the analysis should be plotted. Default is TRUE.
}
  \item{localTheme}{
General ggplo2 theme with which the plot is made, see \code{?ggplot2::theme} for more info. Default is theme_bw(base_size = 14).
}
  \item{returnResult}{
A logic indicating whether the analysis should be returned as a data.frame. Default is FALSE.
}
 \item{returnSummary}{
A logic indicating whether to return the statistical summary (if TRUE) or the underlying data (if FALSE). Default is TRUE.
}
}

\details{
The significance test is performed with R's build in \code{prop.test()} with default parameters and resulting p-values are corrected via p.adjust() using FDR (Benjamini-Hochberg).
}

\value{
If \code{returnResult=TRUE} a data.frame with the statisitcal summary for each oposing consequences in each comparison. If \code{plot=TRUE} a plot summarizing the proportions is also created of switches with specific consequences is created.
}

\references{
Vitting-Seerup et al. The Landscape of Isoform Switches in Human Cancers. Mol. Cancer Res. (2017).
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
