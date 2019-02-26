\name{extractConsequenceEnrichmentComparison}
\alias{extractConsequenceEnrichmentComparison}
\title{
Compare enrichment of specific consequences between comparisons
}
\description{
This function compares the enrichment of a consequnces (f.x. domain gain) between two comparisons (ctrl vs ko1 compared to ctrl vs ko2) and repports whether there is a signifcant difference between the comparisons. It other words it compares the output of \code{extractConsequenceEnrichment}.
}
\usage{
extractConsequenceEnrichmentComparison(
    switchAnalyzeRlist,
    consequencesToAnalyze = 'all',
    alpha=0.05,
    dIFcutoff = 0.1,
    countGenes = TRUE,
    analysisOppositeConsequence=FALSE,
    plot=TRUE,
    localTheme = theme_bw(base_size = 14),
    minEventsForPlotting = 10,
    returnResult=TRUE
)
}

\arguments{
  \item{switchAnalyzeRlist}{
A \code{switchAnalyzeRlist} object where \code{analyzeSwitchConsequences()} have been run to identify consequences of isoform switches
}
\item{consequencesToAnalyze}{
A string indicating which consequnces should be considered. See details for describtion (note it is idencital to the strings used with \code{analyzeSwitchConsequences}). Default is all consequences analyzed with \code{analyzeSwitchConsequences}
}
\item{alpha}{
The cutoff which the (callibrated) fdr correct p-values must be smaller than for calling significant switches. Default is 0.05.
}
\item{dIFcutoff}{
The cutoff which the changes in (absolute) isoform usage must be larger than before an isoform is considered switching. This cutoff can remove cases where isoforms with (very) low dIF values are deemed significant and thereby included in the downstream analysis. This cutoff is analogous to having a cutoff on log2 fold change in a normal differential expression analysis of genes to ensure the genes have a certain effect size. Default is 0.1 (10\%).
}

\item{countGenes}{
A logic indicating whether it is the number of genes (if TRUE) or isoform switches (if FALSE) which primarly result in gain/loss that are counted. Default is TRUE.
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

 \item{minEventsForPlotting}{
The minimum number of events (total gain/loss) must be pressent before the result is visualized. Default is 10.
}

  \item{returnResult}{
A logic indicating whether the analysis should be returned as a data.frame. Default is FALSE.
}
}

\details{
The significance test is performed with R's build in \code{prop.test()} with default parameters and resulting p-values are corrected via p.adjust() using FDR (Benjamini-Hochberg).
}

\value{
If \code{returnResult=TRUE} a data.frame with the statisitcal summary for each oposing consequences in each comparison. If \code{plot=TRUE} a plot summarizing the proportions is also created of switches with specific consequences is created.
}

\references{
\itemize{
    \item{Vitting-Seerup et al. The Landscape of Isoform Switches in Human Cancers. Mol. Cancer Res. (2017).}
    \item{Vitting-Seerup et al. IsoformSwitchAnalyzeR: Analysis of changes in genome-wide patterns of alternative splicing and its functional consequences. bioRxiv (2018).}
}
}
\author{
Kristoffer Vitting-Seerup
}
\seealso{
\code{\link{analyzeSwitchConsequences}}\cr
\code{\link{extractSwitchSummary}} \cr
\code{\link{extractConsequenceEnrichment}} \cr
\code{\link{extractConsequenceGenomeWide}}
}
\examples{
### Load exampled data
data("exampleSwitchListAnalyzed")

extractConsequenceEnrichmentComparison( exampleSwitchListAnalyzed)
}
