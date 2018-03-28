\name{extractSplicingEnrichmentComparison}
\alias{extractSplicingEnrichmentComparison}
\title{
Compare enrichment of specific type of alternative splicing between comparisons
}
\description{
This function compares the enrichment of alternative splicing (f.x. exon skipping) between two comparisons (ctrl vs ko1 compared to ctrl vs ko2) and repports whether there is a signifcant difference between the comparisons. It other words it compares the output of \code{extractSplicingEnrichment}.
}
\usage{
extractSplicingEnrichmentComparison(
    switchAnalyzeRlist,
    splicingToAnalyze = 'all',
    alpha = 0.05,
    dIFcutoff = 0.1,
    onlySigIsoforms = FALSE,
    plot = TRUE,
    localTheme = theme_bw(base_size = 14),
    returnResult=FALSE
)
}

\arguments{
  \item{switchAnalyzeRlist}{
A \code{switchAnalyzeRlist} object where \code{analyzeSwitchConsequences()} have been run to identify consequences of isoform switches
}
\item{splicingToAnalyze}{
A string indicating which consequnces should be considered. See details for description. Default is all.
}
\item{alpha}{
The cutoff which the (callibrated) fdr correct p-values must be smaller than for calling significant switches. Default is 0.05.
}
\item{dIFcutoff}{
The cutoff which the changes in (absolute) isoform usage must be larger than before an isoform is considered switching. This cutoff can remove cases where isoforms with (very) low dIF values are deemed significant and thereby included in the downstream analysis. This cutoff is analogous to having a cutoff on log2 fold change in a normal differential expression analysis of genes to ensure the genes have a certain effect size. Default is 0.1 (10\%).
}
\item{onlySigIsoforms}{ A logic indicating whether to only consider significant isoforms, meaning only analyzing genes where at least two isoforms which both have significant usage changes in opposite direction (quite strict) Naturally this only works if the isoform switch test used have isoform resolution (which the build in \link{isoformSwitchTest} has). If FALSE all isoforms with an absolute dIF value larger than \code{dIFcutoff} in a gene with significant switches (defined by \code{alpha} and \code{dIFcutoff}) are included in the pairwise comparison. Default is FALSE (non significant isoforms are also considered based on the logic that if one isoform changes it contribution - there must be an equivalent opposite change in usage in the other isoforms from that gene).
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
}

\details{
The classification of alternative splicing is always compared to the hypothetical pre-mRNA constructed by concatenating all exons from isoforms of the same gene. \cr

The alternative splicing types, which can be passed to \code{splicingToAnalyze} must be a combination of:
\itemize{
    \item \code{all} : All of the alternative splicing types indicated below.

    \item \code{IR} : Intron Retention.

    \item \code{A5} : Alternative 5' donor site (changes in the 5'end of the upstream exon).

    \item \code{A3} : Alternative 3' acceptor site (changes in the 3'end of the downstream exon).

    \item \code{ATSS} : Alternative Transcription Start Site.

    \item \code{ATTS} : Alternative Transcription Termination Site.

    \item \code{ES} : Exon Skipping.

    \item \code{MES} : Multiple Exon Skipping. Skipping of >1 consequtive exons.

    \item \code{MEE} : Mutually Exclusive Exons.
}

The significance test is performed with R's build in \code{fisher.test()} with default parameters and resulting p-values are corrected via p.adjust() using FDR (Benjamini-Hochberg).
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
\code{\link{analyzeAlternativeSplicing}}\cr
\code{\link{extractSplicingSummary}}\cr
\code{\link{extractSplicingEnrichment}}\cr
\code{\link{extractSplicingGenomeWide}}
}

\examples{
### Load example data
data("exampleSwitchListAnalyzed")

extractSplicingEnrichmentComparison( exampleSwitchListAnalyzed )
}
