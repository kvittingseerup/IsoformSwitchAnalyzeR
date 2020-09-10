\name{extractSplicingEnrichment}
\alias{extractSplicingEnrichment}
\title{
Analyze data for enrichment of specific type of alternative splicing
}
\description{
This functions function analyzes (the number of and) enrichment of specific splice events by for each set of opposing event (e.g.. exon skipping gain vs loss), by analyzing the fraction of events belonging to each type of consequence. Please note this summarizes the differences between the isoforms in a switch - for an overview of the total number of AS events please use \link{extractSplicingSummary}.
}
\usage{
extractSplicingEnrichment(
    switchAnalyzeRlist,
    splicingToAnalyze = 'all',
    alpha = 0.05,
    dIFcutoff = 0.1,
    onlySigIsoforms = FALSE,
    countGenes = TRUE,
    plot = TRUE,
    localTheme = theme_bw(base_size = 14),
    minEventsForPlotting = 10,
    returnResult=TRUE,
    returnSummary=TRUE
)
}

\arguments{
  \item{switchAnalyzeRlist}{
A \code{switchAnalyzeRlist} object where \code{analyzeSwitchConsequences()} have been run to identify consequences of isoform switches
}
\item{splicingToAnalyze}{
A string indicating which consequences should be considered. See details for description. Default is all.
}
\item{alpha}{
The cutoff which the FDR correct p-values must be smaller than for calling significant switches. Default is 0.05.
}
\item{dIFcutoff}{
The cutoff which the changes in (absolute) isoform usage must be larger than before an isoform is considered switching. This cutoff can remove cases where isoforms with (very) low dIF values are deemed significant and thereby included in the downstream analysis. This cutoff is analogous to having a cutoff on log2 fold change in a normal differential expression analysis of genes to ensure the genes have a certain effect size. Default is 0.1 (10\%).
}
\item{onlySigIsoforms}{ A logic indicating whether to only consider significant isoforms, meaning only analyzing genes where at least two isoforms which both have significant usage changes in opposite direction (quite strict) Naturally this only works if the isoform switch test used have isoform resolution (which the build in \link{isoformSwitchTestDEXSeq} has). If FALSE all isoforms with an absolute dIF value larger than \code{dIFcutoff} in a gene with significant switches (defined by \code{alpha} and \code{dIFcutoff}) are included in the pairwise comparison. Default is FALSE (non significant isoforms are also considered based on the logic that if one isoform changes it contribution - there must be an equivalent opposite change in usage in the other isoforms from that gene).
}

\item{countGenes}{
A logic indicating whether it is the number of genes (if TRUE) or isoform switches (if FALSE) which primary result in gain/loss that are counted. Default is TRUE.
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
A logic indicating whether to return the statistical summary (if TRUE) or the underlying data (if FALSE). If FALSE (and \code{plot=TRUE}) the ggplot2 object will be returned instead. Default is TRUE.
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

    \item \code{ES} : Exon Skipping (EI means Exon Inclusion).

    \item \code{MES} : Multiple Exon Skipping. Skipping of >1 consecutive exons. (MEI means Multiple Exon Inclusion).

    \item \code{MEE} : Mutually Exclusive Exons.
}

For details of how to interpret the splice events see the \code{details} section of \code{\link{analyzeAlternativeSplicing}}. \cr


The significance test is performed with R's build in \code{prop.test()} with default parameters and resulting p-values are corrected via p.adjust() using FDR (Benjamini-Hochberg).
}

\value{
If \code{plot=TRUE} a plot summarizing the proportions is also created of switches with specific consequences is created.


If \code{returnResult=TRUE} a data.frame with the statistical summary for each opposing consequences in each comparison. This data.frame will have the following collumns:
\itemize{
    \item{condition_1: Condition 1.}
    \item{condition_2: Condition 2.}
    \item{AStype: The type of splicing considered.}
    \item{nUp: The number of genes with a gain of the splicing type described in the AStype column.}
    \item{nDown: The number of genes with a loss of the splicing type described in the AStype column.}
    \item{propUp: Proportion of total number of genes (of genes with either loss or gain of the splice type described in the AStype column) being having a gain.}
    \item{propUpCiLo: The lower boundary of the confidence interval of the propUp.}
    \item{propUpCiLo: The high boundary of the confidence interval of the propUp.}
    \item{propUpPval: The p-value associated with the null hypothesis that propUp is 0.5.}
    \item{propUpQval: The q-values resulting when p-values are corrected via p.adjust() using FDR (Benjamini-Hochberg).}
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
\code{\link{isoformSwitchTestDEXSeq}}\cr
\code{\link{isoformSwitchTestDRIMSeq}}\cr
\code{\link{analyzeAlternativeSplicing}}\cr
\code{\link{extractSplicingSummary}}\cr
\code{\link{extractSplicingEnrichmentComparison}}\cr
\code{\link{extractSplicingGenomeWide}}
}

\examples{
### Load example data
data("exampleSwitchListAnalyzed")

extractSplicingEnrichment( exampleSwitchListAnalyzed )
}
