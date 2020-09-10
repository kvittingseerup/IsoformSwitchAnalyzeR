\name{extractSplicingGenomeWide}
\alias{extractSplicingGenomeWide}

\title{
Genome wide Analysis of alternative splicing
}

\description{
This function enables a genome wide analysis of changes in isoform usage of isoforms with a common annotation.

Specifically this function extract isoforms of interest and for each splicing type (such as exon skipping) the global distribution of IF (measuring isoform usage) are plotted for each subset of features in that category (e.g with exons skipping vs without exon skipping). This enables a global analysis of isoforms with a common annotation.

The isoforms of interest can either be defined by isoforms form gene differentially expressed, isoform that are differential expressed or isoforms from genes with isoform switching - as controlled by \code{featureToExtract}. Please note that the \link{extractSplicingEnrichment} function probably more relevant than using \code{featureToExtract='isoformUsage'} since it directly uses the paired information from switches.

This function offers both visualization of the result as well as analysis via summary statistics of the comparisons.
}

\usage{
extractSplicingGenomeWide(
    switchAnalyzeRlist,
    featureToExtract = 'isoformUsage',
    splicingToAnalyze = 'all',
    alpha=0.05,
    dIFcutoff = 0.1,
    log2FCcutoff = 1,
    violinPlot=TRUE,
    alphas=c(0.05, 0.001),
    localTheme=theme_bw(),
    plot=TRUE,
    returnResult=TRUE
)
}

\arguments{
\item{switchAnalyzeRlist}{A \code{switchAnalyzeRlist} object containing the result of an isoform switch analysis (such as the one provided by \code{isoformSwitchTestDEXSeq()}) as well as additional annotation data for the isoforms.}

\item{featureToExtract}{
This argument, given as a string, defines the set isoforms which should be analyzed. The available options are:
\itemize{
    \item {\code{'isoformUsage'} (Default): Analyze a subset of isoforms defined by change in isoform usage (controlled by \code{dIFcutoff}) and the significance of the change in isoform expression (controlled by \code{alpha}). Please note that the \link{extractSplicingEnrichment} function probably more relevant than using \code{featureToExtract='isoformUsage'} since it directly uses the paired information from switches.}
    \item {\code{'isoformExp'} :Analyze a subset of isoforms defined by change in isoform expression (controlled by \code{log2FCcutoff}) and the significance of the change in isoform expression (controlled by \code{alpha})}
    \item {\code{'geneExp'} :Analyze all isoforms from a subset of genes defined by by change in gene expression (controlled by \code{log2FCcutoff}) and the significance of the change in gene expression (controlled by \code{alpha})}
    \item {\code{'all'} : Analyze all isoforms stored in the switchAnalyzeRlist (note that this is highly depending on the parameter \code{reduceToSwitchingGenes} in \link{isoformSwitchTestDEXSeq} - which should be set to FALSE (default is TRUE) if the 'all' option should be used here).}
}
}

\item{splicingToAnalyze}{
A string indicating which consequences should be considered. See details for description. Default is all.
}

\item{alpha}{The cutoff which the FDR correct p-values (q-values) must be smaller than for calling significant switches. Default is 0.05.}

\item{dIFcutoff}{The cutoff which the changes in (absolute) isoform usage must be larger than before an isoform is considered switching. This cutoff can remove cases where isoforms with (very) low dIF values are deemed significant and thereby included in the downstream analysis. This cutoff is analogous to having a cutoff on log2 fold change in a normal differential expression analysis of genes to ensure the genes have a certain effect size. Default is 0.1 (10\%).}

\item{log2FCcutoff}{The cutoff which the changes in (absolute) isoform or gene expression must be larger than before an isoform is considered for inclusion.}

\item{violinPlot}{
A logical indicating whether to make a violin plots (if TRUE) or boxplots (if FALSE). Violin plots will always have added 3 black dots, one of each of the 25th, 50th (median) and 75th percentile of the data. Default is TRUE.
}

\item{alphas}{
A numeric vector of length two giving the significance levels represented in plots. The numbers indicate the q-value cutoff for significant (*) and highly significant (***) respectively. Default 0.05 and 0.001 which should be interpret as q<0.05 and q<0.001 respectively). If q-values are higher than this they will be annotated as 'ns' (not significant).
}

  \item{localTheme}{
General ggplo2 theme with which the plot is made, see \code{?ggplot2::theme} for more info. Default is theme_bw().
}

\item{plot}{
A logic indicting whether the analysis should be plotted. If TRUE and \code{returnResult = FALSE} the ggplot2 object will be returned instead. Default is TRUE.
}

\item{returnResult}{
A logical indicating whether to return a data.frame with summary statistics of the comparisons (if TRUE) or not (if FALSE). If FALSE (and \code{plot=TRUE}) the ggplot2 object will be returned instead. Default is TRUE.
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

    \item \code{MES} : Multiple Exon Skipping. Skipping of >1 consecutive exons.

    \item \code{MEE} : Mutually Exclusive Exons.
}


The significance test is performed with R's build in \code{wilcox.test()} (aka 'Mann-Whitney-U') with default parameters and resulting p-values are corrected via p.adjust() using FDR (Benjamini-Hochberg).
}

\value{
If \code{plot=TRUE}: A plot of the distribution of IF values as a function of the annotation and condition compared.
If \code{returnResult=TRUE}: A data.frame with the summary statistics from the comparison of the two conditions with a Wilcox.test.
}

\references{
Vitting-Seerup et al. The Landscape of Isoform Switches in Human Cancers. Mol. Cancer Res. (2017).
}

\author{
Kristoffer Vitting-Seerup
}

\seealso{
\code{\link{isoformSwitchTestDEXSeq}}\cr
\code{\link{isoformSwitchTestDRIMSeq}}\cr
\code{\link{analyzeAlternativeSplicing}}\cr
\code{\link{extractSplicingSummary}}\cr
\code{\link{extractSplicingEnrichment}}\cr
\code{\link{extractSplicingEnrichmentComparison}}
}

\examples{
### Load example data
data("exampleSwitchListAnalyzed")

extractSplicingGenomeWide( exampleSwitchListAnalyzed )
}
