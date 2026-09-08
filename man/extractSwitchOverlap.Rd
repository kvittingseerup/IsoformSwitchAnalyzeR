\name{extractSwitchOverlap}
\alias{extractSwitchOverlap}

\title{
Visualize Switch Overlap
}
\description{
This function produces two Venn diagrams respectively showing the overlap in switching isoforms and genes.
}
\usage{
extractSwitchOverlap(
    switchAnalyzeRlist,
    filterForConsequences = FALSE,
    alpha = 0.05,
    dIFcutoff = 0.1,
    scaleVennIfPossible=TRUE,
    plotIsoforms = TRUE,
    plotSwitches = TRUE,
    plotGenes = TRUE
)
}

\arguments{
\item{switchAnalyzeRlist}{
A \code{switchAnalyzeRlist} object.
}
\item{filterForConsequences}{
A logical indicating whether to filter for genes with functional consequences. Requires that analyzeSwitchConsequences() have been run on the switchAnalyzeRlist. The output will then be the number of significant genes and isoforms originating from genes with predicted consequences. Default is FALSE.
}
  \item{alpha}{
The cutoff which the FDR correct p-values must be smaller than for calling significant switches. Default is 0.05.
}
\item{dIFcutoff}{
The cutoff which the changes in (absolute) isoform usage must be larger than before an isoform is considered switching. This cutoff can remove cases where isoforms with (very) low dIF values are deemed significant and thereby included in the downstream analysis. This cutoff is analogous to having a cutoff on log2 fold change in a normal differential expression analysis of genes to ensure the genes have a certain effect size. Default is 0.1 (10\%).
}
\item{scaleVennIfPossible}{
A logic indicating whether the Venn diagram should be scaled (so the circle area and overlap size reflect the number of features) if possible. Only available for 2- and 3-way Venn Diagrams. Default is TRUE.
}

\item{plotIsoforms}{
A logic indicating whether the Venn diagram of differentially used isoforms should be plotted. Default is TRUE.
}

\item{plotSwitches}{
A logic indicating whether the Venn diagram of identified isoform switches should be plotted. Default is TRUE.
}

\item{plotGenes}{
A logic indicating whether the Venn diagram of genes containing differentially used isoforms should be plotted. Default is TRUE.
}

}

\value{
A Venn diagram which shows the number of isoforms and genes with a isoform switch.
}

\references{
\itemize{
    \item{Han C, Gilis J, Iriondo Delgado E, Clement L, Vitting-Seerup K. IsoformSwitchAnalyzeR v2: analysis of functional isoform changes in long-read and single-cell sequencing data. NAR Genom Bioinform (2026). doi:10.1093/nargab/lqag098.}
}
}

\author{
Kristoffer Vitting-Seerup
}
\seealso{
\code{\link{preFilter}}\cr
\code{\link{isoformSwitchTestDEXSeq}}\cr
\code{\link{isoformSwitchTestSatuRn}}\cr
\code{\link{extractTopSwitches}}\cr
\code{\link{extractSwitchSummary}}\cr
\code{\link{analyzeSwitchConsequences}}
}

\examples{
# Load example data and prefilter
data("exampleSwitchListAnalyzed")

extractSwitchOverlap(exampleSwitchListAnalyzed)

}
