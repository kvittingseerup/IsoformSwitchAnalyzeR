\name{switchPlotTopSwitches}
\alias{switchPlotTopSwitches}

\title{
Creating the Isoform Switch Analysis Plot for the Top Switches
}
\description{
This function outputs the top n (defined by \code{n}) Isoform Switch Analysis Plot (see \link{switchPlot}) for genes with significant isoform switches (as defined by \code{alpha} and \code{dIFcutoff}) to a specific folder (controlled by \code{pathToOutput}. The plots are automatically sorted by decreasing significance or switch size (as controlled by \code{sortByQvals}). The plots can furthermore be created in sub-folders based both which conditions are compared and whether any consequences of the switch have been predicted. In summary it facilitates an easy and prioritized, (but comprehensive), manual analysis of isoform switches.
}
\usage{
switchPlotTopSwitches(
    switchAnalyzeRlist,
    alpha = 0.05,
    dIFcutoff = 0.1,
    n=10,
    sortByQvals=TRUE,
    filterForConsequences = FALSE,
    pathToOutput = getwd(),
    splitComparison=TRUE,
    splitFunctionalConsequences = TRUE,
    IFcutoff=0.05,
    fileType = "pdf",
    additionalArguments=list(),
    quiet=FALSE
)

}
\arguments{
  \item{switchAnalyzeRlist}{
A \code{switchAnalyzeRlist} containing all the annotation for the isoforms.
}

  \item{alpha}{
The cutoff which the FDR correct p-values must be smaller than for calling significant switches. Default is 0.05.
}

\item{dIFcutoff}{
The cutoff which the changes in (absolute) isoform usage must be larger than before an isoform is considered switching. This cutoff can remove cases where isoforms with (very) low dIF values are deemed significant and thereby included in the downstream analysis. This cutoff is analogous to having a cutoff on log2 fold change in a normal differential expression analysis of genes to ensure the genes have a certain effect size. Default is 0.1 (10\%).
}

\item{n}{
The number of top genes (after filtering and sorted according to \code{sortByQvals}) that should be generated in each sub-folder indicated by \code{splitComparison} and \code{splitFunctionalConsequences}. Use Inf to create all (NA will internally be converted to Inf for backward comparability). Default is 10.
}

\item{sortByQvals}{
A logic indicating whether to the top \code{n} features are sorted by decreasing significance (increasing q-values) (if \code{sortByQvals=TURE}) or decreasing switch size (absolute dIF, which are still significant as defined by \code{alpha}) (if \code{sortByQvals=FALSE}). The dIF values for genes are considered as the total change within the gene calculated as \code{sum(abs(dIF))} for each gene. Default is TRUE (sort by p-values).
}

  \item{filterForConsequences}{
A logic indicating whether to only plot gene with predicted consequences of the isoform switch. Requires that predicted consequences have been annotated (via \link{analyzeSwitchConsequences}. Default is FALSE.
}

  \item{pathToOutput}{
A path to the folder in which the plots should be made. Default is working directory ( getwd() ).
}

\item{splitComparison}{
A logic indicating whether to create a sub-folder for each comparison. If splitComparison is TRUE the sub-folders will be created else all isoform switch analyzer plots will saved in the same folder. Default is TRUE.
}

  \item{splitFunctionalConsequences}{
A logic indicating whether to create a sub-folder for those switches with predicted consequences and another sub-folder for those without. Requires that \link{analyzeSwitchConsequences} have been run. If \code{splitComparison=TRUE} the sub-folders from this argument will be created within the comparison sub-folders. Default is TRUE.
}

\item{IFcutoff}{
The cutoff used for the minimum contribution to gene expression (in at least one condition) an isoforms must have to be plotted (measured as Isoform Fraction (IF) values). Default is 0 (which removes isoforms not contributing in any of the conditions).
}

  \item{fileType}{
A string indicating which file type is generated. Available are options are \'pdf\' and \'png\'. Default is pdf.
}

  \item{additionalArguments}{
A named list arguments passed to the \code{switchPlot} function which creates the individual Isoform Switch Analysis Plots. The name of the list entries must correspond to the corresponding argument in the \code{switchPlot} function.
}

\item{quiet}{ A logic indicating whether to avoid printing progress messages. Default is FALSE}
}


\details{
Changes in isoform usage are measure as the difference in isoform fraction (dIF) values, where isoform fraction (IF) values are calculated as <isoform_exp> / <gene_exp>.\cr

For a list of the top switching genes see ?extractTopSwitches.
}

\value{
An Isoform Switch Analysis Plot (as produce by \code{switchPlot}) for each of the top n switches in each comparison where a gene have a significant isoform switch is generated in the folder supplied by \code{pathToOutput}
}
\references{
Vitting-Seerup et al. The Landscape of Isoform Switches in Human Cancers. Mol. Cancer Res. (2017).
}
\author{
Kristoffer Vitting-Seerup
}
\seealso{
\code{\link{switchPlot}}\cr
\code{\link{analyzeSwitchConsequences}}
}
