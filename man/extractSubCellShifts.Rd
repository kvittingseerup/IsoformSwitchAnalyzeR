\name{extractSubCellShifts}
\alias{extractSubCellShifts}
\title{
Global overview of sub-cellular location changes
}
\description{
Extract summary of how sub-celluar locations have changed due to isoform switches.
}
\usage{
extractSubCellShifts(
    switchAnalyzeRlist,
    plotGenes = TRUE,
    returnResult = FALSE,
    localTheme = theme_bw()
)
}

\arguments{
  \item{switchAnalyzeRlist}{
A \code{switchAnalyzeRlist} object where \code{analyzeSwitchConsequences()} have been run to identify consequences of isoform switches
}

\item{plotGenes}{
A logic indicating whether to plot number of genes (if TRUE) or isoform switches (if FALSE). Default is TRUE.
}

  \item{returnResult}{
A logic indicating whether the analysis should be returned as a data.frame. If FALSE a ggplot2 object will be returned instead. Default is FALSE.
}

  \item{localTheme}{
General ggplot2 theme with which the plot is made, see \code{?ggplot2::theme} for more info. Default is theme_bw().
}

}


\value{
If \code{returnResult=FALSE} a plot summarizing switches.

If \code{returnResult=TRUE} a data.frame with the summary for each comparison. This data.frame will have the following collumns:
\itemize{
    \item{condition_1: Condition 1.}
    \item{condition_2: Condition 2.}

    \item{location_gain: The location gained. Is "None if no location was lost"}
    \item{location_loss: The location Lost. Is "None if no location was lost"}

    \item{n_switch: The number of isoform switches with the location change.}
    \item{n_genes: The number of genes containing isoform switche(s) with the location change.}

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
}
