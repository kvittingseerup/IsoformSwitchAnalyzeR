\name{extractGeneExpression}
\alias{extractGeneExpression}
\title{
Extract gene counts or abundances from a switchAnalyzeRlist object.
}
\description{
Extract replicate gene counts or expression from a switchAnalyzeRlist object using all the annotation fixes employed in creating the switchAnalyzeRlist.
}
\usage{
extractGeneExpression(
    switchAnalyzeRlist,
    extractCounts = TRUE
)
}

\arguments{
\item{switchAnalyzeRlist}{
A \code{switchAnalyzeRlist} object.
}

\item{extractCounts}{
A logic to indicate whether to extract counts (if TRUE, default) or expression estimates (if FALSE).
}

}

\value{
A data.frame with the replicate count/abundance estimates.
}

\references{
Vitting-Seerup et al. The Landscape of Isoform Switches in Human Cancers. Mol. Cancer Res. (2017).
}
\author{
Kristoffer Vitting-Seerup
}

\examples{
data("exampleSwitchList")

geneCountMatrix <- extractGeneExpression(
    exampleSwitchList,
    extractCounts = TRUE
)

geneExpresionMatrix <- extractGeneExpression(
    exampleSwitchList,
    extractCounts = FALSE
)
}
