\name{subsetSwitchAnalyzeRlist}
\alias{subsetSwitchAnalyzeRlist}

\title{
A function which subset all entries in a switchAnalyzeRlist.
}

\description{
This function allows the user to remove data from all entries in a switchAnalyzeRlist about isoforms that are no longer of interest. Note that it retain replicate isoforms information for all isoforms associated with genes containing isoforms in the subset (to enable correction for confounding factors when testing with isoformSwitchTestDEXSeq()).
}

\usage{
subsetSwitchAnalyzeRlist(
    switchAnalyzeRlist,
    subset
)
}


\arguments{
\item{switchAnalyzeRlist}{A \code{switchAnalyzeRlist} object.}

\item{subset}{logical expression indicating which rows in the \code{isoformFeatures} entry should be keep. The rest of the switchAnalyzeRlist is then reduced to only contain the matching information.}
}

\value{
A SwitchAnalyzeRlist only containing information about the isoforms (in their specific comparisons) indicated with TRUE in the .
}

\references{
Vitting-Seerup et al. The Landscape of Isoform Switches in Human Cancers. Mol. Cancer Res. (2017).
}
\author{
Kristoffer Vitting-Seerup
}

\seealso{
\code{\link{createSwitchAnalyzeRlist}}\cr
\code{\link{preFilter}}
}

\examples{
data("exampleSwitchList")

subsetSwitchAnalyzeRlist(
    switchAnalyzeRlist = exampleSwitchList,
    subset = exampleSwitchList$isoformFeatures$gene_overall_mean > 10
)
}


