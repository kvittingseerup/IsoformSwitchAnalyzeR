\name{subsetSwitchAnalyzeRlist}
\alias{subsetSwitchAnalyzeRlist}

\title{
A function which subset all enteries in a switchAnalyzeRlist.
}

\description{
This function allows the user to remove data from all entereis in a switchAnalyzeRlist about isoforms that are no longer of interest.
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
A SwitchAnalyzeRlist only containing information about the isoforms (in their specific comparisons) indeicated with TRUE in the .
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

