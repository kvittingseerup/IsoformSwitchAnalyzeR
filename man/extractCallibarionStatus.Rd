\name{extractCalibrationStatus}
\alias{extractCalibrationStatus}

\title{
Extract Whether P-values were Corrected.
}
\description{
This function allows the user to extract information about whether the p-values were calibrated as described in Ferguson et al or not.
}
\usage{
extractCalibrationStatus(
    switchAnalyzeRlist
)
}

\arguments{
\item{switchAnalyzeRlist}{
A \code{switchAnalyzeRlist} object which have been analyzed with isoformSwitchTest.
}
}

\value{
A data.frame where it for each comparison was indicated (via a logical) whether the p-value correction was performed or not.
}

\references{
\itemize{Manuscript in preperation.}
\itemize{Ferguson JP, Palejev D: P-value calibration fro multiple testing problems in genomics. Stat. Appl. Genet. Mol. Biol. 2014, 13:659-673.}
}

\author{
Kristoffer Vitting-Seerup
}

\seealso{
\code{\link{isoformSwitchTest}}
}
\examples{
# Load example data and prefilter
data("exampleSwitchList")
exampleSwitchList <- preFilter(exampleSwitchList)

# Perfom test
exampleSwitchListAnalyzed <- isoformSwitchTest(exampleSwitchList)

# extract whether p-value callibration was performed
extractCalibrationStatus(exampleSwitchListAnalyzed)
}
