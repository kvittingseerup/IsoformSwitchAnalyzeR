\name{importIsoformExpression}
\alias{importIsoformExpression}

\title{
Import expression data from Kallisto, Salmon or RSEM into R.
}

\description{
A general-purpose import function which imports isoform expression data from Kallisto, Salmon or RSEM into R as an expression matrix. This function is meant to be used to import the data and afterwards a switchAnalyzeRlist can be created with \code{importRdata}.
}

\usage{
importIsoformExpression(
    parentDir,
    featureToImport='counts',
    showProgress=TRUE,
    quiet=FALSE
)
}


\arguments{
\item{parentDir}{ Parrent directory where each quantified sample is in a sub-directory.}

\item{featureToImport}{A string controling the normalization of the feature. Supported options are 'counts' and 'TxPM'. See details.}

\item{showProgress}{ A logic indicating whether to make a progress bar (if TRUE) or not (if FALSE). Default is FALSE.}

\item{quiet}{ A logic indicating whether to avoid printing progress messages (incl. progress bar). Default is FALSE}
}

\details{
This function requires all data that should be imported is in a directory (as indicated by \code{parentDir}) where each quantified sample is in a seperate sub-directory.

Following feature can be extracted:
\itemize{
    \item{\code{counts}: The estimated number of reads assigned to a specific transcript. The is what is prefrered for differential expression due to the extra information about the uncertainty that counts carry.}
    \item{\code{TxPM}: \code{Transcripts Per Million} measrung the relative abundance of a specic transcript. An equivilent to RPKM/FPKM except it has been adjusted for as all the biases being modeled by the tools used for the quantifictation including the fragment length distribution and sequence-specific bias as well as GC-fragment bias (this is specific to each tool and how it was run so you need to look up the specific tool). The TxPM is optimal for expression comparison of abundances since most biases will be taking into account. This is the value called TPM by RSEM, Kallisto and Salmon but here refered to as TxPM to avoid confusion with the commenly used Tags Per Million.}
}

For Kallisto the bias estimation is enabled by adding '--bias' to the function call. For Salmon the bias estimation is enabled by adding '--seqBias' and '--gcBias' to the function call. For RSEM the bias estimation is enabled by adding '--estimate-rspd' to the function call.
}

\value{
A \code{data.frame} where the first column contains the isoform_ids and the rest are the replicate expression (as defined by the \code{featureToImport} argument).
}


\references{
Manuscript in preperation.
}
\author{
Kristoffer Vitting-Seerup
}

\seealso{
\code{\link{importRdata}}\cr
\code{\link{createSwitchAnalyzeRlist}}\cr
\code{\link{preFilter}}
}

