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

\item{featureToImport}{A string controling the normalization of the feature. Supported options are 'counts' (default),'cpm','rpkm' and 'rpkmEffective'. See details.}

\item{showProgress}{ A logic indicating whether to make a progress bar (if TRUE) or not (if FALSE). Default is FALSE.}

\item{quiet}{ A logic indicating whether to avoid printing progress messages (incl. progress bar). Default is FALSE}
}

\details{
This function requires all data that should be imported is in a directory (as indicated by \code{parentDir}) where each quantified sample is in a seperate sub-directory.

Following feature can be extracted:
\itemize{
    \item{\code{counts}: The estimated number of reads assigned to a specific transcript.}
    \item{\code{cpm}: Counts Per Million. Simply \code{coints} normalized by library size in millions. Enables comparison of a transcript accross conditions.}
    \item{\code{rpkm}: Reads Per Kilobase transcript per Million. Simply \code{cpm} normalized by isoform lenght. Enables comparison of expression of different isoforms.}
    \item{\code{rpkmEffective}: \code{rpkm} calculated using the effective isoform legnth. Only supported by Kallisto and Salmon. Note that the effecive length correction also includes correction of all factors being modeled including the fragment length distribution and sequence-specific bias as well as GC-fragment bias (if used when running Kallisto and Salmon). Optimal for expression comparison since most biases will be taking into account.}
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

