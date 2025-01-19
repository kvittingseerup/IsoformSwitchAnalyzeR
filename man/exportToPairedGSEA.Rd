\name{exportToPairedGSEA}
\alias{exportToPairedGSEA}
\title{
Export Data for pairedGSEA Analysis
}
\description{
This function prepares input data for paired paired DGE/DGS analysis by extracting a count matrix and corresponding metadata (design matrix) from a \code{switchAnalyzeRlist} object. The results can be saved as an RDS file for downstream analysis.
}
\usage{
exportToPairedGSEA(
    switchAnalyzeRlist,
    pathToOutput = getwd(),
    writeToFile = FALSE,
    rdsFileName = "pairedGSEA_input.rds",
    overwrite = FALSE,
    quiet = FALSE
)
}

\arguments{
  \item{switchAnalyzeRlist}{
A \code{switchAnalyzeRlist} object}
  \item{pathToOutput}{
The directory where the RDS file should be saved if \code{writeToFile=TRUE}. Default is the current working directory.
}
  \item{writeToFile}{
A logical value indicating whether to save the results as an RDS file. Default is \code{FALSE}.
}
  \item{rdsFileName}{
The name of the RDS file to be saved when \code{writeToFile=TRUE}. Default is \code{"pairedGSEA_input.rds"}.
}
  \item{overwrite}{
A logical value indicating whether to overwrite an existing RDS file with the same name. Default is \code{FALSE}.
}
  \item{quiet}{
A logical value indicating whether to suppress progress messages. Default is \code{FALSE}.
}
}


\details{
This function extracts a count matrix and design matrix from the \code{switchAnalyzeRlist} and prepares them for pairedGSEA analysis, specifically paired_diff function. The count matrix is modified such that row names follow the format \code{gene_id:isoform_id}. If the \code{writeToFile} argument is set to \code{TRUE}, the result is saved as an RDS file in the specified directory.
}

\value{
A list containing the following elements:
\itemize{
    \item \code{count_matrix}: A modified count matrix where row names follow the format \code{gene_id:isoform_id}.
    \item \code{metadata}: The design matrix extracted from the \code{switchAnalyzeRlist}.
}
If \code{writeToFile=TRUE}, the list is also saved as an RDS file in the specified directory.
}

\examples{
### Example usage
# Load example data
data("exampleSwitchList")

# Export data for paired GSEA without saving to file
pairedGSEAList <- exportToPairedGSEA(
    switchAnalyzeRlist = exampleSwitchList,
    writeToFile = FALSE
)
}

\author{
Chunxu Han, Kristoffer Vitting-Seerup
}

\seealso{
\code{\link{switchAnalyzeRlist}}}
}