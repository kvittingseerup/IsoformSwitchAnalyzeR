\name{prepareSalmonFileDataFrame}
\alias{prepareSalmonFileDataFrame}

\title{
Prepare data.frame needed run \link{importSalmonData}.
}

\description{
An easy to use wrapper for creating the "salmonFileDataFrame" data.frame needed to run \link{importSalmonData}.
}

\usage{
prepareSalmonFileDataFrame(
    ### Core arguments
    parentDir,

    ### Advanced arguments
    pattern='',
    invertPattern=FALSE,
    ignore.case=FALSE,
    quiet = FALSE
)
}

\arguments{
\item{parentDir}{ Parent directory where each quantified sample is in a sub-directory. The function will then look for files containing the (suffix) of the default files names for the quantification tools. The suffixes identified are 'abundance.tsv' for Kallisto, 'quant.sf' for Salmon, 'isoforms.results' for RSEM and 't_data.ctab' for StringTie. This is an alternative to \code{sampleVector} (aka only one of them should be used).}

\item{pattern}{ A character string containing a regular expression for which files to import (applied to full path). Default is "" corresponding to all. See base::grepl for more details.}

\item{invertPattern}{ A Logical. If TRUE only use files which do not match the \code{pattern} argument.}

\item{ignore.case}{ A logical. If TRUE case is ignored duing matching with the \code{pattern} argument. If FALSE the matching with the \code{pattern} argument is case sensitive.}

\item{quiet}{ A logic indicating whether to avoid printing progress messages (incl. progress bar). Default is FALSE}
}


\value{
The data.frame with 3 columns.
\itemize{
    \item{\code{Column 1: "files".} Contains the file each found in subdirectiories of the \code{parentDir} directory.}
    \item{\code{Column 2: "names".} The name of the subdirectory.}
    \item{\code{Column 3: "condition".} Set to NA as the function does not attemp to guess conditions. To use \link{importSalmonData} you will need to add these manually.}
}
}


\references{
Vitting-Seerup et al. The Landscape of Isoform Switches in Human Cancers. Mol. Cancer Res. (2017).
}

\author{
Kristoffer Vitting-Seerup
}

\seealso{
\code{\link{importSalmonData}}
}

\examples{
### Please note
# The way of importing files in the following example with
#   "system.file('pathToFile', package="IsoformSwitchAnalyzeR") is
#   specialized way of accessing the example data in the IsoformSwitchAnalyzeR package
#   and not something you need to do - just supply the string e.g.
#   parentDir = "individual_quantifications_in_subdir/" to the functions
#   path (e.g. "myAnnotation/isoformsQuantified.gtf") to the isoformExonAnnoation argument

### Prepare data.frame with quant file info
salmonDf <- prepareSalmonFileDataFrame(
    system.file("extdata/drosophila", package="IsoformSwitchAnalyzeR")
)

### Add conditions
salmonDf$condition <- c('wt','wt','ko','ko')

### Create switchAnalyzeRlist
aSwitchList <- importSalmonData(salmonDf)
}
