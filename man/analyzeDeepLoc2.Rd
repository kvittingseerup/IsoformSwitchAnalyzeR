\name{analyzeDeepLoc2}
\alias{analyzeDeepLoc2}

\title{
Import Result of DeepLoc2 Analysis
}
\description{
Allows for easy integration of the result of DeepLoc2 (external sequence analysis of signal peptides) in the IsoformSwitchAnalyzeR workflow. Please note that due to the 'removeNoncodinORFs' option in \code{analyzeCPAT} and \code{analyzeCPC2} we recommend using analyzeCPC2/analyzeCPAT before using analyzeDeepLoc2, analyzeSignalP, analyzeNetSurfP2, analyzePFAM if you have predicted the ORFs with \code{analyzeORF}.
}
\usage{
analyzeDeepLoc2(
    switchAnalyzeRlist,
    pathToDeepLoc2resultFile,
    enforceProbabilityCutoff = TRUE,
    probabilityCutoff = NULL,
    quiet = FALSE
)
}

\arguments{
  \item{switchAnalyzeRlist}{
A \code{switchAnalyzeRlist} object
}
  \item{pathToDeepLoc2resultFile}{
A string indicating the full path to the DeepLoc2 result file(s). If multiple result files were created (multiple web-server runs) just supply all the paths as a vector of strings. See \code{details} for suggestion of how to run and obtain the result of the DeepLoc2 tool.
}

\item{enforceProbabilityCutoff}{ A logic indicating whether the location specific probability cutoff developed for DeepLoc2 should be strictly enforced (see details). This filter works independently from \code{probabilityCutoff}. If FALSE the localization called internally by DeepLoc2 are used. Default is NULL.}

\item{probabilityCutoff}{ A numeric between 0 and 1 indicating the minimum probability for calling a sub-cellular localization. This filter works independently from \code{enforceProbabilityCutoff}. If NULL the localization called internally by DeepLoc2 are used. Default is NULL.}


\item{quiet}{ A logic indicating whether to avoid printing progress messages (incl. progress bar). Default is FALSE}
}

\details{
Subcellular localization are extremely important for biological functions. The performance of DeepLoc2 have been optimized by using subcellular location specific (aka class specific) probability cutoffs. However if no location is above these tresholds the location with the largest probability is assigned. The \code{enforceProbabilityCutoff} argument turns of the assignment of class if the location specific tresholds are not met.

The DeepLoc2 web-server is stringent with regards to the number of sequences in the files uploaded so we suggest using the files containing subsets. See \link{extractSequence} for info on how to split the amino acid fasta files.

Notes for how to run the external tools:
NA

Please note that the \code{analyzeDeepLoc2()} function will automatically only import the DeepLoc2 results from the isoforms stored in the switchAnalyzeRlist - even if many more are stored in the result file.
}

\value{
A column called 'sub_cell_location' is added to \code{isoformFeatures} containing the predicted subcellular localization. If multiple locations are predicted they are seperated by comma. Furthermore the data.frame 'subCellLocationAnalysis' is added to the \code{switchAnalyzeRlist} containing the details of the sub-cellular location analysis.

The data.frame added have one row pr isoform and contains 12 columns:
\itemize{
    \item{\code{isoform_id}: The name of the isoform analyzed. Matches the 'isoform_id' entry in the 'isoformFeatures' entry of the switchAnalyzeRlist}
    \item{\code{Localizations}: A text string indicating the subcellular localization. If multiple locations are predicted they are seperated by comma. This string is derived as described for the \code{probabilityCutoff} argument.}
    \item{\code{rest}: One column for each sub-cellular location containing the predicted probability of the isoform beining in that location.}
}
}

\references{
\itemize{
    \item{\code{This function} : Vitting-Seerup et al. The Landscape of Isoform Switches in Human Cancers. Mol. Cancer Res. (2017).}
}
}

\author{
Kristoffer Vitting-Seerup
}
\seealso{
\code{\link{createSwitchAnalyzeRlist}}\cr
\code{\link{extractSequence}}\cr
\code{\link{analyzePFAM}}\cr
\code{\link{analyzeNetSurfP2}}\cr
\code{\link{analyzeCPAT}}\cr
\code{\link{analyzeSignalP}}\cr
\code{\link{analyzeSwitchConsequences}}
}

\examples{
### Please note the way of importing files in the following example with
# "system.file('pathToFile', package="IsoformSwitchAnalyzeR") is
# specialized way of accessing the example data in the IsoformSwitchAnalyzeR package
# and not something you need to do - just supply the string e.g.
# "myAnnotation/predicted_annoation.txt" to the function.

### Load example data (matching the result files also store in IsoformSwitchAnalyzeR)
data("exampleSwitchListIntermediary")
exampleSwitchListIntermediary

### Add sub-cellular location analysis
exampleSwitchListAnalyzed <- analyzeDeepLoc2(
    switchAnalyzeRlist = exampleSwitchListIntermediary,
    pathToDeepLoc2resultFile = system.file("extdata/deeploc2.csv", package = "IsoformSwitchAnalyzeR"),
    quiet = FALSE
)


exampleSwitchListAnalyzed
}
