\name{analyzeDeepTMHMM}
\alias{analyzeDeepTMHMM}

\title{
Import Result of a DeepTMHMM analysis
}
\description{
Allows for easy integration of the result of DeepTMHMM (performing external sequence analysis of isoform topology) in the IsoformSwitchAnalyzeR workflow. Please note that due to the 'removeNoncodinORFs' option in \code{analyzeCPAT} and \code{analyzeCPC2} we recommend using analyzeCPC2/analyzeCPAT before using analyzeTopcons2, analyzeDeepTMHMM, analyzeNetSurfP2, analyzePFAM and analyzeSignalP if you have predicted the ORFs with \code{analyzeORF}.
}
\usage{
analyzeDeepTMHMM(
    switchAnalyzeRlist,
    pathToDeepTMHMMresultFile,
    ignoreAfterBar = TRUE,
    ignoreAfterSpace = TRUE,
    ignoreAfterPeriod = FALSE,
    showProgress = TRUE,
    quiet = FALSE
)
}

\arguments{
  \item{switchAnalyzeRlist}{
A \code{switchAnalyzeRlist} object
}
  \item{pathToDeepTMHMMresultFile}{
A string indicating the full path to the DeepTMHMM result file. Can be gziped. If multiple result files were created (multiple web-server runs) just supply all the paths as a vector of strings.
}
\item{ignoreAfterBar}{A logic indicating whether to subset the isoform ids by ignoring everything after the first bar ("|"). Useful for analysis of GENCODE data. Default is TRUE.}
\item{ignoreAfterSpace}{A logic indicating whether to subset the isoform ids by ignoring everything after the first space (" "). Useful for analysis of gffutils generated GTF files. Default is TRUE.}
\item{ignoreAfterPeriod}{ A logic indicating whether to subset the gene/isoform is by ignoring everything after the first period ("."). Should be used with care. Default is FALSE.}

\item{showProgress}{
A logic indicating whether to make a progress bar (if TRUE) or not (if FALSE). Default is TRUE.
}
\item{quiet}{ A logic indicating whether to avoid printing progress messages (incl. progress bar). Default is FALSE}
}

\details{
The topological structure of a protein is the predition/annoation of which parts of a membrane associated protein are on the inside, within and on the outside of the cell membrane. This is very important knowleadge when designing drugs or trying to understand intercellular communication.

DeepTMHMM can be run from from \url{https://biolib.com/DTU/DeepTMHMM} and afterwards all files can be downloaded as a "gff3 format" file can be used as input to this function.
}

\value{
A data.frame 'topologyAnalysis' is added to the \code{switchAnalyzeRlist} containing the type of region(s) as well as positional data of that region for each isoform.

The data.frame added have one row per topological region of an isoform and contains the columns:
\itemize{
    \item{\code{isoform_id}: The name of the isoform analyzed. Matches the 'isoform_id' entry in the 'isoformFeatures' entry of the switchAnalyzeRlist}
    \item{\code{region_type}: A text string indicating the location of the region compared to the membrane.}
    \item{\code{orf_aa_start}: The start coordinate given as amino acid position (of the ORF).}
    \item{\code{orf_aa_end}: The end coordinate given as amino acid position (of the ORF).}
    \item{\code{transcriptStart}: The transcript coordinate of the start of the IDR.}
    \item{\code{transcriptEnd}: The transcript coordinate of the end of the IDR.}
    \item{\code{regionStarExon}: The exon index in which the start of the IDR is located.}
    \item{\code{regionEndExon}: The exon index in which the end of the IDR is located.}
    \item{\code{regionStartGenomic}: The genomic coordinate of the start of the IDR.}
    \item{\code{regionEndGenomic}:  The genomic coordinate of the end of the IDR.}
}
}

\references{
\itemize{
    \item{\code{This function} : Vitting-Seerup et al. The Landscape of Isoform Switches in Human Cancers. Mol. Cancer Res. (2017).}
    \item{\code{DeepTMHMM} : Hallgren et al: In prep.}
}
}

\author{
Kristoffer Vitting-Seerup
}
\seealso{
\code{\link{createSwitchAnalyzeRlist}}\cr
\code{\link{extractSequence}}\cr
\code{\link{analyzeCPAT}}\cr
\code{\link{analyzeSignalP}}\cr
\code{\link{analyzePFAM}}\cr
\code{\link{analyzeIUPred2A}}\cr
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

### Add toplogical analysis
exampleSwitchListAnalyzed <- analyzeDeepTMHMM(
    switchAnalyzeRlist   = exampleSwitchListIntermediary,
    pathToDeepTMHMMresultFile = system.file("extdata/DeepTMHMM.gff3", package = "IsoformSwitchAnalyzeR"),
    showProgress=FALSE
)

exampleSwitchListAnalyzed
}
