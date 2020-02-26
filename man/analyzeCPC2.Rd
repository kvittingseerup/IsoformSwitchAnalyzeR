\name{analyzeCPC2}
\alias{analyzeCPC2}

\title{
Import Result of External Sequence Analysis
}
\description{
Allows for easy integration of the result of CPC2 (external sequence analysis of coding potential) in the IsoformSwitchAnalyzeR workflow. Please note that due to the 'removeNoncodinORFs' option we recommend using \code{analyzeCPC2} before \code{analyzePFAM} and \code{analyzeSignalP} if you have predicted the ORFs with \code{analyzeORF}. This is an alternative to analyzing CPAT results with \code{analyzeCPAT}.
}
\usage{
analyzeCPC2(
    switchAnalyzeRlist,
    pathToCPC2resultFile,
    codingCutoff = 0.5,
    removeNoncodinORFs,
    quiet=FALSE
)
}

\arguments{
    \item{switchAnalyzeRlist}{:A \code{switchAnalyzeRlist} object}
    \item{pathToCPC2resultFile}{: A string indicating the full path to the CPC2 result file. See \code{details} for suggestion of how to run and obtain the result of the CPAT tool.}
    \item{codingCutoff}{: Numeric indicating the cutoff used by CPC2 for distinguishing between coding and non-coding transcripts. The cutoff appears to be species independent. Default is 0.5.}
    \item{removeNoncodinORFs}{: A logic indicating whether to remove ORF information from the isoforms which the CPC2 analysis classifies as non-coding. This can be particular useful if the isoform (and ORF) was predicted de-novo but is not recommended if ORFs was imported from a GTF file. This will affect all downstream analysis and plots as both analysis of domains and signal peptides requires that ORFs are annotated (e.g. analyzeSwitchConsequences will not consider the domains (potentially) found by Pfam if the ORF have been removed).}
    \item{quiet}{: A logic indicating whether to avoid printing progress messages (incl. progress bar). Default is FALSE}
}

\details{
Notes for how to run the external tools:
Use default parameters and if required select the most similar species. If the [webserver](\url{http://cpc2.cbi.pku.edu.cn/batch.php}) (batch submission) was used, download the tab-delimited result file (via the "Download the result" button). If a stand-alone version was just just supply the path to the result file.

Please note that the \code{analyzeCPC2()} function will automatically only import the CPC2 results from the isoforms stored in the switchAnalyzeRlist - even if many more are stored in the result file.
}

\value{
Two columns are added to \code{isoformFeatures}: 'codingPotentialValue' and 'codingPotential' containing the predicted coding potential values and a logic indicating whether the isoform is coding or not respectively (based on the supplied cutoff).
}

\references{
    \itemize{
        \item{\code{This function} : Vitting-Seerup et al. The Landscape of Isoform Switches in Human Cancers. Mol. Cancer Res. (2017).}
        \item{\code{CPC2} : Kang et al CPC2: a fast and accurate coding potential calculator based on sequence intrinsic features. Nucleic Acids Res. 2017}
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
\code{\link{analyzeSignalP}}\cr
\code{\link{analyzeCPAT}}\cr
\code{\link{analyzeSwitchConsequences}}
}

\examples{
### Load example data (matching the result files also store in IsoformSwitchAnalyzeR)
data("exampleSwitchListIntermediary")
exampleSwitchListIntermediary

### Add CPC2 analysis
exampleSwitchListAnalyzed <- analyzeCPC2(
    switchAnalyzeRlist   = exampleSwitchListIntermediary,
    pathToCPC2resultFile = system.file("extdata/cpc2_result.txt", package = "IsoformSwitchAnalyzeR"),
    removeNoncodinORFs   = TRUE   # because ORF was predicted de novo
)

exampleSwitchListAnalyzed
}
