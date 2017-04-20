\name{analyzeSignalP}
\alias{analyzeSignalP}

\title{
Import Result of SignalP Analysis
}
\description{
Allows for easy integration of the result of SignalP (external sequence analysis of signal peptides) in the IsoformSwitchAnalyzeR workflow. Please note that due to the 'removeNoncodinORFs' of the \code{analyzeCPAT} argument we recommend using \code{analyzeCPAT} before \code{analyzePFAM} and \code{analyzeSignalP} if you have predicted the ORFs with \code{analyzeORF}.
}
\usage{
analyzeSignalP(
    switchAnalyzeRlist,
    pathToSignalPresultFile,
    quiet=FALSE
)
}

\arguments{
  \item{switchAnalyzeRlist}{
A \code{switchAnalyzeRlist} object
}
  \item{pathToSignalPresultFile}{
A string indicating the full path to the summary SignalP result file. See \code{details} for suggestion of how to run and obtain the result of the SignalP tool.
}
\item{quiet}{ A logic indicating whether to avoid printing progress messages (incl. progress bar). Default is FALSE}
}

\details{
Notes for how to run the external tools:
If using the web-server (\url{http://www.cbs.dtu.dk/services/SignalP/}) SignalP should be run with the paramter "standard" under "Output format" and "No graphics" under "Graphics output". When using a stand-alone version SignalP should be run with the '-f summary' option. If using the web-server the results should be copy pasted (from the web-page) into an empty plain text document (notepat, sublimetext TextEdit or similar (not word)) and save that to a txt file. This file is then used as input to the function. If a stand-alone version was just supply the path to the summary result file
}

\value{
A collum called 'signal_peptide_identified' is added to \code{isoformFeatures} containing a binary indiaction (yes/no) of whether a transcript contains a signal peptide or not. Furthermore the data.frame 'signalPeptideAnalysis' is added to the \code{switchAnalyzeRlist} containing the details of the signal peptide analysis.

The data.frame added have one row pr isoform and contains 6 columns:
\itemize{
    \item{\code{isoform_id}: The name of the isoform analyzed. Matches the 'isoform_id' entry in the 'isoformFeatures' entry of the switchAnalyzeRlist}
    \item{\code{has_signal_peptide}: A text string indicating whether there is a signal peptide or not. Can be yes or no}
    \item{\code{network_used}: A text string indicating whether SignalP used the Neural Network (NN) optimized for proteins with trans-membrane sections (string='TM') or proteins without trans-membrane sections (string='noTM'). Per default, SignalP 4.1 uses the NN with TM as a preprocessor to determine whether to use TM or noTM in the final prediction (if 4 or more positions are predicted to be in a transmembrane state, TM is used, otherwise SignalP-noTM). Refrence: \url{http://www.cbs.dtu.dk/services/SignalP/instructions.php#method}}
    \item{\code{aa_removed}: A integer giving the number of amino acids removed when the signal peptide is cleaved off.}
    \item{\code{transcriptClevageAfter}: The transcript position of the last nucleotide in the isoform which is removed when the signal peptide is cleaved off.}
    \item{\code{genomicClevageAfter}: The genomic position of the last nucleotide in the isoform which is removed when the signal peptide is cleaved off.}
}
}

\references{
\itemize{
    \item{\code{This function} : Manuscript in preperation}
    \item{\code{SignalP} : Petersen et al. SignalP 4.0: discriminating signal peptides from transmembrane regions. Nature Methods, 8:785-786, 2011}
}
}

\author{
Kristoffer Vitting-Seerup
}
\seealso{
\code{\link{createSwitchAnalyzeRlist}}
\code{\link{extractSequence}}
\code{\link{analyzeSwitchConsequences}}
}
\examples{
### Load example data (matching the result files also store in IsoformSwitchAnalyzeR)
data("exampleSwitchListIntermediary")
exampleSwitchListIntermediary

### Add SignalP analysis
exampleSwitchListAnalyzed <- analyzeSignalP(
    switchAnalyzeRlist       = exampleSwitchListIntermediary,
    pathToSignalPresultFile = system.file("extdata/signalP_results.txt", package = "IsoformSwitchAnalyzeR")
    )

exampleSwitchListAnalyzed
}
