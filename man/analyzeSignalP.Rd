\name{analyzeSignalP}
\alias{analyzeSignalP}

\title{
Import Result of SignalP Analysis
}
\description{
Allows for easy integration of the result of SignalP (external sequence analysis of signal peptides) in the IsoformSwitchAnalyzeR workflow. Please note that due to the 'removeNoncodinORFs' option in \code{analyzeCPAT} and \code{analyzeCPC2} we recommend using analyzeCPC2/analyzeCPAT before using analyzeSignalP, analyzeNetSurfP2, analyzePFAM if you have predicted the ORFs with \code{analyzeORF}.
}
\usage{
analyzeSignalP(
    switchAnalyzeRlist,
    pathToSignalPresultFile,
    minSignalPeptideProbability = 0.5,
    quiet=FALSE
)
}

\arguments{
  \item{switchAnalyzeRlist}{
A \code{switchAnalyzeRlist} object
}
  \item{pathToSignalPresultFile}{
A string indicating the full path to the summary SignalP result file(s). If multiple result files were created (multiple web-server runs) just supply all the paths as a vector of strings. See \code{details} for suggestion of how to run and obtain the result of the SignalP tool.
}
\item{minSignalPeptideProbability}{ A numeric between 0 and 1 indicating the minimum probability for calling a signal peptide. Default is 0.5}

\item{quiet}{ A logic indicating whether to avoid printing progress messages (incl. progress bar). Default is FALSE}
}

\details{
A signal peptide is a short peptide sequence which indicate a protein is destined towards the secrotry pathway.

The SignalP web-server is less stringent than PFAM with regards to the number of sequences in the files uploaded so we sugget trying the combined fasta file first - and if that does not work try the files containg subsets. See \link{extractSequence} for info on how to split the amino acid fasta files.

Notes for how to run the external tools:
If using the web-server (\url{http://www.cbs.dtu.dk/services/SignalP/}) SignalP should be run with the paramter "Short output (no figures)" under "Output format" and one should sepect the appropriate "Organism group". When using a stand-alone version SignalP should be run with the '-f summary' option. If using the web-server the results can be downloaded using the "Downloads" botton in the top-right conor where the user should select "Prediction summary" and supply the path to the resulting file to the pathToSignalPresultFile argument. If a stand-alone version was just supply the path to the summary result file.


Please note that the \code{analyzeSignalP()} function will automatically only import the SignalP results from the isoforms stored in the switchAnalyzeRlist - even if many more are stored in the result file.
}

\value{
A collum called 'signal_peptide_identified' is added to \code{isoformFeatures} containing a binary indiaction (yes/no) of whether a transcript contains a signal peptide or not. Furthermore the data.frame 'signalPeptideAnalysis' is added to the \code{switchAnalyzeRlist} containing the details of the signal peptide analysis.

The data.frame added have one row pr isoform and contains 6 columns:
\itemize{
    \item{\code{isoform_id}: The name of the isoform analyzed. Matches the 'isoform_id' entry in the 'isoformFeatures' entry of the switchAnalyzeRlist}
    \item{\code{has_signal_peptide}: A text string indicating whether there is a signal peptide or not. Can be yes or no}
    \item{\code{network_used}: A text string indicating whether SignalP used the Neural Network (NN) optimized for proteins with trans-membrane sections (string='TM') or proteins without trans-membrane sections (string='noTM'). Per default, SignalP 4.1 uses the NN with TM as a preprocessor to determine whether to use TM or noTM in the final prediction (if 4 or more positions are predicted to be in a transmembrane state, TM is used, otherwise SignalP-noTM). Refrence: \url{http://www.cbs.dtu.dk/services/SignalP/instructions.php}}
    \item{\code{aa_removed}: A integer giving the number of amino acids removed when the signal peptide is cleaved off.}
    \item{\code{transcriptClevageAfter}: The transcript position of the last nucleotide in the isoform which is removed when the signal peptide is cleaved off.}
    \item{\code{genomicClevageAfter}: The genomic position of the last nucleotide in the isoform which is removed when the signal peptide is cleaved off.}
}
}

\references{
\itemize{
    \item{\code{This function} : Vitting-Seerup et al. The Landscape of Isoform Switches in Human Cancers. Mol. Cancer Res. (2017).}
    \item{\code{SignalP} : Almagro et al. SignalP 5.0 improves signal peptide predictions using deep neural networks. Nat. Biotechnol (2019).}
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
\code{\link{analyzeSwitchConsequences}}
}
\examples{
### Load example data
data("exampleSwitchListIntermediary")
exampleSwitchListIntermediary

### Add SignalP analysis
exampleSwitchListAnalyzed <- analyzeSignalP(
    switchAnalyzeRlist       = exampleSwitchListIntermediary,
    pathToSignalPresultFile = system.file(
        "extdata/signalP_results.txt",
        package = "IsoformSwitchAnalyzeR")
    )

exampleSwitchListAnalyzed
}
