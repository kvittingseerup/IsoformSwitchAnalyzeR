\name{analyzeIUPred2A}
\alias{analyzeIUPred2A}

\title{
Import Result of IUPred2A or IUPred3 analysis
}
\description{
Allows for easy integration of the result of IUPred2A or IUPred3 (performing external sequence analysis of Intrinsically Disordered Regions (IDR) and Intrinsically Disordered Binding Regions (IDBR) ) in the IsoformSwitchAnalyzeR workflow. This function also supports using a sliding window to extract IDRs. Please note that due to the 'removeNoncodinORFs' option in \code{analyzeCPAT} and \code{analyzeCPC2} we recommend using analyzeCPC2/analyzeCPAT before using analyzeIUPred2A, analyzePFAM and analyzeSignalP if you have predicted the ORFs with \code{analyzeORF}.
}
\usage{
analyzeIUPred2A(
    switchAnalyzeRlist,
    pathToIUPred2AresultFile,
    smoothingWindowSize = 11,
    probabilityCutoff = 0.5,
    minIdrSize = 30,
    annotateBindingSites = TRUE,
    minIdrBindingSize = 15,
    minIdrBindingOverlapFrac = 0.8,
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
  \item{pathToIUPred2AresultFile}{
A string indicating the full path to the IUPred2A result file. If multiple result files were created (multiple web-server runs) just supply all the paths as a vector of strings. Can both be gziped or unpacked.
}

\item{smoothingWindowSize}{
An integer indicating how large a sliding window should be used to calculate a smoothed (via mean) disordered probability score of a particular position in a peptide. This has as a smoothing effect which prevents IDRs from not being detected (or from being split into sub-IDRs) by a single residue with low probability. The trade off is worse accuracy of detecting the exact edges of the IDRs. To turn of smoothing simply set to 1. Default is 5 amino acids.
}

\item{probabilityCutoff}{
A double indicating the cutoff applied to the (smoothed) disordered probability score (see "smoothingWindowSize" argument above) for calling a residue as "disordered". The default, 30 amino acids, is an accepted standard for long IDRs.
}

\item{minIdrSize}{
An integer indicating how long a stretch of disordered amino acid constitute the "region" part of the Intrinsically Disordered Region (IDR) definition. The default, 30 amino acids, is an accepted standard for long IDRs.
}

\item{annotateBindingSites}{
An logic indicating whether to also integrate the ANCHOR2 prediction of Intrinsically Disordered Binding Regions (IDBRs). See details for more info. Default is TRUE.
}

\item{minIdrBindingSize}{
An integer indicating how long a stretch of binding site the "region" part of the Intrinsically Disordered Binding Regions (IDBR) is defined as. Default is 15 AA.
}

\item{minIdrBindingOverlapFrac}{
An numeric indicating the min fraction of a predicted IDBR must also be within a IDR before the IDR is considered as a an IDR with a binding region. See details for more info. Default is 0.8 (aka 80\%).
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
Intrinsically Disordered Regions (IDR) are regions of a protein which does not have a fixed three-dimensional structure (opposite protein domains). Such regions are thought to play important roles in all aspects of biology (and when it goes wrong) through multiple different functional aspects. One such functional aspect is facilitating protein interactions via regions called Intrinsically Disordered Binding Regions (IDBR).

The IUPred2A webserver is somewhat strict with regards to the number of sequences in the files uploaded so we suggest multiple runs each with one of the the files contain subsets. See \link{extractSequence} for info on how to split the amino acid fasta files.

Notes for how to run the webserver: \cr
1) Go to \url{https://iupred2a.elte.hu}
2) Upload the amino avoid file (_AA) created with \code{extractSequence}.
3) Add your email (you will receive a notification when the job is done).
4) Use default parameters ("IUPred2 long disorder" as prediction type and "ANCHOR2" for context dependent prediction):
4) In the email you receive when the results are done use the link given after "The text file can be found here:" and save the result (right click on a blank space and use "save as" or use the keyboard shortcut "Ctrl/cmd + s" (or use wget etc to download in the first place))
5) Supply a string indicating the path to the downloaded file directly to the "pathToIUPred2AresultFile" argument. If multiple files are created (multiple web-server runs) just supply the path to all of them as a string.

IDR are only added to isoforms annotated as having an ORF even if other isoforms exists in the result file. This means if you quantify the same isoform many times (many different IsoformSwitchAnalyzeR workflows on the same organism) you can just run IUPred2A once on all isoforms and then supply the entire file to \code{pathToIUPred2AresultFile}.

Please note that the \code{analyzeIUPred2A()} function will automatically only import the IUPred2A results from the isoforms stored in the switchAnalyzeRlist - even if many more are stored in the result file.

Notes on Intrinsically Disordered Binding Regions (IDBR): As the prediction of IDR and IDBR are done by two different tools we require two things to annotate the Intrinsically Disordered Binding Regions (IDBR). Firstly the IDBR (predicted by ANCHOR2) must be a region of at least \code{minIdrBindingSize} amino acids (after the smoothing) - note this is different from the \code{minIdrSize} parameter. Secondly the fraction of the IDBR which overlaps the IDR predictions (done by IUPred2, again after smoothing) must be at least \code{minIdrBindingOverlapFrac}. When that is the case the IDR type will be annotated as "IDR_w_binding_region" instead of just "IDR". The current default parameters have not been rigorously tested and should be considered experimental.
}

\value{
A column called 'idr_identified' is added to \code{isoformFeatures} containing a binary indication (yes/no) of whether a transcript contains any IDR regions or not. Furthermore the data.frame 'idrAnalysis' is added to the \code{switchAnalyzeRlist} containing positional data of each IDR identified.

The data.frame added have one row per isoform and contains the columns:
\itemize{
    \item{\code{isoform_id}: The name of the isoform analyzed. Matches the 'isoform_id' entry in the 'isoformFeatures' entry of the switchAnalyzeRlist}
    \item{\code{orf_aa_start}: The start coordinate given as amino acid position (of the ORF).}
    \item{\code{orf_aa_end}: The end coordinate given as amino acid position (of the ORF).}
    \item{\code{idr_type}: A text string indicating the IDR type (one of 'IDR' or 'IDR_w_binding_region'.}
    \item{\code{nr_idr_binding_sites_overlapping}: (Only if \code{annotateBindingSites=TRUE}). An integer indicating the number of IDBRs predicted within the IDR.}
    \item{\code{max_fraction_of_idr_binding_sites_overlapping}: (Only if \code{annotateBindingSites=TRUE}). A fraction indicating the the largest overlap any IDBR had with the IDR.}
    \item{\code{transcriptStart}: The transcript coordinate of the start of the IDR.}
    \item{\code{transcriptEnd}: The transcript coordinate of the end of the IDR.}
    \item{\code{idrStarExon}: The exon index in which the start of the IDR is located.}
    \item{\code{idrEndExon}: The exon index in which the end of the IDR is located.}
    \item{\code{idrStartGenomic}: The genomic coordinate of the start of the IDR.}
    \item{\code{idrEndGenomic}:  The genomic coordinate of the end of the IDR.}
}
}

\references{
\itemize{
    \item{\code{This function} : Vitting-Seerup et al. The Landscape of Isoform Switches in Human Cancers. Mol. Cancer Res. (2017).}
    \item{\code{IUPred2A} : Meszaros et al. IUPred2A: Context-dependent prediction of protein disorder as a function of redox state and protein binding. Nucleic Acids Res (2018).}
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
\code{\link{analyzeNetSurfP2}}\cr
\code{\link{analyzeSwitchConsequences}}
}

\examples{
### Please note
# The way of importing files in the following example with
#   "system.file('pathToFile', package="IsoformSwitchAnalyzeR") is
#   specialized way of accessing the example data in the IsoformSwitchAnalyzeR package
#   and not something you need to do - just supply the string e.g.
#   "myAnnotation/upred2aResultFile.txt" to the functions

data("exampleSwitchListIntermediary")

upred2aResultFile <- system.file(
    "extdata/iupred2a_result.txt.gz",
    package = "IsoformSwitchAnalyzeR"
)

exampleSwitchListAnalyzed <- analyzeIUPred2A(
    switchAnalyzeRlist       = exampleSwitchListIntermediary,
    pathToIUPred2AresultFile = upred2aResultFile,
    showProgress = FALSE
)

}
