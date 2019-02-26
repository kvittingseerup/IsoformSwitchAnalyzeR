\name{analyzeNetSurfP2}
\alias{analyzeNetSurfP2}

\title{
Import Result of NetSurfP2 analysis
}
\description{
Allows for easy integration of the result of NetSurfP2 (performing external sequence analysis which include Intrinsically Disordered Regions (IDR)) in the IsoformSwitchAnalyzeR workflow. This function also supports using a sliding window to extract IDRs. Please note that due to the 'removeNoncodinORFs' option in \code{analyzeCPAT} and \code{analyzeCPC2} we recommend using analyzeCPC2/analyzeCPAT before using analyzeNetSurfP2, analyzePFAM and analyzeSignalP if you have predicted the ORFs with \code{analyzeORF}.
}
\usage{
analyzeNetSurfP2(
    switchAnalyzeRlist,
    pathToNetSurfP2resultFile,
    smoothingWindowSize = 5,
    probabilityCutoff = 0.5,
    minIdrSize = 30,
    showProgress = TRUE,
    quiet = FALSE
)
}

\arguments{
  \item{switchAnalyzeRlist}{
A \code{switchAnalyzeRlist} object
}
  \item{pathToNetSurfP2resultFile}{
A string indicating the full path to the NetSurfP-2 result file. Can be gziped.
}

\item{smoothingWindowSize}{
An integer indicating how large a sliding window should be used to calculate a smoothed (via mean) disordered probability score of a particular position in a peptide. This has as a smoothing effect which prevents IDRs from not being detected (or from being split into sub-IDRs) by a single residue with low probability. The tradeoff is worse accuracy of detecting the exact edges of the IDRs. To turn of smoothing simply set to 1. Default is 5 amino acids.
}

\item{probabilityCutoff}{
A doubble indicating the cutoff applied to the (smoothed) disordered probability score (see "smoothingWindowSize" argument above) for calling a residue as "disordered". The default, 30 amino acids, is an accepted standard for long IDRs.
}

\item{minIdrSize}{
An integer indicating how long a strech of disordered amino acid constiture the "region" part of the Intrinsically Disordered Region definition. The default, 30 amino acids, is an accepted standard for long IDRs.
}

\item{showProgress}{
A logic indicating whether to make a progress bar (if TRUE) or not (if FALSE). Default is TRUE.
}
\item{quiet}{ A logic indicating whether to avoid printing progress messages (incl. progress bar). Default is FALSE}
}

\details{
Intrinsically Disordered Regions (IDR) are regions of a protein which does not have a fixed three-dimentional structure (oposite protein domains). Such regions are thought to play important roles in all aspects of biolgiy (and when it goes wrong) through multiple different functional aspectcs - including facilitating protein interactions.

The NetSurfP web-server currently have no restriction on the number of sequences in the file uploaded so we sugget using the combined aa fasta file. See \link{extractSequence} for info on how to split the amino acid fasta files.

Notes for how to run the external tools: \cr
Use default paramters. If you want to use the webserver it is easily done as follows:.
1) Go to \url{http://www.cbs.dtu.dk/services/NetSurfP-2.0/}
2) Upload the amino acoid file (_AA) created with \code{extractSequence}.
3) Submit your job.
4) Wait till job is finished (if you submit your email you will recieve a notification).
5) In the top-right conor of the result site use the "Export All" botton to download the results as a CNV file.
6) Supply a string indicating the path to the downloaded cnv file directly to the "pathToNetSurfP2resultFile" argument.

IDR are only added to isoforms annotated as having an ORF even if other isoforms exists in the file. This means if you quantify the same isoform many times you can just run NetSurfP2 once on all isoforms and then supply the entire file to \code{analyzeNetSurfP2}.
}

\value{
A colum called 'idr_identified' is added to \code{isoformFeatures} containing a binary indication (yes/no) of whether a transcript contains any protein domains or not. Furthermore the data.frame 'idrAnalysis' is added to the \code{switchAnalyzeRlist} containing positional data of each IDR identified.

The data.frame added have one row per isoform and contains the columns:
\itemize{
    \item{\code{isoform_id}: The name of the isoform analyzed. Matches the 'isoform_id' entry in the 'isoformFeatures' entry of the switchAnalyzeRlist}
    \item{\code{orf_aa_start}: The start coordinate given as amino acid position (of the ORF).}
    \item{\code{orf_aa_end}: The end coordinate given as amino acid position (of the ORF).}
    \item{\code{transcriptStart}: The transcript coordinate of the start of the IDR.}
    \item{\code{transcriptEnd}: The transcript coordinate of the end of the IDR.}
    \item{\code{idrStarExon}: The exon index in which the start of the IDR is located.}
    \item{\code{idrEndExon}: The exon index in which the end of the IDR is located.}
    \item{\code{idrStartGenomic}: The genomic coordinat of the start of the IDR.}
    \item{\code{idrEndGenomic}:  The genomic coordinat of the end of the IDR.}
}
}

\references{
\itemize{
    \item{\code{This function} : Vitting-Seerup et al. The Landscape of Isoform Switches in Human Cancers. Mol. Cancer Res. (2017).}
    \item{\code{NetSurfP-2} : Klausen et al: NetSurfP-2.0: improved prediction of protein structural features by integrated deep learning. BioRxiv (2018).}
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
\code{\link{analyzeSwitchConsequences}}
}
\examples{
### Load example data (matching the result files also store in IsoformSwitchAnalyzeR)
data("exampleSwitchListIntermediary")
exampleSwitchListIntermediary

### Add NetSurfP-2 analysis
exampleSwitchListAnalyzed <- analyzeNetSurfP2(
    switchAnalyzeRlist = exampleSwitchListIntermediary,
    pathToNetSurfP2resultFile = system.file("extdata/netsurfp2_results.csv.gz", package = "IsoformSwitchAnalyzeR")
)

exampleSwitchListAnalyzed
}
