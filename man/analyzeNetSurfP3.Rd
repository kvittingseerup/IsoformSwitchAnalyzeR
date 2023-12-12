\name{analyzeNetSurfP3}
\alias{analyzeNetSurfP3}

\title{
Import Result of NetSurfP3 analysis
}
\description{
Allows for easy integration of the result of NetSurfP3 (performing external sequence analysis which include Intrinsically Disordered Regions (IDR)) in the IsoformSwitchAnalyzeR workflow. This function also supports using a sliding window to extract IDRs. Please note that due to the 'removeNoncodinORFs' option in \code{analyzeCPAT} and \code{analyzeCPC2} we recommend using analyzeCPC2/analyzeCPAT before using analyzeNetSurfP3, analyzePFAM and analyzeSignalP if you have predicted the ORFs with \code{analyzeORF}.
}
\usage{
analyzeNetSurfP3(
    switchAnalyzeRlist,
    pathToNetSurfP3resultFile,
    smoothingWindowSize = 5,
    probabilityCutoff = 0.5,
    minIdrSize = 30,
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
  \item{pathToNetSurfP3resultFile}{
A string indicating the full path to the NetSurfP-3 result file. Can be gziped. If multiple result files were created (multiple web-server runs) just supply all the paths as a vector of strings.
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

\item{ignoreAfterBar}{A logic indicating whether to subset the isoform ids by ignoring everything after the first bar ("|"). Useful for analysis of GENCODE data. Default is TRUE.}

\item{ignoreAfterSpace}{A logic indicating whether to subset the isoform ids by ignoring everything after the first space (" "). Useful for analysis of gffutils generated GTF files. Default is TRUE.}

\item{ignoreAfterPeriod}{ A logic indicating whether to subset the gene/isoform is by ignoring everything after the first period ("."). Should be used with care. Default is FALSE.}

\item{showProgress}{
A logic indicating whether to make a progress bar (if TRUE) or not (if FALSE). Default is TRUE.
}
\item{quiet}{ A logic indicating whether to avoid printing progress messages (incl. progress bar). Default is FALSE}
}

\details{
Intrinsically Disordered Regions (IDR) are regions of a protein which does not have a fixed three-dimensional structure (opposite protein domains). Such regions are thought to play important roles in all aspects of biology (and when it goes wrong) through multiple different functional aspects - including facilitating protein interactions.

The NetSurfP-3.0 web-server currently have a restriction of max 10,000 sequences in the file uploaded. If you have more than that (one for each isoform in \code{summary( switchAnalyzeRlist)}) we recommend multiple runs each with one of the the files containing subsets - else you can just run the combined fasta file. See \link{extractSequence} for info on how to split the amino acid fasta files.

Notes for how to run the external tools: \cr
Use default parameters. If you want to use the webserver it is easily done as follows:
1) Go to \url{https://services.healthtech.dtu.dk/services/NetSurfP-3.0/}
2) Upload the amino avoid file (_AA) created with \code{extractSequence}.
3) Submit your job.
4) Wait till job is finished (if you submit your email you will receive a notification).
5) In the top-right corner of the result site use the "Export All" bottom to download the results as a CSV file.
6) Supply a string indicating the path to the downloaded cnv file directly to the "pathToNetSurfP3resultFile" argument.

If you runNetSurfP-3.0 locally just use the "--csv" argument and provide the resulting csv file to the \code{pathToNetSurfP3resultFile} argument.


IDR are only added to isoforms annotated as having an ORF even if other isoforms exists in the file. This means if you quantify the same isoform many times you can just run NetSurfP3 once on all isoforms and then supply the entire file to \code{analyzeNetSurfP3}.

Please note that the \code{analyzeNetSurfP3()} function will automatically only import the NetSurfP-3.0 results from the isoforms stored in the switchAnalyzeRlist - even if many more are stored in the result file.
}

\value{
A column called 'idr_identified' is added to \code{isoformFeatures} containing a binary indication (yes/no) of whether a transcript contains any protein domains or not. Furthermore the data.frame 'idrAnalysis' is added to the \code{switchAnalyzeRlist} containing positional data of each IDR identified.

The data.frame added have one row per isoform and contains the columns:
\itemize{
    \item{\code{isoform_id}: The name of the isoform analyzed. Matches the 'isoform_id' entry in the 'isoformFeatures' entry of the switchAnalyzeRlist}
    \item{\code{orf_aa_start}: The start coordinate given as amino acid position (of the ORF).}
    \item{\code{orf_aa_end}: The end coordinate given as amino acid position (of the ORF).}
    \item{\code{idr_type}: A text string indicating the IDR type (one of 'IDR' or 'IDR_w_binding_region'.}
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
    \item{\code{NetSurfP-3} : Høie et al: NetSurfP-3.0: accurate and fast prediction of protein structural features by protein language models and deep learning. Nucleic Acids Research (2022).}
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
