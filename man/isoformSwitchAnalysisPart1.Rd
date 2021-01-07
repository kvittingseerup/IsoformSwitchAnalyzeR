\name{isoformSwitchAnalysisPart1}
\alias{isoformSwitchAnalysisPart1}
\title{
Isoform Switch Analysis Workflow Part 1: Extract Isoform Switches and Their Bio-sequences
}
\description{
This high-level function takes a pre-existing switchAnalyzeRlist as input (see \code{importRdata}). Then part 1 of the workflow is performed. Specifically it is filtered to remove low expression, isoform switches are identified via the statistical methods (unless \code{switchTestMethod='none'}) and ORF are predicted if not already annotated. Lastly the function extracts the nucleotide sequence and the ORF AA sequences of the isoforms involved in isoform switches. To enable external and internal sequence analysis these sequences are both saved to the computer (as fasta files) and added to the switchAnalyzeRlist.

This function is meant to be used as part 1 of the isoform switch analysis workflow, which can be followed by the second step via \code{isoformSwitchAnalysisPart2}.
}

\usage{
isoformSwitchAnalysisPart1(
    ### Core arguments
    switchAnalyzeRlist,

    ### Annotation arguments
    genomeObject = NULL,
    pathToGTF = NULL,

    ### Output arguments
    prepareForWebServers,
    pathToOutput = getwd(),
    outputSequences = TRUE,

    ### Other arguments
    quiet = FALSE
)
}

\arguments{
    \item{switchAnalyzeRlist}{
A \code{switchAnalyzeRlist}.
}

  \item{genomeObject}{
A \code{BSgenome} object (for example Hsapiens for Homo sapiens).
}
\item{pathToGTF}{
A string indicating the full path to the (gziped or unpacked) GTF file which contains the the known annotation (aka from a official source) which was used to guided the transcript assembly (isoform deconvolution).
}


\item{prepareForWebServers}{
A logical indicating whether the amino acid fasta files saved (if \code{outputSequences=TRUE}) should be prepared for the online web-services currently supported (as they have some limitations on what can submitted). See details. Default is FALSE (for backward compatibility).
}

  \item{pathToOutput}{
A path to the folder in which the plots should be made. Default is working directory ( getwd() ).
}

  \item{outputSequences}{
A logical indicating whether transcript nucleotide and amino acid sequences should be outputted to \code{pathToOutput}. Default is TRUE.
}


\item{quiet}{ A logical indicating whether to avoid printing progress messages (incl. progress bar). Default is FALSE}
}


\details{
This function performs the first part of a Isoform Analysis Workflow by
\enumerate{
    \item Remove low-expressed isoforms and single-isoform genes (see \link{preFilter})
    \item Identify isoform switches.
    \item If no ORFs are annotated the isoforms are analyzed for open reading frames (ORFs, see \link{analyzeORF})
    \item The isoform nucleotide and ORF amino acid sequences are extracted and saved to fasta files as well as added to the switchAnalyzeRlist enabling external sequence analysis with CPAT, Pfam and SignalP (see vignette for more info).
}

if \code{prepareForWebServers=TRUE} both the "removeLongAAseq" and "alsoSplitFastaFile" will be enabled in the \code{extractSequence} function.
}

\value{
This function have two outputs. It returns a \code{switchAnalyzeRlist} object where information about the isoform switch test, ORF prediction and nt and aa sequences have been added. Secondly (if \code{outputSequences=TRUE}) the nucleotide and amino acid sequence of transcripts involved in switches are also save as fasta files enabling external sequence analysis.
}

\references{
Vitting-Seerup et al. The Landscape of Isoform Switches in Human Cancers. Mol. Cancer Res. (2017).
}
\author{
Kristoffer Vitting-Seerup
}
\seealso{
\code{\link{preFilter}}\cr
\code{\link{isoformSwitchTestDEXSeq}}\cr
\code{\link{isoformSwitchTestDRIMSeq}}\cr
\code{\link{analyzeORF}}\cr
\code{\link{extractSequence}}

}
\examples{
### load example data
data("exampleSwitchList")

### Subset for quick runtime
exampleSwitchList <- subsetSwitchAnalyzeRlist(
    exampleSwitchList,
    abs(exampleSwitchList$isoformFeatures$dIF) > 0.4
)

### Show summary
summary(exampleSwitchList)

### Run Part 1
exampleSwitchList <- isoformSwitchAnalysisPart1(
    switchAnalyzeRlist=exampleSwitchList,
    prepareForWebServers = FALSE,
    outputSequences = FALSE # keeps the function from outputting the fasta files from this example
)

### Show summary
summary(exampleSwitchList)
}
