\name{isoformSwitchAnalysisPart1}
\alias{isoformSwitchAnalysisPart1}
\title{
Isoform Switch Analysis Workflow Part 1: Extract Isoform Switches and Their Sequences
}
\description{
This high-level function takes either a CuffSet object or a pre-existing switchAnalyzeRlist as input. If the input is a CuffSet object a switchAnalyzeRlist is created or else the function uses the provieded switchAnalyzeRlist. Then isoform switches are predicted (unless \code{switchTestMethod='none'}) and ORF are predicted if not already annotated. Lastly the function extracts the nucleotide sequence and the ORF AA sequences of the isoforms involved in isoform switches. Thes sequences are both saved to external files and added to the switchAnalyzeRlist to enable external and internal sequence analysis respectively.

This function is meant to be used as part 1 of the isoform switch analysis workflow, which can be followed by the second step via \code{isoformSwitchAnalysisPart2}.
}

\usage{
isoformSwitchAnalysisPart1(
    switchAnalyzeRlist,
    alpha = 0.05,
    dIFcutoff = 0.1,
    switchTestMethod='DEXSeq',
    orfMethod = "longest",
    genomeObject = NULL,
    cds = NULL,
    pathToOutput = getwd(),
    outputSequences = TRUE,
    prepareForWebServers = FALSE,
    overwriteORF=FALSE,
    quiet=FALSE
)
}

\arguments{
    \item{switchAnalyzeRlist}{
A \code{switchAnalyzeRlist}.
}
   \item{alpha}{ The cutoff which the (calibrated) fdr correct p-values must be smaller than for calling significant switches. Defualit is 0.05.
}
    \item{dIFcutoff}{ The cutoff which the changes in (absolute) isoform usage must be larger than before an isoform is considered switching. This cutoff can remove cases where isoforms with (very) low dIF values are deemed significant and thereby included in the downstream analysis. This cutoff is analogous to having a cutoff on log2 fold change in a normal differential expression analysis of genes to ensure the genes have a certain effect size. Default is 0.1 (10\%).}

    \item{switchTestMethod}{ A sting indicating which statistical method should be used for testing differential isoform usage. The following options are avilable:
        \itemize{
            \item {\code{'DEXSeq'} : Uses DEXSeq to perform the statiscal test. See \link{isoformSwitchTestDEXSeq}. Default}
            \item {\code{'DRIMSeq'} : Uses the DRIMSeq package to perform the statiscal test. See \link{isoformSwitchTestDRIMSeq}.}
            \item {\code{'none'} : No statistical test is performed. Should only be used if a test have already been performed and should not be overwritten (e.g when importing cuffdiff data).}
        }
    }

  \item{orfMethod}{
  A string indicating which of the 4 ORF identification methods should be used. The methods are:
  \itemize{
    \item {\code{longest} : Identifies the longest ORF in the transcript. This approach is similar to what the CPAT tool uses in it's analysis of coding potential}
    \item {\code{longestAnnotated} : Identifies the longest ORF downstream of an annoated translation start site (supplied via the \code{cds} argument)}
    \item {\code{mostUpstreamAnnoated} : Identifies the ORF downstream of the most unstream overlapping annoated translation start site (supplied via the \code{cds} argument)}
  }
  Default is \code{longest}.
}
  \item{genomeObject}{
A \code{BSgenome} object (for example Hsapiens for Homo sapiens).
}
  \item{pathToOutput}{
A path to the folder in which the plots should be made. Default is working directory ( getwd() ).
}
  \item{cds}{
A CDSSet object containing annoated coding regions, see ?CDSSet and ?getCDS for more information. Only necessary if \'orfType\' arguments is \'longestAnnotated\' or \'mostUpstreamAnnoated\'.
}

\item{overwriteORF}{
A logical indicating whether to overwrite the ORF analysis already stored in the supplied switchAnalyzeRlist. Default is FALSE.
}
  \item{outputSequences}{
A logical indicating whether transcript nucleotide and amino acid sequences should be outputted to \code{pathToOutput}. Default is TRUE.
}
\item{prepareForWebServers}{
A logical indicating whether the amino acid fasta files saved (if \code{outputSequences=TRUE}) should be prepared for the online web-servises currently supported (as they have some limitations on what can submitted). See details. Default is FALSE (for backward compatability).
}

\item{quiet}{ A logical indicating whether to avoid printing progress messages (incl. progress bar). Default is FALSE}
}


\details{
This function performs the first part of a Isoform Analysis Workflow by
\enumerate{
    \item Remove non-expressed isoforms and single-isoform genes (see \link{preFilter})
    \item Predict isoform switches unless \code{switchTestMethod} is set to 'none'.
    \item If no ORFs are annotated the isoforms are analyzed for open reading frames (ORFs, see \link{analyzeORF})
    \item The isoform nucleotide and ORF amino acid sequences are extracted and saved to fasta files as well as added to the switchAnalyzeRlist enabling external sequence analysis with CPAT, Pfam and SignalP (see vignette for more info).
}

if \code{prepareForWebServers=TRUE} both the "filterAALength" and "alsoSplitFastaFile" will be enabled in the \code{extractSequence} function.
}

\value{
This function have two outputs. It returns a \code{switchAnalyzeRlist} object where information about the isoform switch test, ORF prediction and nt and aa sequences have been added. Secondly (if \code{outputSequences} is TRUE) the nucleotide and amino acid sequence of transcripts involved in switches are also save as fasta files enabling external sequence analysis.
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
data("exampleSwitchList")
exampleSwitchList

exampleSwitchList <- isoformSwitchAnalysisPart1(
    switchAnalyzeRlist=exampleSwitchList,
    dIFcutoff = 0.4,        # Set high for short runtime in example data
    outputSequences = FALSE # keeps the function from outputting the fasta files from this example
)

exampleSwitchList
}
