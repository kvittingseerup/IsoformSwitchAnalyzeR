\name{isoformSwitchAnalysisCombined}
\alias{isoformSwitchAnalysisCombined}
\title{
Isoform Switch Analysis Workflow: Extract, Annotate and Visualize all Significant Isoform Switches
}
\description{
This high-level function uses a pre-existing switchAnalyzeRlist as input. Then isoform switches are identified, annotated with ORF and intron retention. Then functional consequences are identified and isoform switch analysis plots are generated for the top \code{n} isoform switches. Lastly a plot summarizing the global effect of isoform switches with functional consequences is generated. If external analysis of protein domains (Pfam), coding potential (CPAT) or signal peptides (SignalP) should be incorporated please use the combination of \code{isoformSwitchAnalysisPart1} and \code{isoformSwitchAnalysisPart2} instead.
}
\usage{
isoformSwitchAnalysisCombined(
    ### Core arguments
    switchAnalyzeRlist,

    ### Annotation arguments
    genomeObject = NULL,
    pathToGTF = NULL,

    ### Analysis and output arguments
    n = Inf,
    consequencesToAnalyze = c('intron_retention', 'ORF_seq_similarity', 'NMD_status'),
    pathToOutput = getwd(),
    fileType = 'pdf',
    outputPlots = TRUE,

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

\item{n}{
The number of top genes (after filtering and sorted according to \code{sortByQvals}) that should be saved to each sub-folder indicated by \code{splitConditions, splitFunctionalConsequences}. Use Inf to create all. Default is Inf (all).
}

  \item{consequencesToAnalyze}{
    A vector of strings indicating what type of functional consequences to analyze. Do note that there is bound to be some differences between transcripts (else there would be identical). See details in \link{analyzeSwitchConsequences} for full list of usable strings and their meaning. Default is c('intron_retention','coding_potential','ORF_seq_similarity','NMD_status','domains_identified','signal_peptide_identified') (corresponding to analyze: intron retention, CPAT result, ORF AA sequence similarity, NMD status, PFAM domains annotated and signal peptides annotated by Pfam).
}

  \item{pathToOutput}{
A path to the folder in which the plots should be made. Default is working directory ( getwd() ).
}
  \item{fileType}{
A string indicating which file type is generated. Available options are \'pdf\' and \'png\'. Default is pdf.
}

\item{outputPlots}{
A logic indicating whether all isoform switches as well as the summary of functional consequences should be saved in the directory specified by \code{pathToOutput} argument. Default is TRUE.
}
\item{quiet}{ A logic indicating whether to avoid printing progress messages (incl. progress bar). Default is FALSE}

}

\details{
This function performs the full Isoform Analysis Workflow by
\enumerate{
    \item Remove non-expressed isoforms and single-isoform genes (see \link{preFilter})
    \item predict swithces (only if switches is not already annotated, see \link{isoformSwitchTestDEXSeq})
    \item Analyzing the isoforms for open reading frames (ORFs, see \link{analyzeORF})
    \item Output fasta files containing the nucleotide and amino acid sequences which enables external sequence analysis with CPAT, Pfam and SignalP (see \link{extractSequence})
    \item Predict functional consequences of switching (see \link{analyzeSwitchConsequences})
    \item Ouput Isoform Switch Analysis plots for all genes with a signicant switch (see \link{switchPlot})
    \item Ouput a visualization of general consequences of isoform switches.
}
}

\value{
This function outputs:
\enumerate{
    \item The supplied \code{switchAnalyzeRlist} now annotated with all the analysis described above
    \item One folder per comparison of condition containing the isoform switch analysis plot of all significant isoforms.
    \item A plot summarizing the overall consequences off all the isoform switches.
}
}

\references{
Vitting-Seerup et al. The Landscape of Isoform Switches in Human Cancers. Mol. Cancer Res. (2017).
}
\author{
Kristoffer Vitting-Seerup
}
\seealso{
\code{\link{isoformSwitchAnalysisPart1}}\cr
\code{\link{isoformSwitchAnalysisPart2}}\cr
\code{\link{preFilter}}\cr
\code{\link{isoformSwitchTestDEXSeq}}\cr
\code{\link{isoformSwitchTestSatuRn}}\cr
\code{\link{analyzeORF}}\cr
\code{\link{extractSwitchSummary}}\cr
\code{\link{analyzeSwitchConsequences}}\cr
\code{\link{switchPlotTopSwitches}}
}
\examples{
data("exampleSwitchList")
exampleSwitchList

library(BSgenome.Hsapiens.UCSC.hg19)
exampleSwitchListAnalyzed <- isoformSwitchAnalysisCombined(
    switchAnalyzeRlist=exampleSwitchList,
    outputPlots = FALSE      # keeps the function from outputting the Isoform Switch AnalyzeR Plots from this example
)

exampleSwitchListAnalyzed
}
