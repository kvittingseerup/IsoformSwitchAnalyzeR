\name{isoformSwitchAnalysisPart2}
\alias{isoformSwitchAnalysisPart2}
\title{
Isoform Switch Analysis Workflow Part 2: Plot All Isoform Switches and Their Annotation
}
\description{
This high-level function adds the results of the external sequence analysis supplied (if any), then proceeds to analyze alternative splicing. Then functional consequences of the isoform switches are identified and isoform switch analysis plots are created for the top \code{n} isoform switches. Lastly a plot summarizing the functional consequences is created. This function is meant to be used after \link{isoformSwitchAnalysisPart1} have been used.
}
\usage{
isoformSwitchAnalysisPart2(
    ### Core arguments
    switchAnalyzeRlist,

    ### External annotation arguments
    codingCutoff = NULL,
    removeNoncodinORFs,
    pathToCPATresultFile = NULL,
    pathToCPC2resultFile = NULL,
    pathToPFAMresultFile = NULL,
    pathToIUPred2AresultFile = NULL,
    pathToNetSurfP2resultFile = NULL,
    pathToSignalPresultFile = NULL,

    ### Analysis and output arguments
    n = Inf,
    consequencesToAnalyze = c(
        'intron_retention',
        'coding_potential',
        'ORF_seq_similarity',
        'NMD_status',
        'domains_identified',
        'IDR_identified',
        'IDR_type',
        'signal_peptide_identified'
    ),
    pathToOutput = getwd(),
    fileType = 'pdf',
    outputPlots = TRUE,

    ### Other arguments
    quiet = FALSE
)
}

\arguments{
  \item{switchAnalyzeRlist}{
The \code{switchAnalyzeRlist} object as produced by \link{isoformSwitchAnalysisPart1}
}

  \item{codingCutoff}{
Numeric indicating the cutoff used by CPAT/CPC2 for distinguishing between coding and non-coding transcripts.
\enumerate{
    \item For CPAT: The cutoff is dependent on species analyzed. Our analysis suggest that the optimal cutoff for overlapping coding and noncoding isoforms are 0.725 for human and 0.721 for mouse - HOWEVER the suggested cutoffs from the CPAT article (see references) derived by comparing known genes to random non-coding regions of the genome is 0.364 for human and 0.44 for mouse. No default is used.
    \item For CPC2: The cutoff suggested is 0.5 for all species, and this cutoff will be used if nothing is specified by the user
}
}

 \item{removeNoncodinORFs}{
A logic indicating whether to remove ORF information from the isoforms which the CPAT analysis classifies as non-coding. This can be particular useful if the isoform (and ORF) was predicted de-novo but is not recommended if ORFs was imported from a GTF file. This will affect all downstream analysis and plots as both analysis of domains and signal peptides requires that ORFs are annotated (e.g. analyzeSwitchConsequences will not consider the domains (potentially) found by Pfam if the ORF have been removed).
}

  \item{pathToCPATresultFile}{
Path to the CPAT result file. If the webserver is used please download the tab-delimited file from the bottom of the result page and give that as input, else simply supply the result file. See \link{analyzeCPAT} for details.
}
\item{pathToCPC2resultFile}{
Path to the CPC2 result file. If the webserver is used please download the tab-delimited file from the bottom of the result page and give that as input, else simply supply the result file. See \link{analyzeCPC2} for details.
}

\item{pathToPFAMresultFile}{
  A string indicating the full path to the Pfam result file(s). If multiple result files were created (multiple web-server runs) just supply all the paths as a vector of strings. If the webserver is used you need to copy paste the result part of the mail you get into a empty plain text document (notepad, sublimetext TextEdit or similar (aka not word)) and save that. See \link{analyzePFAM} for details.
}


  \item{pathToIUPred2AresultFile}{
  A string indicating the full path to the NetSurfP-2 result csv file. See \link{analyzeIUPred2A} for details.
}

  \item{pathToNetSurfP2resultFile}{
  A string indicating the full path to the NetSurfP-2 result csv file. See \link{analyzeNetSurfP2} for details.
}

  \item{pathToSignalPresultFile}{
  A string indicating the full path to the SignalP result file(s). If multiple result files were created (multiple web-server runs) just supply all the paths as a vector of strings. If using the web-server the results should be copy pasted into a empty plain text document (notepad, sublimetext TextEdit or similar (aka not word)) and save that. See \link{analyzeSignalP} for details.
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
This function performs the second part of a Isoform Analysis Workflow by:
\enumerate{
    \item Adding external sequence analysis (see \link{analyzeCPAT}, \link{analyzeCPC2}, \link{analyzePFAM} and \link{analyzeSignalP})
    \item Predict functional consequences of switching (see \link{analyzeSwitchConsequences})
    \item Output Isoform Switch Consequence plots for all genes where there is a significant isoform switch (see \link{switchPlot})
    \item Output a visualization of general consequences of isoform switches.
}
}

\value{
This function
\enumerate{
    \item Returns the supplied \code{switchAnalyzeRlist} now annotated with all the analysis described above
    \item Generate one folder per comparison of conditions containing the isoform switch analysis plot of all genes with significant isoforms switches
    \item Saves 3 plots summarizing the overall consequences of all the isoform switchces.
}
}

\references{
Vitting-Seerup et al. The Landscape of Isoform Switches in Human Cancers. Mol. Cancer Res. (2017).
}
\author{
Kristoffer Vitting-Seerup
}

\seealso{
\code{\link{analyzeCPAT}}\cr
\code{\link{analyzeCPC2}}\cr
\code{\link{analyzeIUPred2A}}\cr
\code{\link{analyzeNetSurfP2}}\cr
\code{\link{analyzePFAM}}\cr
\code{\link{analyzeSignalP}}\cr
\code{\link{analyzeAlternativeSplicing}}\cr
\code{\link{extractSwitchSummary}}\cr
\code{\link{analyzeSwitchConsequences}}\cr
\code{\link{switchPlotTopSwitches}}
}

\examples{
### Please note
# The way of importing files in the following example with
# "system.file('pathToFile', package="IsoformSwitchAnalyzeR") is
# specialized way of accessing the example data in the IsoformSwitchAnalyzeR package
# and not smoothing you need to do - just supply the string e.g.
# "/path/to/externalAnalysis/toolResult.txt" pointing to the result file.

### Load example data
data("exampleSwitchListIntermediary")

### Subset for quick runtime
exampleSwitchListIntermediary <- subsetSwitchAnalyzeRlist(
    exampleSwitchListIntermediary,
    abs(exampleSwitchListIntermediary$isoformFeatures$dIF) > 0.4
)

### Run part 2
exampleSwitchListAnalyzed <- isoformSwitchAnalysisPart2(
    switchAnalyzeRlist        = exampleSwitchListIntermediary,
    pathToCPC2resultFile      = system.file("extdata/cpc2_result.txt", package = "IsoformSwitchAnalyzeR"),
    pathToPFAMresultFile      = system.file("extdata/pfam_results.txt", package = "IsoformSwitchAnalyzeR"),
    pathToIUPred2AresultFile  = system.file("extdata/iupred2a_result.txt.gz", package = "IsoformSwitchAnalyzeR"),
    pathToSignalPresultFile   = system.file("extdata/signalP_results.txt", package = "IsoformSwitchAnalyzeR"),
    codingCutoff              = 0.5,   # since we are using CPC2
    removeNoncodinORFs        = TRUE,  # Because ORF was predicted de novo
    outputPlots               = FALSE  # keeps the function from outputting the plots from this example code
)
}
