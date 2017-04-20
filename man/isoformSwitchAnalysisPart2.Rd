\name{isoformSwitchAnalysisPart2}
\alias{isoformSwitchAnalysisPart2}
\title{
Isoform Switch Analysis Workflow Part 2: Plot All Isoform Switches and Their Annotation
}
\description{
This high-level function adds the results of the extrenal sequence analysis supplied (if any), then proceeds to annotate intron ration. Then functional consequences of the isoform switches are identified and isoform switch analysis plots are created for the top \code{n} isoform switches. Lastly a plot summarizing the functional consequences is created. This function is meant to be used after \link{isoformSwitchAnalysisPart1} have been used.
}
\usage{
isoformSwitchAnalysisPart2(
    switchAnalyzeRlist,
    alpha=0.05,
    dIFcutoff = 0.1,
    n=NA,
    codingCutoff,
    removeNoncodinORFs,
    pathToCPATresultFile=NULL,
    pathToPFAMresultFile=NULL,
    pathToSignalPresultFile=NULL,
    consequencesToAnalyze=c('intron_retention','coding_potential','ORF_seq_similarity','NMD_status','domains_identified','signal_peptide_identified'),
    pathToOutput=getwd(),
    fileType='pdf',
    asFractionTotal=FALSE,
    outputPlots=TRUE,
    quiet=FALSE
)
}

\arguments{
  \item{switchAnalyzeRlist}{
A \code{switchAnalyzeRlist} object as produced by \link{createSwitchAnalyzeRlist}
}
  \item{alpha}{
The cutoff which the (callibrated) fdr correct p-values must be smaller than for calling significant switches. Defualit is 0.05.
}
\item{dIFcutoff}{
The cutoff which the changes in (absolute) isoform usage must be larger than before an isoform is considered switching. This cutoff can remove cases where isoforms with (very) low dIF values are deemed signicant and thereby included in the downstream analysis. This cutoff is analogours to having a cutoff on log2 fold change in a normal differential expression analysis of genes to ensure the genes have a certain effect size. Default is 0.1 (10\%).
}
\item{n}{
The number of top genes (after filtering and sorted according to \code{sortByQvals}) that should be saved to each subfolder indicated by \code{splitConditions, splitFunctionalConsequences}. Use NA to create all. Default is NA (all).
}
  \item{codingCutoff}{
Numeric indicating the cutoff used by CPAT for distinguishing between coding and non-coding transcripts. The cutoff is dependent on species analyzed. Our analysis suggest that the optimal cutoff for overlapping coding and noncoding isoforms are 0.725 for human and 0.721 for mouse - HOWEVER the suggested cutoffs from the CPAT article (see references) derived by comparing known genes to random non-coding regions of the genome is 0.364 for human and 0.44 for mouse.
}

 \item{removeNoncodinORFs}{
A logic indicating wether to remove ORF information from the isoforms which the CPAT analysis classifies as non-coding. This can be particular useful if the isoform (and ORF) was predicted de-novo but is not recommended if ORFs was imported from a GTF file. This will affect all downstream analysis and plots as both analysis of domains and signal peptides requires that ORFs are annotated (e.g. analyzeSwitchConsequences will not consider the domains (potentially) found by Pfam if the ORF have been removed).
}

  \item{pathToCPATresultFile}{
Path to the CPAT result file. If the webserver is used please dowload the tab-delimited file from the bottom of the result page and give that as input, else simply supply the result file.
}
  \item{pathToPFAMresultFile}{
  Path to the PFAM result file. If the webserver is used you need to copy paste the result part of the mail you get into a empty plain text document (notepad, sublimetext TextEdit or similar (aka not word)) and save that to a file. This file is then used as input to the function.
}
  \item{pathToSignalPresultFile}{
  Path to the summary SignalP result file. If using the web-server the results should be copy pasted into a empty plain text document (notepad, sublimetext TextEdit or similar (aka not word)) and save that to a file. This file is then used as input to the function.
}
  \item{consequencesToAnalyze}{
    A vector of strings indicating what type of functional consequences to analyze. Do note that there is bound to be some differences beteen transcripts (else there would be identical). See details in \link{analyzeSwitchConsequences} for full list of usable strings and their meaning. Default is c('intron_retention','coding_potential','ORF_seq_similarity','NMD_status','domains_identified','signal_peptide_identified') (corresponding to analyze: intron retention, CPAT result, ORF AA sequence similarity, NMD status, PFAM domains annotated and signal peptides annotated by Pfam).
}
  \item{pathToOutput}{
A path to the folder in which the plots should be made. Default is working directory ( getwd() ).
}
  \item{fileType}{
A string indicating which file type is generated. Available options are \'pdf\' and \'png\'. Default is pdf.
}
  \item{asFractionTotal}{
A logic indicating whether the number of consequences should be calculated as numbers (if FALSE) or as a fraction of the total number of switches in the plot summarizing general consequences of all the isoform swithces. Default is FALSE.
}
\item{outputPlots}{
A logic indicating whether all isoform switches as well as the summary of functional consequences should be ouputted in the directory specified by \code{pathToOutput} argument. Default is TRUE.
}
\item{quiet}{ A logic indicating whether to avoid printing progress messages (incl. progress bar). Default is FALSE}

}

\details{
This function performs the second part of a Isoform Analysis Workflow by:
\enumerate{
    \item Adding external sequence analysis (see \link{analyzeCPAT}, \link{analyzePFAM} and \link{analyzeSignalP})
    \item Predict functional consequences of switching (see \link{analyzeSwitchConsequences})
    \item Output Isoform Switch Consequence plots for all genes where there is a significant isoform switch (see \link{switchPlot})
    \item Output a visualization of general consequences of isoform switches.
}
}

\value{
This function ouputs
\enumerate{
    \item The supplied \code{switchAnalyzeRlist} now annotated with all the analysis described above
    \item One folder per comparison of condition containing the isoform switch analysis plot of all genes with significant isoforms switches
    \item A plot summarizing the overall consequences off all the isoform switchces.
}
}

\references{
Manuscript in preperation.
}
\author{
Kristoffer Vitting-Seerup
}

\seealso{
\code{\link{analyzeCPAT}}
\code{\link{analyzePFAM}}
\code{\link{analyzeSignalP}}
\code{\link{extractSwitchSummary}}
\code{\link{analyzeSwitchConsequences}}
\code{\link{switchPlotTopSwitches}}
}

\examples{
data("exampleSwitchListIntermediary")
exampleSwitchListAnalyzed <- isoformSwitchAnalysisPart2(
    switchAnalyzeRlist       = exampleSwitchListIntermediary,
    dIFcutoff                = 0.4,   # Set high for short runtime in example data
    pathToCPATresultFile     = system.file("extdata/cpat_results.txt", package = "IsoformSwitchAnalyzeR"),
    pathToPFAMresultFile     = system.file("extdata/pfam_results.txt", package = "IsoformSwitchAnalyzeR"),
    pathToSignalPresultFile  = system.file("extdata/signalP_results.txt", package = "IsoformSwitchAnalyzeR"),
    codingCutoff             = 0.364,
    removeNoncodinORFs       = TRUE,  # Because ORF was predicted de novo
    outputPlots              = FALSE  # keeps the function from outputting the plots from this example
)
}
