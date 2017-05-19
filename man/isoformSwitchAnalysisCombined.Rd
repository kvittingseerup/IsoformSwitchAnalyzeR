\name{isoformSwitchAnalysisCombined}
\alias{isoformSwitchAnalysisCombined}
\title{
Isoform Switch Analysis Workflow: Extract, Annoate and Visualize all Significant Isoform Swtiches
}
\description{
This high-level function takes a CuffSet object or a pre-existing switchAnalyzeRlist as input. If the input is a CuffSet a switchFindeRList is created and else the function uses the provieded switchAnalyzeRlist.

Then isoform switches are identified, annotated with ORF and intron retion. Then functional consequences are identified and isoform switch analysis plots are generated for the top \code{n} isoform switches. Lastly a plot summarizing the global effect of isoform switches with functional consequences is generated. If external analysis of protein domians (Pfam), coding potential (CPAT) or signal peptides (SignalP) should be incorporated please use the combination of \code{isoformSwitchAnalysisPart1} and \code{isoformSwitchAnalysisPart2} instead.
}
\usage{
isoformSwitchAnalysisCombined(
    input,
    alpha=0.05,
    dIFcutoff = 0.1,
    switchTestMethod='DRIMSeq',
    calibratePvalues=TRUE,
    n=NA,
    pathToOutput=getwd(),
    overwriteORF=FALSE,
    outputSequences=FALSE,
    genomeObject,
    orfMethod='longest',
    cds=NULL,
    consequencesToAnalyze=c('intron_retention','ORF_seq_similarity','NMD_status'),
    fileType='pdf',
    asFractionTotal=FALSE,
    outputPlots=TRUE,
    quiet=FALSE
)
}

\arguments{
\item{input}{
Either a \code{CuffSet} created by readCufflins() or a \code{switchAnalyzeRlist} object.
}
  \item{alpha}{
The cutoff which the (callibrated) fdr correct p-values must be smaller than for calling significant switches. Defualit is 0.05.
}
\item{dIFcutoff}{
The cutoff which the changes in (absolute) isoform usage must be larger than before an isoform is considered switching. This cutoff can remove cases where isoforms with (very) low dIF values are deemed signicant and thereby included in the downstream analysis. This cutoff is analogours to having a cutoff on log2 fold change in a normal differential expression analysis of genes to ensure the genes have a certain effect size. Default is 0.1 (10\%).
}
\item{switchTestMethod}{ A sting indicating which statistical method should be used for testing differential isoform usage. The following options are avilable:
        \itemize{
            \item {\code{'DRIMSeq'} : Uses the DRIMSeq package to perform the statiscal test. See \link{isoformSwitchTestDRIMSeq}. Default}
            \item {\code{'IsoformSwitchAnalyzeR'} : Uses the statistical test implemented in IsoformSwitchAnalyzeR. See \link{isoformSwitchTest}.}
            \item {\code{'none'} : No statistical test is performed. Should only be used if a test have already been performed and should not be overwritten (e.g when importing cuffdiff data).}
        }
    }

    \item{calibratePvalues}{Only applicable if \code{switchTestMethod='IsoformSwitchAnalyzeR'}.A logic indicating whether the p-values should be calibrated accordin to the methods in Ferguson et al (2014), this will only be done if the p-value distribution is very skewed (if the estimated sigma^2 values is < 0.9) as suggested by the authors. See details for more information. Defualt is TRUE.}

\item{n}{
The number of top genes (after filtering and sorted according to \code{sortByQvals}) that should be saved to each subfolder indicated by \code{splitConditions, splitFunctionalConsequences}. Use NA to create all. Default is NA (all).
}
  \item{pathToOutput}{
A path to the folder in which the plots should be made. Default is working directory ( getwd() ).
}

\item{overwriteORF}{
A logical indicating whether to overwrite the ORF analysis already stored in the supplied switchAnalyzeRlist. Default is FALSE.
}
  \item{outputSequences}{
A logic indicating whether transcript nucleotide and amino acid sequences should be outputtet to outputDestination. Default is TRUE.
}
  \item{genomeObject}{
A \code{BSgenome} object (for example Hsapiens for Homo sapiens).
}
  \item{orfMethod}{
A string indicating which of the 3 ORF identification methods should be used. The methods are:
  \itemize{
    \item {\code{longest} : Identifies the longest ORF in the transcript. This approach is similar to what the CPAT tool uses in it's analysis of coding potential}
    \item {\code{longestAnnotated} : Identifies the longest ORF downstream of an annotated translation start site (supplied via the \code{cds} argument)}
    \item {\code{mostUpstreamAnnoated} : Identifies the ORF downstream of the most upstream overlapping annotated translation start site (supplied via the \code{cds} argument)}
  }
  Default is \code{longest}.
}
  \item{cds}{
A CDSSet object containing annotated coding regions, see ?CDSSet and ?getCDS for more information. Only nessesary if \'orfType\' arguments is \'longestAnnotated\' or \'mostUpstreamAnnoated\'.
}
  \item{consequencesToAnalyze}{
A vector of strings indicating what type of functional consequences to analyze. Note there is bound to be some differences between transcripts (else there would be identical). See details in \link{analyzeSwitchConsequences} for full list of usable strings and their meaning. Default is c('intron_retention','ORF_seq_similarity','NMD_status') (corresponding to analyze: intron retention, ORF AA sequence similarity and NMD status).
}
  \item{fileType}{
A string indicating which file type is generated. Available options are \'pdf\' and \'png\'. Default is pdf.
}
   \item{asFractionTotal}{
A logic indicating whether the number of consequences should be calculated as numbers (if FALSE) or as a fraction of the total number of switches in the plot summarizing general consequences of all the isoform swithces. Default is FALSE.
}
\item{outputPlots}{
A logic indicating whether all isoform switches as well as the summary of functional consequences should be outputted in the directory specified by \code{pathToOutput} argument. Default is TRUE.
}
\item{quiet}{ A logic indicating whether to avoid printing progress messages (incl. progress bar). Default is FALSE}

}

\details{
This function performs the full Isoform Analysis Workflow by
\enumerate{
    \item If necessary import the data (see \link{importCufflinksCummeRbund})
    \item Remove non-expressed isoforms and single-isoform genes (see \link{preFilter})
    \item predict swithces (only if switches is not already annotated, see \link{isoformSwitchTest})
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
Manuscript in preperation.
}
\author{
Kristoffer Vitting-Seerup
}
\seealso{
\code{\link{isoformSwitchAnalysisPart1}}\cr
\code{\link{isoformSwitchAnalysisPart2}}\cr
\code{\link{preFilter}}\cr
\code{\link{isoformSwitchTest}}\cr
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
    input=exampleSwitchList,
    genomeObject = Hsapiens,
    dIFcutoff = 0.4,         # Set high for short runtime in example data
    outputSequences = FALSE, # keeps the function from outputting the fasta files from this example
    outputPlots = FALSE      # keeps the function from outputting the Isoform Switch AnalyzeR Plots from this example
)

exampleSwitchListAnalyzed
}
