\name{importPredefinedSwitches}
\alias{importPredefinedSwitches}
\title{
Import Isoform Switches Defined Outside IsoformSwitchAnalyzeR
}
\description{
Function for constructing a \code{switchAnalyzeRlist} from a set of isoform switches that were identified by some external approach (for example a splicing/isoform QTL analysis) rather than by IsoformSwitchAnalyzeR's own statistical testing (\link{isoformSwitchTestDEXSeq} / \link{isoformSwitchTestSatuRn}). This makes it possible to use the rest of the IsoformSwitchAnalyzeR workflow - ORF analysis, external sequence analysis, prediction of functional consequences and transcript-structure visualization (\link{switchPlotTranscript}) - on switches you already know about.

Since these switches were not derived from replicate expression data within IsoformSwitchAnalyzeR, functions that rely on such data (e.g. \link{preFilter}, \link{isoformSwitchTestDEXSeq}, \link{switchPlot}) cannot be used on the resulting \code{switchAnalyzeRlist} - see details.
}
\usage{
importPredefinedSwitches(
    ### Core arguments
    preDefinedSwitches,
    isoformExonAnnoation,
    isoformNtFasta = NULL,
    omicDataType = 'transcriptomics',

    ### Advanced arguments
    ignoreAfterBar = TRUE,
    ignoreAfterSpace = TRUE,
    ignoreAfterPeriod = FALSE,
    removeNonConvensionalChr = FALSE,
    removeTECgenes = TRUE,
    PTCDistance = 50,
    addAnnotatedORFs = TRUE,
    onlyConsiderFullORF = FALSE,
    removeFusionTranscripts = TRUE,
    removeUnstrandedTranscripts = TRUE,
    showProgress = TRUE,
    quiet = FALSE
)
}
\arguments{
\item{preDefinedSwitches}{
A \code{data.frame} with one row per isoform switch and the following required columns:
\itemize{
    \item{\code{isoformsUsedMore}: The isoform_id of the isoform used more (up-regulated) in \code{condition_2} compared to \code{condition_1}.}
    \item{\code{isoformsUsedLess}: The isoform_id of the isoform used less (down-regulated) in \code{condition_2} compared to \code{condition_1}.}
    \item{\code{condition_1}: The name of the first (typically reference/control) condition of the comparison this switch belongs to.}
    \item{\code{condition_2}: The name of the second (typically case/alternative) condition of the comparison this switch belongs to.}
}
Multiple comparisons (different \code{condition_1}/\code{condition_2} pairs) can be supplied in the same call - each unique combination is treated as its own comparison, analogous to the \code{comparisonsToMake} argument used elsewhere in IsoformSwitchAnalyzeR.

An optional \code{switch_q_value} column (the multiple-testing corrected significance of the switch, in the [0,1] interval) can also be supplied and will be stored as \code{isoform_switch_q_value} in the resulting \code{switchAnalyzeRlist}. If not supplied, a sentinel value is used instead (see details).
}

\item{isoformExonAnnoation}{
A string indicating the full path to the (gziped or unpacked) GTF/GFF file containing the exon structure and isoform annotation for the isoforms in \code{preDefinedSwitches}. Passed to \link{importGTF} - see that function for details on supported formats.
}

\item{isoformNtFasta}{
A (vector of) text string(s) providing the path(s) to a fasta file with the nucleotide sequence of the isoforms. Passed to \link{importGTF} - see that function for details. Default is NULL. Cannot be combined with \code{omicDataType="proteomics"} - see that argument for why.
}

\item{omicDataType}{
A string, either \code{"transcriptomics"} (default) or \code{"proteomics"}, stating what kind of data the pre-defined switches were derived from.

With \code{"transcriptomics"}, \link{importGTF} is called with \code{exonFeatureType="exon"} (the normal full exon structure), and the resulting \code{switchAnalyzeRlist} has \code{omicDataType="transcriptomics"}.

With \code{"proteomics"}, \link{importGTF} is called with \code{exonFeatureType="CDS"} instead - so the \code{exons} slot (and everything built from it) is coding-sequence-only, non-coding isoforms are dropped, and NMD/PTC status is not determined (see \code{exonFeatureType} in \link{importGTF} for full details) - and the resulting \code{switchAnalyzeRlist} has \code{omicDataType="proteomics"}. Because of this, \code{"proteomics"} cannot be combined with \code{isoformNtFasta} (same restriction as \code{importGTF(exonFeatureType="CDS")}) - extract sequence afterwards via \code{extractSequence(genomeObject=<a BSgenome>)} instead.
}

\item{ignoreAfterBar}{A logic indicating whether to subset the isoform ids (in both \code{preDefinedSwitches} and the annotation) by ignoring everything after the first bar ("|"). Default is TRUE.}

\item{ignoreAfterSpace}{A logic indicating whether to subset the isoform ids by ignoring everything after the first space (" "). Default is TRUE.}

\item{ignoreAfterPeriod}{A logic indicating whether to subset the isoform ids by ignoring everything after the first period ("."). Should be used with care. Default is FALSE.}

\item{removeNonConvensionalChr}{A logic indicating whether to remove non-conventional chromosomes. Passed to \link{importGTF}. Default is FALSE.}

\item{removeTECgenes}{A logic indicating whether to remove genes marked as "To be Experimentally Confirmed". Passed to \link{importGTF}. Default is TRUE.}

\item{PTCDistance}{ Only considered if \code{addAnnotatedORFs=TRUE}. The premature termination codon-distance passed to \link{importGTF}. Default is 50.}

\item{addAnnotatedORFs}{A logic indicating whether the ORF annotated in the GTF should be added to the \code{switchAnalyzeRlist}. Passed to \link{importGTF}. Default is TRUE.}

\item{onlyConsiderFullORF}{A logic indicating whether only fully annotated ORFs should be added. Passed to \link{importGTF}. Default is FALSE.}

\item{removeFusionTranscripts}{A logic indicating whether to remove genes with cross-chromosome fusion transcripts. Passed to \link{importGTF}. Default is TRUE.}

\item{removeUnstrandedTranscripts}{A logic indicating whether to remove non-stranded isoforms. Passed to \link{importGTF}. Default is TRUE.}

\item{showProgress}{A logic indicating whether to show a progress bar. Currently not used but kept for consistency with other import functions. Default is TRUE.}

\item{quiet}{ A logic indicating whether to avoid printing progress messages. Default is FALSE.}
}

\details{
This function uses \link{importGTF} to build the isoform/gene/exon annotation and then overwrites the statistical fields normally estimated by IsoformSwitchAnalyzeR with sentinel values that mark the switches as pre-defined:
\itemize{
    \item{\code{dIF} is set to \code{Inf} for isoforms in \code{isoformsUsedMore}, \code{-Inf} for isoforms in \code{isoformsUsedLess}, and left as \code{NA} for the gene's other isoforms (which are still retained in the \code{switchAnalyzeRlist} so the full transcript model can be visualized with \link{switchPlotTranscript}).}
    \item{\code{gene_switch_q_value} and (if \code{switch_q_value} was not supplied) \code{isoform_switch_q_value} are set to \code{-1} - a value guaranteed to pass any \code{alpha} significance cutoff.}
}

The resulting \code{switchAnalyzeRlist} has \code{sourceId} set to \code{"preDefinedSwitches"}. Because it does not contain real replicate expression data or a real experimental design, several functions that depend on those (\link{preFilter}, \link{isoformSwitchTestDEXSeq}, \link{isoformSwitchTestSatuRn}, \link{switchPlot}, \link{expressionAnalysisPlot}, \link{extractConsequenceGenomeWide}, \link{extractSplicingGenomeWide}) will refuse to run on it. Use \link{switchPlotTranscript} instead of \link{switchPlot} for visualization, and \link{isoformSwitchAnalysisPart1}/\link{isoformSwitchAnalysisPart2} (which automatically skip the filtering/testing steps for a pre-defined switchAnalyzeRlist) to run the ORF/annotation/consequence part of the workflow.
}
\value{
A \code{switchAnalyzeRlist} with \code{sourceId} set to \code{"preDefinedSwitches"} and \code{omicDataType} set to whatever was supplied via the \code{omicDataType} argument (\code{"transcriptomics"} or \code{"proteomics"}), containing one \code{isoformFeatures} row per isoform per comparison it is part of. See \link{switchAnalyzeRlist} for more details.
}
\references{
Vitting-Seerup et al. The Landscape of Isoform Switches in Human Cancers. Mol. Cancer Res. (2017).
}
\author{
Kristoffer Vitting-Seerup
}
\seealso{
\code{\link{importGTF}}\cr
\code{\link{createSwitchAnalyzeRlist}}\cr
\code{\link{switchPlotTranscript}}
}
\examples{
predefinedSwitches <- data.frame(
    isoformsUsedMore = 'TCONS_00000234',
    isoformsUsedLess = 'TCONS_00000235',
    condition_1 = 'test1',
    condition_2 = 'test2'
)

aSwitchList <- importPredefinedSwitches(
    preDefinedSwitches   = predefinedSwitches,
    isoformExonAnnoation = system.file("extdata/example.gtf.gz", package="IsoformSwitchAnalyzeR")
)
aSwitchList
}
