\name{importPairedGSEA}
\alias{importPairedGSEA}
\title{
Import PairedGSEA Results and GTF File into IsoformSwitchAnalyzeR
}
\description{
Function for importing paired differential expression and splicing results from pairedGSEA, along with a (gziped or unpacked) GTF file into R as a \code{switchAnalyzeRlist}. This function filters the input based on specified significance thresholds and ensures compatibility with IsoformSwitchAnalyzeR.
}
\usage{
importPairedGSEA(
    pairedGSEAResults,
    pathToGTF,
    expressionCutoff = 0.05,
    splicingCutoff = 0.05,
    minLogFCExpression = 0,
    minLogFCSplicing = 0,
    isoformNtFasta = NULL,
    extractAaSeq = FALSE,
    addAnnotatedORFs = TRUE,
    onlyConsiderFullORF = FALSE,
    removeNonConvensionalChr = FALSE,
    ignoreAfterBar = TRUE,
    ignoreAfterSpace = TRUE,
    ignoreAfterPeriod = FALSE,
    removeTECgenes = TRUE,
    PTCDistance = 50,
    removeFusionTranscripts = TRUE,
    removeUnstrandedTranscripts = TRUE,
    quiet = FALSE
)
}
\arguments{
\item{pairedGSEAResults}{
A \code{data.frame} or \code{DFrame} containing the paired differential expression and splicing results from pairedGSEA. It must contain the following columns:
\itemize{
    \item{\code{gene}: Gene identifiers.}
    \item{\code{lfc_expression}: Log fold change for expression.}
    \item{\code{padj_expression}: Adjusted p-value for expression.}
    \item{\code{lfc_splicing}: Log fold change for splicing.}
    \item{\code{padj_splicing}: Adjusted p-value for splicing.}
}
}

\item{pathToGTF}{
Can either be:
\itemize{
    \item{\code{1}: A string indicating the full path to the (gziped or unpacked) GTF file which have been quantified. If supplied the exon structure and isoform annotation will be obtained from the GTF file. An example could be "myAnnotation/myGenome/isoformsQuantified.gtf")}
    \item{\code{2}: A string indicating the full path to the (gziped or unpacked) RefSeq GFF file which have been quantified. If supplied the exon structure and isoform annotation will be obtained from the GFF file. Please note only GFF files from RefSeq downloaded from \url{ftp://ftp.ncbi.nlm.nih.gov/genomes/} are supported (see database FAQ in vignette for more info). An example could be "RefSeq/isoformsQuantified.gff")}
}
}

\item{expressionCutoff}{
A numeric specifying the maximum adjusted p-value for filtering significant differential expression. Default is 0.05.
}

\item{splicingCutoff}{
A numeric specifying the maximum adjusted p-value for filtering significant splicing changes. Default is 0.05.
}

\item{minLogFCExpression}{
A numeric specifying the minimum absolute log fold change for filtering significant differential expression. Default is 0.
}

\item{minLogFCSplicing}{
A numeric specifying the minimum absolute log fold change for filtering significant splicing changes. Default is 0.
}

\item{isoformNtFasta}{
A (vector of) text string(s) providing the path(s) to the a fasta file containing the nucleotide sequence of all isoforms quantified. This is useful for: 
\itemize{
    \item{1) People working with non-model organisms where extracting the sequence from a BSgenome might require extra work.}
    \item{2) Workflow speed-up for people who already have the fasta file (which most people running Salmon, Kallisto, or RSEM for the quantification have as that is used to build the index).}
}
The file will automatically be subsetted to the isoforms found in the GTF file, so additional sequences (such as decoys) do not need to be manually removed. Note this differs from a fasta file with sequences of the entire genome.
}

\item{extractAaSeq}{
A logic indicating whether the nucleotide sequence imported via \code{isoformNtFasta} should be translated to amino acid sequence and stored in the \code{switchAnalyzeRlist}. Requires ORFs to be imported, see \code{addAnnotatedORFs}. Default is FALSE.
}

\item{addAnnotatedORFs}{
A logic indicating whether the ORFs from the GTF should be added to the \code{switchAnalyzeRlist}. ORFs are defined as regions annotated as 'CDS' in the 'type' column (column 3). Default is TRUE.
}

\item{onlyConsiderFullORF}{
 A logic indicating whether the ORFs added should only be added if they are fully annotated. Here fully annotated is defined as those that both have a annotated 'start_codon' and 'stop_codon' in the 'type' column (column 3). This argument is only considered if onlyConsiderFullORF=TRUE. Default is FALSE.
}

\item{removeNonConvensionalChr}{
 A logic indicating whether non-conventional chromosomes, here defined as chromosome names containing either a '_' or a period ('.'). These regions are typically used to annotate regions that cannot be associated to a specific region (such as the human 'chr1_gl000191_random') or regions quite different due to different haplotypes (e.g. the 'chr6_cox_hap2'). Default is FALSE.
}

\item{ignoreAfterBar}{
A logic indicating whether to subset the isoform ids by ignoring everything after the first bar ("|"). Useful for analysis of GENCODE files. Default is TRUE.
}

\item{ignoreAfterSpace}{
A logic indicating whether to subset the isoform ids by ignoring everything after the first space (" "). Useful for analysis of gffutils generated GTF files. Default is TRUE.
}

\item{ignoreAfterPeriod}{
A logic indicating whether to subset the gene/isoform is by ignoring everything after the first period ("."). Should be used with care. Default is FALSE.
}

\item{removeTECgenes}{
A logic indicating whether to remove genes marked as "To be Experimentally Confirmed" (if annotation is available). The default is TRUE aka to remove them which is in line with Gencode recommendations (TEC are not in Gencode annotations). For more info about TEC see \url{https://www.gencodegenes.org/pages/biotypes.html}.
}

\item{PTCDistance}{
Only considered if \code{addAnnotatedORFs=TRUE}. A numeric giving the premature termination codon-distance: The minimum distance from the annotated STOP to the final exon-exon junction, for a transcript to be marked as NMD-sensitive. Default is 50
}

\item{removeFusionTranscripts}{
A logic indicating whether to remove genes with cross-chromosome fusion transcripts as IsoformSwitchAnalyzeR cannot handle them.
}

\item{removeUnstrandedTranscripts}{
A logic indicating whether to remove non-stranded isoforms as the IsoformSwitchAnalyzeR workflow cannot handle them.
}

\item{quiet}{
A logic indicating whether to suppress progress messages. Default is FALSE.
}
}
\details{
This function is designed to import a GTF file (leveraging the \code{importGTF} function) and integrate differential expression and splicing results from pairedGSEA. By combining these datasets, the function filters the transcripts in the GTF file to retain only those associated with significant genes, as defined by user-specified thresholds for p-values and log fold changes.

The primary goal is to streamline downstream analyses by focusing on a reduced set of transcripts relevant to the genes identified as significant in pairedGSEA results.
}
\value{
A \code{switchAnalyzeRlist} containing filtered and annotated gene and isoform information. See ?switchAnalyzeRlist for more details.

If no genes match after filtering, an empty \code{switchAnalyzeRlist} is returned with a warning.
}

\references{
Vitting-Seerup et al. The Landscape of Isoform Switches in Human Cancers. Mol. Cancer Res. (2017).
Dam, S.H., Olsen, L.R. & Vitting-Seerup, K. Expression and splicing mediate distinct biological signals. BMC Biol 21, 220 (2023).
}
\author{
Kristoffer Vitting-Seerup, Chunxu Han
}
\seealso{
\code{\link{createSwitchAnalyzeRlist}}\cr
\code{\link{importGTF}}
}